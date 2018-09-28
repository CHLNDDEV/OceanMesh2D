 function efs = smooth_outer(efs)
% This method takes a cell-aray of edge function class instances 
% and smoothes them together so they blend into each other.
% Relax gradient of outer edgefx with inner edgefx using limgradStruct
    for ii = length(efs)-1:-1:1
        hh_m = efs{ii}.F.Values; found = 0;
        for bn = ii+1:length(efs)
            % smooth with all inner boxes (with buffer)
            x = efs{ii}.F.GridVectors{1}; dx = x(2) - x(1);
            inx = x >= min(efs{bn}.F.GridVectors{1}) - dx & ...
                  x <= max(efs{bn}.F.GridVectors{1}) + dx;
            y = efs{ii}.F.GridVectors{2}; dy = y(2) - y(1);
            iny = y >= min(efs{bn}.F.GridVectors{2}) - dy & ...
                  y <= max(efs{bn}.F.GridVectors{2}) + dy;
            if isempty(find(inx,1)) || isempty(find(iny,1)); continue; end
            found(bn) = 1;
            % Get the grid of coarse one inside the fine one
            [x,y] = ndgrid(x(inx),y(iny));
            % Use fine griddedInterpolant to interpolate fine to coarse
            hh_t = efs{bn}.F(x,y);
            hh_m(inx,iny) = hh_t;
        end
        if all(found == 0); continue; end
        disp(['Relaxing the gradient of #' num2str(ii) ' outer edgefx ' ...
              'using #' num2str(find(found)) ' inner edgefxs']);
        hfun = zeros(numel(efs{ii}.F.Values),1);
        [xg,yg] = ndgrid(efs{ii}.F.GridVectors{1},efs{ii}.F.GridVectors{2});
        hh_m = ConvertToPlanarMetres(xg,yg,hh_m) ; 
        nn = 0;
        for ipos = 1 : efs{ii}.nx
            for jpos = 1 : efs{ii}.ny
                nn = nn + 1;
                hfun(nn,1) = hh_m(ipos,jpos);
            end
        end
        [hfun,flag] = limgradStruct(efs{ii}.ny,efs{ii}.h0,hfun,...
            efs{ii}.g,sqrt(length(hfun)));
        if flag == 1
            disp('Gradient relaxing converged!');
        else
            error(['FATAL: Gradient relaxing did not converge, '
                'please check your edge functions']);
        end
        % reshape it back
        nn = 0;
        for ipos = 1 : efs{ii}.nx
            for jpos = 1 : efs{ii}.ny
                nn = nn+1;
                hh_m(ipos,jpos) = hfun(nn);
            end
        end
        hh_m = ConvertToWGS84(yg,hh_m) ; 
        % Save it back into the interpolant
        efs{ii}.F.Values = hh_m;
    end
 end
        
 function [ hh_m ] = ConvertToPlanarMetres( xg, yg, hh_m ) 
    %CONVERTTOPLANARMETRES
    % Given a structured grid with edgelength's defined in hh_m,
    % convert it to planar metres using the Haversine formula assuming
    % a spherical earth of radius 6378.137km.
    % 
    % INPUTS:
    % XG : A grid of points in WGS84 containing the x-locations of these
    %      points.
    % YG : A grid of points in WGS*4 containing the y-locations of these
    %      points.
    % HH_M : A grid the same size as XG and YG containing the edgelength in
    %      WGS84 degrees 
    % OUTPUTS: 
    % HH_M : A grid the same size as XG and YG containing the edglength in
    %        planar metres assuming the edgelength is orientated in the x-direction at
    %        each (XG,YG) point. 

    % Try to remove problem near pole
    ul = 89;
    yg(yg > ul) = ul; yg(yg < -ul) = -ul;

    ptsx = xg(:)'; 
    ptsy = yg(:)'; 

    ptsx2 = xg(:)'+hh_m(:)'; 
    ptsy2 = yg(:)'; 

    blahx = [ptsx; ptsx2]; 
    blahy = [ptsy; ptsy2]; 

    [temp]=m_lldist(blahx(:),blahy(:),2);

    hh_m = reshape(temp(1:2:end),size(xg,1),size(xg,2))*1e3; 

end
 
function [ hh_m ] = ConvertToWGS84( yg, hh_m )
%CONVERTTOWGS84 Convert a grid of points describing edgelength in space in
% planar metres to WGS84 degrees.
% Given a structured grid with edgelength's defined in hh_m,
% convert from planar metres to WGS84 degrees using inverse Haversine
% 
% INPUTS:
% YG : A grid of points in WGS*4 containing the y-locations of these
%      points.
% HH_M : A grid the same size as XG and YG containing the edgelength in
%      planar metres 
% OUTPUTS: 
% HH_M : A grid the same size as XG and YG containing the edglength in
%        WGS84 degrees assuming the edgelength is orientated in the x-direction at
%        each (XG,YG) point. 

% Ensure to remove problem at pole
ul = 89;
yg(yg > ul) = ul; yg(yg < -ul) = -ul;

% Convert to projected coordinates 
Re = 6378.137e3; % <-radius of Earth

% We use a simple inverse Harvesine formula by assuming that the hh_m 
% is applied along a latitude parallel, thus the latitude is constant
% between two points and we only need to solve for the difference in
% longitude (the new hh_m)
el = sin(hh_m./(2*Re));
hh_m = 2*asind(el./cosd(yg));

end


 