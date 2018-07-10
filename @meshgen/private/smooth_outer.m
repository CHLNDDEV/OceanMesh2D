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
        if found == 0; continue; end
        disp(['Relaxing the gradient of #' num2str(ii) ' outer edgefx ' ...
              'using #' num2str(find(found)) ' inner edgefxs']);
        hfun = zeros(numel(efs{ii}.F.Values),1);
        nn = 0;
        for ipos = 1 : efs{ii}.nx
            for jpos = 1 : efs{ii}.ny
                nn = nn + 1;
                hfun(nn,1) = hh_m(ipos,jpos);
            end
        end
        [hfun,flag] = limgradStruct(efs{ii}.ny,efs{ii}.gridspace,hfun,...
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
        % Save it back into the interpolant
        efs{ii}.F.Values = hh_m;
    end
 end
        