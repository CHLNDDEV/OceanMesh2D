function obj = Make_f11(obj,D_filename,file_suffixes,bathyfile)
% obj = Make_f11(obj,D_filename,file_suffixes)
% Input a msh class object get the values of the density over the depth 
% based on D_filename, and computes the depth-averaged value which
% populates the f11 struct in the msh class object.
%
%  Author:      William Pringle                                 
%  Created:     March 19 2018                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(obj.b)
    error('No bathymetry data in grid to calculate the depth-averaged value')
end
rho0 = 1000;
%g = 9.807;

% %% Do some projection for obj.p
% proj = 'Mercator';
% R = 6378206.4; %[m]
% m_proj(proj,'lon',[ min(obj.p(:,1)) max(obj.p(:,1)) ],...
%             'lat',[ min(obj.p(:,2)) max(obj.p(:,2))]) 
% [xx,yy] = m_ll2xy(obj.p(:,1),obj.p(:,2));
% xx = R*xx; yy= R*yy; 
% 
% %% Get the element areas and slopes, connectivity table
% A = polyarea(xx(obj.t(:,1:3))',yy(obj.t(:,1:3))')'; 
% % Get x differences
% a = [ xx(obj.t(:,3)) - xx(obj.t(:,2))
%       xx(obj.t(:,1)) - xx(obj.t(:,3))   
%       xx(obj.t(:,2)) - xx(obj.t(:,1)) ] ;
% a = reshape(a,[],3);
% 
% % Get y differences
% b = [ yy(obj.t(:,2)) - yy(obj.t(:,3))
%       yy(obj.t(:,3)) - yy(obj.t(:,1))   
%       yy(obj.t(:,1)) - yy(obj.t(:,2)) ];    
% b = reshape(b,[],3);
% 
% % Get the vertex to element table
% vtoe = VertToEle(obj.t);
% vtoe(vtoe == 0) = length(obj.t) + 1;
% A(end+1) = 0; 
% An = sum(A(vtoe))';

% nearest neighbour extrapolation of T and S (0 for no extrapolation)?
FillNaN = 0;

[direc , title , ext ] = fileparts( D_filename);

if isempty(title)
    D_filename = dir(direc);
    D_filename(1:2) = [];
    idx = zeros(length(D_filename),1);
    for jj = 1:size(file_suffixes,1)
        for ii = 1:length(D_filename)
            if idx(ii) == 0
                idx(ii) = contains(D_filename(ii).name,file_suffixes(jj,:));
            end
        end
    end
    D_filename = D_filename(find(idx));
else
   name = D_filename; D_filename = [];
   D_filename.folder = direc;
   D_filename.name = [title ext];
end

Bx = NaN(length(obj.b),size(D_filename,1));
By = NaN(length(obj.b),size(D_filename,1));
tic
for ii = 1:size(D_filename,1)
    Dn = [D_filename(ii).folder '/' D_filename(ii).name];
    [~ , ~ , ext ] = fileparts( Dn);
    %% Read the grid and density data 
    if strcmp(ext,'.nc')
        try
            title = ncreadatt(Dn,'/','title'); % Reading some 
            lon = ncread(Dn,'lon');
            lat = ncread(Dn,'lat');
            z   = ncread(Dn,'depth');
            sigma_t = ncread(Dn,'I_an');
            time_all = 0;
        catch
            title = ncreadatt(Dn,'/','History'); % Reading some 
            [time_all(ii),lon,lat,z,~,~,dpx,dpy] = Make_Gridded_rho(Dn,...
         {'time','depth','lon','lat','salinity','water_temp','surf_el'}); %,...
                                                                %bathyfile);
            disp(['Read time ' datestr(time_all(ii))])
            %if ii == 1
            %    sigma_t = zeros([size(rho) size(D_filename,1)]);
            %end
            %sigma_t(:,:,:,ii) = rho - 1000; clear rho 
        end
    elseif strcmp(ext,'.mat')
        load(Dn);
        %sigma_t = rho - 1000; clear rho
    else
        error('does not understand file extension')
    end
%     % zeta on the finite-element mesh
%     [Lon,Lat] = ndgrid(lon,lat);
%     F = griddedInterpolant(Lon,Lat,zeta,'linear','none');
%     zeta = F(obj.p);
%     % Get rho on the finite-element mesh
%     [~,~,~,rho] = Gridded_to_Mesh_SeaBed_DepthAveraged(...
%             obj.p(:,1),obj.p(:,2),obj.b,z,rho,lon,lat,FillNaN);  
%         
%     % Compute bx, by for each element (loop over to reduce memory issues)
%     bxe = zeros(length(obj.t),length(z));
%     bye = zeros(length(obj.t),length(z));
%     for k = 1:length(z)
%         if k == 1
%             BCPe = zeta.*rho{k};
%         else
%             BCPe = rho{k} - rho0;
%         end
%         BCPe = BCPe(obj.t);
%         bxe(:,k) = 0.5 * sum(BCPe.*b,2); 
%         bye(:,k) = 0.5 * sum(BCPe.*a,2);
%     end
%     bxe(end+1,:) = 0; bye(end+1,:) = 0;
%     bx = zeros(length(obj.b),length(z));
%     by = zeros(length(obj.b),length(z));
%     for k = length(z):-1:1
%         if k > 1
%             bxe(:,k) = trapz(z(1:k),bxe(:,1:k),2);
%             bye(:,k) = trapz(z(1:k),bye(:,1:k),2);
%         end
%         temp = sum(reshape(bxe(vtoe,k),size(vtoe,1),[]),'omitnan')';
%         bx(:,k) = temp./An;
%         by(:,k) = temp./An;
%     end

    % Interpolate the gradient at each depth
    dpx(isnan(dpy)) = NaN; dpy(isnan(dpx)) = NaN;
    [lonN,latN] = ndgrid(lon,lat);
    bx = NaN(length(obj.p),length(z));
    by = NaN(length(obj.p),length(z));
    for kk = 1:length(z)
        if kk == 1
            Fx = griddedInterpolant(lonN,latN,dpx(:,:,kk),'linear','none');
            Fy = griddedInterpolant(lonN,latN,dpy(:,:,kk),'linear','none');
        else
            Fx.Values = dpx(:,:,kk);
            Fy.Values = dpy(:,:,kk);
        end
        bb = obj.b >= z(kk);
        pt = obj.p(bb,:); 
        bxt = Fx(obj.p(bb,:));
        byt = Fy(obj.p(bb,:));
        if ~isempty(find(isnan(bxt), 1)) 
            idn = find(isnan(bxt));
            idg = find(~isnan(bxt));
            idx = knnsearch(pt(idg,:),pt(idn,:));
            bxt(idn) = bxt(idg(idx));
        end
        if ~isempty(find(isnan(byt), 1)) 
            idn1 = find(isnan(byt));
            idg = find(~isnan(byt));
            if sum(idn - idn1) ~= 0
                idx = knnsearch(pt(idg,:),pt(idn1,:));
            end
            byt(idn1) = byt(idg(idx));
        end
        bx(bb,kk) = bxt;
        by(bb,kk) = byt;
    end
    % Now do the integration over the depth taking into account actual
    % depth on the computational mesh
    for kk = length(z):-1:1
        bb = isnan(Bx(:,ii)) & obj.b > z(kk) & ...
            ~isnan(bx(:,kk)) & ~isnan(by(:,kk));
        Bx(bb,ii) = trapz(z(1:kk),bx(bb,1:kk),2) + ...
                          bx(bb,kk).*(obj.b(bb) - z(kk));
        By(bb,ii) = trapz(z(1:kk),by(bb,1:kk),2) + ...
                          by(bb,kk).*(obj.b(bb) - z(kk)); 
    end
    % Get the depth averaged
    Bx(:,ii) = Bx(:,ii)./obj.b/rho0;
    By(:,ii) = By(:,ii)./obj.b/rho0;

%     figure;
%     fastscatter(obj.p(:,1),obj.p(:,2),9.807*hypot(Bx,By))
%     colormap(cmocean('speed'))
%     caxis([0 1e-4])
%     colorbar;
end
toc
disp(time_all)
DTIMINC = round(seconds(median(diff(time_all))));

%
% tic 
% Bx = zeros(length(obj.b),size(sigma_t,4));
% By = zeros(length(obj.b),size(sigma_t,4));
% for t = 1:size(sigma_t,4)
%     disp(['Computing ' datestr(time_all(t))])
%     %% Calculate the depth-averaged density
%     [~,Sigma_tm,~,Sigma_3D] = Gridded_to_Mesh_SeaBed_DepthAveraged(...
%             obj.p(:,1),obj.p(:,2),obj.b,z,sigma_t(:,:,:,t),lon,lat,FillNaN);                   
% 
%     %% Calculate the integrated barcolinic pressure's
%     Sigma_3D = reshape(cell2mat(Sigma_3D),[],length(Sigma_3D));
%     % Make sure derivative with "ground" returns NaN;
%     BCP = NaN(length(obj.b),length(z));
%     for k = 1:length(z)
%         bb = obj.b >= z(k);
%         BCP(bb,k) = trapz(z(1:k),Sigma_3D(bb,1:k),2);
%     end
%     BCP = BCP/rho0;
% 
%     %% Calculate the bx/by's for each element
%     % Compute bx, by for each element (loop over to reduce memory issues)
%     bxe = zeros(length(obj.t),length(z));
%     bye = zeros(length(obj.t),length(z));
%     for k = 1:length(z)
%         BCPe = BCP(:,k);
%         BCPe = BCPe(obj.t);
%         bxe(:,k) = 0.5 * sum(BCPe.*b,2); 
%         bye(:,k) = 0.5 * sum(BCPe.*a,2);
%     end
%     %% Sum the contributions of each element connected to a node
%     % Add in a zero ghost element and refer to the zero indices of vtoe to this
%     bxe(end+1,:) = 0; bye(end+1,:) = 0;
% 
%     % Sum up all element contributions to each node and take the integral
%     bx = zeros(length(obj.b),length(z));
%     by = zeros(length(obj.b),length(z));
%     for k = 1:length(z)
%         bb = obj.b >= z(k);
%         temp = sum(reshape(bxe(vtoe,k),size(vtoe,1),[]),'omitnan')';
%         bx(bb,k) = temp(bb)./An(bb);
%         temp = sum(reshape(bye(vtoe,k),size(vtoe,1),[]),'omitnan')';
%         by(bb,k) = temp(bb)./An(bb);
%     end
%     Bx(:,t) = trapz(z,bx,2)./obj.b;
%     By(:,t) = trapz(z,by,2)./obj.b;
%     Bx(obj.b == 0,t) = NaN; By(obj.b == 0,t) = NaN;
% %     % On the structured grid
%     %% Calculate the integrated barcolinic pressure's
%     Sigma_3D = sigma_t(:,:,:,t);
%     % Make sure derivative with "ground" returns NaN;
%     BCP = NaN(size(Sigma_3D));
%     for k = 1:length(z)
%         BCP(:,:,k) = trapz(z(1:k),Sigma_3D(:,:,1:k),3);
%     end
%     BCP = BCP/rho0;
%     
%     bx = zeros(size(Sigma_3D));
%     by = zeros(size(Sigma_3D));
%     for k = 1:length(z)
%         [bx(:,:,k),by(:,:,k)] = gradient(BCP(:,:,k),(lon(2)-lon(1))*111e3,(lat(2)-lat(1))*111e3);
%     end
%     bxv = reshape(bx,[],length(z));
%     byv = reshape(by,[],length(z));
%     Bx_struc = NaN(size(bxv,1),1);
%     By_struc = NaN(size(bxv,1),1);
%     for k = 1:length(z)
%         bb = ~isnan(bxv(:,k));
%         Bx_struc(bb) = trapz(z(1:k),bxv(bb,1:k),2)./z(k);
%         bb = ~isnan(byv(:,k));
%         By_struc(bb) = trapz(z(1:k),byv(bb,1:k),2)./z(k);
%     end
%     Bx_struc(Bx_struc == 0) = NaN;
%     By_struc(By_struc == 0) = NaN;
%     Bx_struc = reshape(Bx_struc,size(Sigma_3D,1),[])';
%     By_struc = reshape(By_struc,size(Sigma_3D,1),[])';
%     figure;
%     pcolor(hypot(Bx_struc,By_struc))
%     shading interp
%     colormap(cmocean('ice'))
%     caxis([0 2e-6])
%     colorbar;
%end
%toc
% Calculate the depth-averaged
% rhoxe = sum(Sigma_tm(obj.t).*b,2)/rho0; 
% rhoye = sum(Sigma_tm(obj.t).*a,2)/rho0;
% rhoxe(end+1,:) = 0; rhoye(end+1,:) = 0;
% rhox = sum(rhoxe(vtoe))'./An;
% rhoy = sum(rhoye(vtoe))'./An;
% Bx2D = 0.5*obj.b.^2.*rhox;
% By2D = 0.5*obj.b.^2.*rhoy;


%% Make into f11 struct
%summary = ncreadatt(D_filename,'/','summary'); % global attributes
obj.f11.DataTitle = title;
%[~,name,~] = fileparts(D_filename);
%obj.f11.DataSubTitle = [summary ': ' name];
obj.f11.DTIMINC = DTIMINC;
obj.f11.TimeVec = time_all;
obj.f11.NumOfNodes = length(obj.b);
obj.f11.Val = cell(size(Bx,2),1); %zeros(3,length(obj.b),size(sigma_t,4));
for t = 1:size(Bx,2)
    idx = find(~isnan(Bx(:,t)) & ~isnan(By(:,t)));
    obj.f11.Val{t} = [idx'; Bx(idx,t)'; By(idx,t)']; % Sigma_tm'; 
end
%EOF
end


