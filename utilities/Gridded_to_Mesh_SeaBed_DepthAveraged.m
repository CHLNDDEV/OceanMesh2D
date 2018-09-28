function [Nb,Nm,Nmw,N_interp] = Gridded_to_Mesh_SeaBed_DepthAveraged(...
                              lon_M,lat_M,B,zcontour,N,lon_N,lat_N,fillNaN)
% Gridded_to_Mesh_SeaBed_DepthAveraged: Interpolate gridded data at specified 
% depth contours, such as buoyancy frequencies or densities, to an unstructured
% grid, then compute that values at the seabed (Nb), the depth-averaged
% value (Nm) and the weighted depth-averaged value (Nmw) 
%
% Input : lon_M    - longitude points of nodes in mesh
%         lat_M    - latitude points of nodes in mesh
%         B        - depths of nodes in mesh
%         zcontour - the depth contours where we have values of N
%         N        - matrix of N  (lon,lat,z)
%         lon_N    - vector of lon 
%         lat_N    - vector of lat
%
% Output : Nb     - Data (e.g. Buoyancy frequency) at seabed
%          Nm     - Depth-averaged value
%          Nmw    - Weighted depth-average value (linearly
%                   decreasing weights from bottom to surface)
%
% Author: William Pringle, CHL, Notre Dame University
% Created: 2017-9-28

if nargin == 7
    fillNaN = 1;
end

% Get the ndgrid of lon_N and lat_N
[Lon,Lat] = ndgrid(lon_N,lat_N);

%% Do the interpolation at each depth contour
% initialisation
Nb = zeros(size(B)); Nm = zeros(size(B)); Nmw = zeros(size(B)); 
% do the interpolation onto the mesh
N_interp = cell(length(zcontour),1);

zcontour = double(zcontour);
%N = inpaintn(N);
for zvalue = 1:length(zcontour)
    Nnow = squeeze(N(:,:,zvalue));
    % Make the interpolant using griddedInterpolant spline (need to inpaint
    % nans first if do this)
    %Nnow = inpaint_nans(Nnow);
    %F = griddedInterpolant(Lon,Lat,Nnow,'spline','none');
    % Make the interpolant using griddedInterpolant linear
    F = griddedInterpolant(Lon,Lat,Nnow,'linear','none');
    
    % Find all nodes less than current depth
    J = find( B >= zcontour(zvalue));
    
    % Make the interp cell and interpolate into it
    N_interp{zvalue}    = NaN(size(B));
    N_interp{zvalue}(J) = F(lon_M(J),lat_M(J));
    
    if fillNaN
        % Nearest neighbour fill for NaN results
        while ~isempty(find(isnan(N_interp{zvalue}(J)), 1))
            N_interp{zvalue} = fillmissing(N_interp{zvalue},'nearest',1);
            if size(N_interp{zvalue},2) > 1
                N_interp{zvalue} = fillmissing(N_interp{zvalue},'nearest',2);
            end
        end
    end
end

if nargout == 4
    return;
end

%% Do the calculation over the depth
zcontour(end+1) = Inf; % add on an end for testing purposes
for zvalue = 1:length(zcontour)-1
    % Test for data above current depth but smaller than next contour
    J = find( B > zcontour(zvalue) &  B <= zcontour(zvalue+1)) ;
    if ~isempty(J)
        % For Nb just set equal to the last available contour value
        dz = B(J) - zcontour(zvalue);
        Nb(J) = N_interp{zvalue}(J);

        % For Nm do the integral sum
        % For depths within the current range
        Nm(J) = Nm(J) + Nb(J).*dz;
        Weight = 0.5*(zcontour(zvalue)+B(J))./B(J);
        Nmw(J) = Nmw(J) + Weight.*Nb(J).*dz;
    end
    % For depths larger than current range
    J = find( B > zcontour(zvalue+1)); 
    DZ = zcontour(zvalue+1)- zcontour(zvalue);
    if ~isempty(J)
        Nm(J) = Nm(J) + 0.5*(N_interp{zvalue}(J)+N_interp{zvalue+1}(J))*DZ; 
        Weight = 0.5*(zcontour(zvalue+1)+zcontour(zvalue))./B(J);
        Nmw(J) = Nmw(J) + 0.5*Weight.*...
                        (N_interp{zvalue}(J)+N_interp{zvalue+1}(J))*DZ; 
    end
end
% Divide our integrated values by the depth, and eliminate NaNs
Nm = Nm./(B-zcontour(1));
Nmw = Nmw./(B-zcontour(1));
Nb(isnan(Nb)) = 0;
Nm(isnan(Nm)) = 0;
Nmw(isnan(Nmw)) = 0;
%EOF
end