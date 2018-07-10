function [Hx,Hy] = Unstruc_Bath_Slope( EToV,xx,yy,B)
% [Hx,Hy] = Unstruc_Bath_Slope( EToV,xx,yy,B)
% Unstruc_Bath_Slope : Gets the bathymetric slopes Hx and Hy at each node 
%                      on an unstructured triangular grid
%
% Inputs   : EToV - NE x 3 array of triangular elements (node indices)
%            xx   - N x 1 vector of x points (Cartesian)
%            yy   - N x 1 vector of y points (Cartesian)
%            B    - N x 1 vector of bathymetric depths
%
% Outputs  : Hx   - N x 1 vector of bathymetric slope in x direction
%            Hy   - N x 1 vector of bathymetric slope in y direction
%
% Requires : VertToEle function (included in 'sub/' folder)
%
% Author   : William Pringle, Dec 12 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic
%% Get the element areas
A = polyarea(xx(EToV(:,1:3))',yy(EToV(:,1:3))')'; 
%% Compute the slopes for each element
% Get x differences
a = [ xx(EToV(:,3)) - xx(EToV(:,2))
      xx(EToV(:,1)) - xx(EToV(:,3))   
      xx(EToV(:,2)) - xx(EToV(:,1)) ] ;
a = reshape(a,[],3);

% Get y differences
b = [ yy(EToV(:,2)) - yy(EToV(:,3))
      yy(EToV(:,3)) - yy(EToV(:,1))   
      yy(EToV(:,1)) - yy(EToV(:,2)) ];    
b = reshape(b,[],3);

% Compute Hx for each element
Hxe = 0.5 * sum(B(EToV).*b,2); 
% Compute Hy for each element
Hye = 0.5 * sum(B(EToV).*a,2);

%% Sum the contributions of each element connected to a node
% Get the vertex to element table
vtoe = VertToEle(EToV);
% Add in a zero ghost element and refer to the zero indices of vtoe to this
A(end+1) = 0; Hxe(end+1) = 0; Hye(end+1) = 0;
vtoe(vtoe == 0) = length(EToV) + 1;
% Sum up all element contributions to each node
Hx = sum(Hxe(vtoe))';
Hy = sum(Hye(vtoe))';
An = sum(A(vtoe))';
% Do the division by the area
Hx = Hx./An;
Hy = Hy./An;
%toc
end