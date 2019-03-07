function[h]=lansey(N)
%LANSEY  The Lansey modification of Cynthia Brewer's "Spectral" colormap.
%   __________________________________________________________________
%
%   *|* lansey.png --- The lansey colormap versus parula and jet.  
%   Type 'jhelp lansey' to view this image. *|*
%   __________________________________________________________________
%
%   LANSEY(M) returns an M-by-3 matrix containing the Lansey colormap. 
%
%   LANSEY with no arguments returns a colormap having the same number of
%   colors as the colormap of the current figure.
%
%   This colormap is Jonathan Lansey's modification of the 11-division 
%   version of Cynthia Brewer's "Spectral" colormap.
%
%   LINSPECER itself is available from
%
%       http://www.mathworks.com/matlabcentral/fileexchange/42673.
%
%   This product includes color specifications and designs developed by 
%   Cynthia Brewer (http://colorbrewer.org/).
%
%   LANSEY is a simplified version of LINSPECER by Jonathan Lansey,
%   modified and redistributed in accordance with the copyright policies
%   of LINSPECER and ColorBrewer.org, see LANSEY_COPYRIGHT for details. 
%
%   To make LANSEY your default colormap, add to your startup.m file the 
%   line "set(0,'DefaultFigureColormap',lansey)".
% 
%   'lansey --f' generates the figure shown above.
%
%   Usage: h=lansey(M);
%          colormap lansey
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details
 

if nargin>0
    if strcmpi(N, '--f')
       type makefigs_lansey
       makefigs_lansey;
       return
    end
end


if nargin < 1, N = size(get(gcf,'colormap'),1); end

h = cmap2linspecer(colorm(N));
h = cell2mat(h);

function[vOut] = cmap2linspecer(vIn) % changes the format from a double array to a cell array with the right format
vOut = cell(size(vIn,1),1);
for ii=1:size(vIn,1)
    vOut{ii} = vIn(ii,:);
end

% Jonathan Lansay's notes:
% colorm returns a colormap which is really good for creating informative
% heatmap style figures.
% No particular color stands out and it doesn't do too badly for colorblind people either.
% It works by interpolating the data from the
% 'spectral' setting on http://colorbrewer2.org/ set to 11 colors
% It is modified a little to make the brightest yellow a little less bright.
function[cmap] = colorm(varargin)
n = varargin{1};
%frac=1;
frac=.95; % Slight modification from colorbrewer here to make the yellows in the center just a bit darker
cmapp = [158, 1, 66; 213, 62, 79; 244, 109, 67; 253, 174, 97; 254, 224, 139; 255*frac, 255*frac, 191*frac; 230, 245, 152; 171, 221, 164; 102, 194, 165; 50, 136, 189; 94, 79, 162];
x = linspace(1,n,size(cmapp,1));
xi = 1:n;
cmap = zeros(n,3);
for ii=1:3
    cmap(:,ii) = pchip(x,cmapp(:,ii),xi);
end
cmap = flipud(cmap/255);



