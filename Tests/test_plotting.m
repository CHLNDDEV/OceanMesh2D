% Test plotting 
clearvars; clc;

addpath('..')
addpath(genpath('../utilities/'))
addpath(genpath('../datasets/'))
addpath(genpath('../m_map/'))

bbox = [166 176;		% lon_min lon_max
        -48 -40]; 		% lat_min lat_max
min_el    = 1e3;  		% minimum resolution in meters.
max_el    = 100e3; 		% maximum resolution in meters. 
max_el_ns = 5e3;        % maximum resolution nearshore in meters.
grade     = 0.35; 		% mesh grade in decimal percent.
R         = 3;    		% number of elements to resolve feature width.
coastline = 'GSHHS_f_L1';
gdat = geodata('shp',coastline,'bbox',bbox,'h0',min_el);
fh = edgefx('geodata',gdat,...
            'fs',R,'max_el_ns',max_el_ns,...
            'max_el',max_el,'g',grade);
mshopts = meshgen('ef',fh,'bou',gdat,'plot_on',1,'nscreen',5,'proj','trans');
mshopts = mshopts.build; 
m = mshopts.grd;
m = make_bc(m,'auto',gdat,'distance'); % make the boundary conditions


bbox_s =  [172   176;
           -42   -39];
       
close all; 
%%
m.b = rand(length(m.p(:,1)),1); 
m.bx = rand(length(m.p(:,1)),1); 
m.by = rand(length(m.p(:,1)),1); 

% start testing all types
types = {'tri','bd','ob','b','slp','reso','resodx','qual'}; 
add   = {'notri','earth','log','mesh'}; 
projs = {}; 
for i = 1 : length(types)
    for j = 1 : length(add)
        test = strcat(types{i},add{j}); 
        disp(['Testing plotting option: ',test]); 
        plot(m,'type',types{i});
        pause; disp('Look right?'); close all; 
    end
end
