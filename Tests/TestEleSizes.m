% Test size bounds with varying mesh gradation rates.

clearvars; clc; close all

RESO_TOL = 95; %percentage of resolution in bounds

bbox = [
    166 176;    % lon_min lon_max
    -48 -40     % lat_min lat_max
    ];
min_el = 1e3;   % minimum resolution in meters.
max_el = 100e3; % maximum resolution in meters.
max_el_ns = 5e3;    % maximum resolution nearshore in meters.
grade = [0.15; 0.25; 0.35]; % mesh grade in decimal percent.
R = 3;          % number of elements to resolve feature width.

coastline = 'GSHHS_f_L1';
gdat = geodata('shp',coastline,'bbox',bbox,'h0',min_el);

for i = 1:3 % for each grade
    fh = edgefx('geodata',gdat,...
        'fs',R,'max_el_ns',max_el_ns,...
        'max_el',max_el,'g',grade(i));
    mshopts = meshgen('ef',fh,'bou',gdat,'plot_on',0,'nscreen',5,'proj','trans');
    mshopts = mshopts.build;
    m1 = mshopts.grd;
    
    [bars,barlen] = GetBarLengths(m1,0);
    % sort bar lengths in descending order
    [barlen,IA] = sort(barlen,'descend');
    bars = bars(IA,:);
    % get the minimum bar length for each node
    [B1,IB] = unique(bars(:,1),'last');
    [B2,IC] = unique(bars(:,2),'last');
    d1 = NaN*m1.p(:,1); d2 = NaN*m1.p(:,1);
    d1(B1) = barlen(IB); d2(B2) = barlen(IC);
    reso = min(d1,d2);
    
    reso_in_bounds = 100*sum(reso > min_el & reso < max_el)/length(reso);
    if reso_in_bounds < RESO_TOL
        error(['Resolution bounds does not match for grade %4.2f. '...
            'Got %4.2f%% of vertices with resolution in bounds.'],...
            grade(i),reso_in_bounds);
        exit(1)
    end
    fprintf('Passed for %4.2f: %4.2f%% of vertices have resolution in bounds.\n',...
        grade(i),reso_in_bounds);
end
