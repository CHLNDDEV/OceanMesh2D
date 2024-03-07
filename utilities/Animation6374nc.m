function Animation6374nc(varargin)
% Anim6374nc generate an animatiom of an ADCIRC 63 netCDF file pair
% Some idears come from https://github.com/BrianOBlanton/adcirc_util
%
% P/V pairs:
%     nc63           - Elevation Time Series at All Nodes in the Model Grid
%     nc74           - Wind Stress Time Series at All Nodes in the Model Grid
%     BestTrack      - best track dataset, JTWC or other format files
%     VideoOutput    - whether or not to output animation (def=True)
%     VideoWriteDir  - Path to save the Video
%     VideoProfile   - mp4 use 'MPEG-4', but only in Windows and macOS
%     FrameRate      - Rate of video playback in frames per second
%     ImageWriteDir  - Path to save each image/frame for every time step
%     ImageOutput    - whether or not to output images/frames (def=False)
%     ImageReso      - '-q95'(range 0-100,def=95), a larger val produces 
%                      higher quality and lower compression
%     AxisLims       - plot axis limits (def=mesh lims)
%     Title          - title string as a cell array (def={''})
%     IterStart      - starting iteration number (def=1)
%     IterStride     - iteration stride (skip) (def=1)
%     IterStop       - stopping iteration number (def=-1, represent end)
%     ColorMin       - min (def=min(zeta))
%     ColorMax       - max (def=max(zeta))
%     ColorMap       - colormap to use (def=jet(32))
%     Base_date      - datenum of starttime
%     FigurePosition - here use normalized position, [left bottom width height]
%     Projection     - see projections in m_map
%     FontSize       - default value is set as 16
%     N_Structural_G - dimension of structural grid to interpolate wind (def=100)
%     AddYourLogo    - add your logo to the plot at the bottom left corner 
%     
%     export_fig can be found here: https://github.com/altmany/export_fig
%     
%     Example:
%               Animation6374nc('nc63','./fort.63.nc', ...
%                               'nc74','./fort.74.nc', ...
%                               'BestTrack','./bio062007.txt', ...
%                               'Base_date','2007-11-10 06:00:00', ...
%                               'IterStart',1, ...
%                               'IterStop',-1, ...
%                               'VideoOutput',1, ...
%                               'VideoWriteDir','./', ...
%                               'ImageOutput',1, ...
%                               'ImageWriteDir','./')
%
%  Author:      Jiangchao Qiu, (MIT/ESSG; email:qiujch24@mit.edu)
%  Created:     May 14 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.nc63           = '';
p.nc74           = '';
p.BestTrack      = '';
p.VideoOutput    = 1;
p.VideoWriteDir  = './';
p.VideoProfile   = 'MPEG-4';
p.FrameRate      = 5;
p.ImageWriteDir  = './';
p.ImageOutput    = 0;
p.ImageReso      = '-q95';
p.AxisLims       = [];  
p.Title          = {};     
p.IterStart      = 1; 
p.IterStride     = 1;
p.IterStop       = -1;  
p.ColorMin       = [];  
p.ColorMax       = [];  
p.ColorMap       = jet(32);  
p.Base_date      = [];
p.FigurePosition = [0.1,0.25,0.7,0.75];
p.Projection     = 'Miller';
p.FontSize       = 16;
p.N_Structural_G = 90;
p.AddYourLogo    = 'Plot by MIT/ESSG';

p = parse_pv_pairs(p,varargin);

%% read into data from nc format file
x            =  ncread(p.nc63,'x');
y            =  ncread(p.nc63,'y');
element      =  ncread(p.nc63,'element')';
zeta         =  ncread(p.nc63,'zeta');
time         =  ncread(p.nc63,'time');
windx        =  ncread(p.nc74,'windx');
windy        =  ncread(p.nc74,'windy');
% find the max and min zeta value as the up and down limit of colormap
min_zeta     =  min(min(zeta)); p.ColorMin = min_zeta;
max_zeta     =  max(max(zeta)); p.ColorMax = max_zeta;
% convert seconds since base date into hours since base date 
% base_date shoule be set as the cold start date in fort.15 file
time         =  time/3600; % hour time series since the start time
Base_date    =  datetime(p.Base_date,'InputFormat','yyyy-MM-dd HH:mm:ss');

nTimes = length(time);

if p.IterStop==-1
    p.IterStop=nTimes;
end

%% read into JTWC best track data
tc = readtable(p.BestTrack);
tc_lon = table2cell(tc(:,8));
tc_lat = table2cell(tc(:,7));
tc_lon_lat = zeros(56,2);
for i = 1:length(tc_lon)
    tc_lon_t = str2double(cell2mat(extractBefore(tc_lon(i,1),"E")))/10;
    tc_lat_t = str2double(cell2mat(extractBefore(tc_lat(i,1),"N")))/10;
    tc_lon_lat(i,1) = tc_lon_t;
    tc_lon_lat(i,2) = tc_lat_t;
end

%% we use VideoWriter to make the Animation
Animation = VideoWriter([p.VideoWriteDir,'video']);
Animation.FrameRate = p.FrameRate;
open(Animation)

disp("############")
disp("Start to Generate Image for Each Time Step:")
disp("############")

%% plot image for each timestep and read as one frame
% loop for each time step
fig = figure;
for idx = p.IterStart:p.IterStride:p.IterStop
    % read into data for current date
    current_date  = Base_date + time(idx,:)/24;
    disp(['Image is Ploting for Date:',datestr(current_date)])
    current_ele = zeta(:,idx); 
    current_windx = windx(:,idx);
    current_windy = windy(:,idx);
    % set figure properties
    set(gcf,'unit','normalized','position',p.FigurePosition)
    set(gca,'FontSize',p.FontSize);
    m_proj(p.Projection,'lon',[min(x)-0.25 max(x)+0.25],'lat',[min(y)-0.25 max(y)+0.25]);
    m_grid('linestyle','none','box','fancy','tickdir','in'); % 
    %set(gcf, 'color', 'none');    
    %set(gca, 'color', 'none');
    %m_grid('linestyle','none','box','off','backgroundcolor','none','ticklen',0,'xticklabels',[],'yticklabels',[]);
    title(p.Title,FontSize=p.FontSize+2);
    m_northarrow(82,21,1,'type',2);
    %vecscl=0.1;
    %m_vec(vecscl,81.3,19.8,-0.08,0,'k','shaftwidth',1,'headlength',5,'key',{'40 m/s','wind stress'});
    %add your logo
    m_text(81,10,p.AddYourLogo,'fontsize',p.FontSize,'color','k');
    caxis([p.ColorMin p.ColorMax]);
    colormap(p.ColorMap);
    cb = colorbar;
    set(cb.XLabel,{'String','Rotation'},{'m',0}); 
    %%
    m_trisurf(element,x,y,current_ele);
    m_line(tc_lon_lat(:,1),tc_lon_lat(:,2),'linestyle','--','linewi',3,'color',"#A2142F");
    %% interpolate unstructural wind values into structural/uniform grid and plot as quiver
    [x_grid_indom,y_grid_indom,windx_grid_indom] = interpolant_nc74(x,y,current_windx,p.N_Structural_G);
    [~           ,~           ,windy_grid_indom] = interpolant_nc74(x,y,current_windy,p.N_Structural_G);
    m_quiver(x_grid_indom,y_grid_indom,windx_grid_indom,windy_grid_indom,1,'color','k','linewidth',0.08);
    % whether or not make animation/video (default is true)
    if p.VideoOutput
        frame = getframe(gcf);
        writeVideo(Animation,frame);
    end
    % whether or not write out each image/frame 
    if p.ImageOutput
        export_fig([p.ImageWriteDir,'time_',num2str(idx)],'-png','-transparent',p.ImageReso);
    end
    clf;
end

if p.VideoOutput
    disp("Animation done!")
end

if p.ImageOutput
    disp("Image saved!")
end

close(fig)
close(Animation);

end 
