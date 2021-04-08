function N=m_windrose(long,lat,thet,spd,varargin)
% M_WINDROSE Rose Histogram plot for wind statistics
%      N=M_WINDROSE(LONG,LAT,ANGLE,SPEED) plots a  rose histogram plot for 
%      wind speed and direction vectors ANGLE (degrees) and SPEED on a 
%      map at location (LONG,LAT). 
%
%      M_WINDROSE(...,PARAM,VALUE) lets you specific additional properties
%      of the windrose. These include  
%
%     'nspeeds': The number of speed bins. Scalar or vector (default 
%                [0:4:24]). If a scalar then this sets the number of 
%                bins between 0 and the maximum speed. If the lower
%                limit is >0, then a 'calm' circle is placed in the
%                center. WRMING: Speeds greater than MAX(NSPEEDS) are not 
%                included in the binning.
%     'ndirs'  : The number of direction bins. Scalar (default 36)
%     'size'  : The fraction of the map width or height covered by the
%               background rings. Scalar (fraction, default 0.075)
%     'nrings' : The number of background rings. Scalar or vector 
%                (default [2 4 6]). If scalar, set to the number of rings
%                spaced 2% apart.
%     'labelrings':  logical or angle
%     'barstyle'  : 'colour' (default),'checkered'.
%     'alpha'  : transparency of ring (0 to 1, default 0.4). Useful to see 
%                coastlines behind, or when stations are crowded together.
%     'parts'  : 'bars', 'background', 'all' (default). Not usually needed
%                for a user call, but it may be useful to draw all
%                'background' rings first, then all 'bars' on top, to improve
%                visibility when roses overlap.
%
%    Although in the simplest call, LONG,LAT are scalars and ANGLE, SPEED
%    are vectors, it is also possible to vectorize calls with:
%        - LON,LAT - vectors of length N. ANGLE,SPEED as MxN matrices (this
%          can be used if all wind time series have the same length)    
%        - LON,LAT,ANGLE,SPEED - cell arrays of length N (this can be
%          useful if wind time series have different lengths)
%
 

% Rich Pawlowicz Feb/2020


% Deform the circle to fit the map projection
global MAP_PROJECTION MAP_VAR_LIST

if isempty(MAP_PROJECTION)
   disp('No Map Projection initialized - call M_PROJ first!');
   return;
end


% Handle struct array inputs
N=length(long);
if iscell(long)
    if N==length(lat) && N==length(thet) && N==length(spd)  
       for k=1:N
           m_windrose(long{k},lat{k},thet{k},spd{k},varargin{:},'part','background');
       end
       for k=1:N
           m_windrose(long{k},lat{k},thet{k},spd{k},varargin{:},'part','bar');
       end
       return;
    else
        error('Not all input parameters are cell arrays of the same size');
    end
elseif isvector(long) & N>1
    % Handle matrix inputs by calling (scalar) handling function in a loop
    if  N==length(lat) && N==size(thet,2) && N==size(spd,2)  
       for k=1:N
           m_windrose(long(k),lat(k),thet(:,k),spd(:,k),varargin{:},'part','background');
       end
       for k=1:N
           m_windrose(long(k),lat(k),thet(:,k),spd(:,k),varargin{:},'part','bar');
       end
       return;
    else
        error('Not all input parameters are arrays with the same number of columns');
    end
end

if ~( length(long)==1 && length(lat)==1 && length(thet)==length(spd) )
    error('Input parameters LONG/LAT/ANGLE/SPEED have non-conforming sizes');
end
    
% Defaults

ndirs=36;
nspds=[0:4:20];
rsize=.075;
nrings=[2 4 6];
labelrings=false;
barstyle='col';
alpha=.4;
parts='all';

k=1;
while ~isempty(varargin)
    switch lower(varargin{k}(1:3))
        case 'nsp'
            nspds=varargin{k+1};
            if length(nspds)==1
              nspds=[0:nspds-1]'/(nspds-1)*ceil(max(spd(:)));
            end
        case 'ndi'
            ndirs=varargin{k+1};
        case 'siz'
            rsize=varargin{k+1};
        case 'nri'
            nrings=varargin{k+1};
        case 'lab'
            labelrings=varargin{k+1};
        case 'bar'
            barstyle=varargin{k+1};
        case 'alp'
            alpha=varargin{k+1};
        case 'par'
            parts=varargin{k+1};
        otherwise
            error(['Unknown input option ''' varargin{k} '''']);
    end
    varargin([k k+1])=[];
end
        
if length(nrings)==1
    nrings=[0:nrings]*2;
end

if length(nspds)==1
    nspds=[0:nspds]*nanmax(spd);
else
    if max(spd)>max(nspds)
        percnt=nansum(spd>max(nspds))/length(spd)*100;
        disp(sprintf('WARNING: %.1f%% of speeds in time series exceed the binned range',percnt));
    end
end

% If we have data with a speed of zero, it often also
% has a direction of zero which can get weird as it *may*
% leave a large spike in the zero degree direction. So,
% spread out the data in direction space for speeds of zero.
if nspds(1)==0
    ii=spd==0;
    if any(ii)
        thet(ii)=[0:sum(ii)-1]/sum(ii)*360;
    end
end


ang=360/ndirs;
ndirs=[-ang/2:ang:360]; % Make the intervals bracket 0 degrees (north)

% Manipulate the angles so they fall in the bins we have set.
ii=thet<ndirs(1);
if any(ii)
    thet(ii)=thet(ii)+360;
end
ii=thet>ndirs(end);
if any(ii)
    thet(ii)=thet(ii)-360;
end;


%ndirs=[0:ang:360];
cdirs=ndirs(1:end-1)+diff(ndirs)/2;
cspds=nspds(1:end-1)+diff(nspds)/2;

% Get center and local x/y axis vectors
[X,Y]=m_ll2xy(long,lat,'clip','off');
[XN ,YN ]=m_ll2xy([long long]',[lat lat+.001]','clip','off');
[XE ,YE ]=m_ll2xy([long long+(.001)./cos(lat*pi/180)]',[lat lat]','clip','off');

% Sets the size of the ring
maxrad=min(diff(MAP_VAR_LIST.xlims),diff(MAP_VAR_LIST.ylims))*rsize/2;

%%%%%%%%%%%%%

if strcmp(parts(1:3),'all') | strcmp(parts(1:3),'bac')
    % This is for the circular rings
    u=maxrad*sind([0:5:360]);
    v=maxrad*cosd([0:5:360]);
    mU=u*diff(XE)*1000 + v*diff(XN)*1000;
    mV=u*diff(YE)*1000 + v*diff(YN)*1000;
    scfacB=(mU+i*mV);
    scfacB=scfacB/max(abs(scfacB))*maxrad;

    % Draw the background stuff first

    % The white circle
    patch(X+real(scfacB),Y+imag(scfacB),ones(size(scfacB)),'w','edgecolor',[1 1 1]*(1-alpha),...
           'facealpha',alpha,'edgealpha',alpha,'clip','off');
    % The specified rings
    for k=1:length(nrings)-1
        line(X+real(scfacB)*nrings(k)/nrings(end),Y+imag(scfacB)*nrings(k)/nrings(end),'color',[1 1 1]*(1-alpha) ,'clip','off');
    end
    % Lines through the center of the rings
    for k=0:5
       line(X+real(scfacB([1 37]+6*k)),Y+imag(scfacB([1 37]+6*k)),'color',[1 1 1]*(1-alpha),'clip','off' );
    end

end

%%%%%%%%%%

if strcmp(parts(1:3),'all') | strcmp(parts(1:3),'bar')

    % This is for the data bins
    u=maxrad*sind(cdirs);
    v=maxrad*cosd(cdirs);

    mU=u*diff(XE)*1000 + v*diff(XN)*1000;
    mV=u*diff(YE)*1000 + v*diff(YN)*1000;
    scfacS=(mU+i*mV);
    scfacS=scfacS/max(abs(scfacS))*maxrad;

    % Pick out the "calm" winds
    if nspds(1)>0
        kk=spd<nspds(1);
        Ncalm=sum(kk);
    else
        Ncalm=0;
    end

    % Bin all the data
    [N,Xedges,Yedges]=histcounts2(spd,thet,nspds,ndirs);

    [Ns,Na]=size(N);

    % Now, we juggle the lowermost bin to spread it's contents out equally
    % over all angles.
    %N(1,:)=sum(N(1,:))/size(N,2)*ones(1,size(N,2));

    % Convert to a cumulative frequency in %
    N=cumsum([repmat(Ncalm/Na,1,Na);N])/length(spd)*100;

    % Now, N is a matrix in which speeds are binned into the various row, and angles
    % are binned into the various columns. What we have to do now is to construct
    % a large matrix - each column of which represent the patch data associated with
    % a given element of N.

    % First the angles, which are the same for all rows.
    scfac=repmat(scfacS,Ns,1);

    % Now the widening of bars, which is the same for all columns, and scaled 
    % down slightly so that things look "nice"

    yfac=  tand(ang/2)*.9 ;

    % Each bar goes from the first element to the next.
    N1=N(1:end-1,:);
    N2=N(2:end,:);

    % Now strip N matrix-size into a long row.
    N1=N1(:)';
    N2=N2(:)';
    scfac=scfac(:).';
    yfac=yfac(:)';

    % For each element in row, draw a rectangle (complex math with scfac gives rotations)

    zz=([N1;N2;N2;N1;N1]+i*[-N1;-N2;N2;N1;-N1].*yfac  ).*scfac(ones(5,1),:)/nrings(end);

    % Now draw all the bars

    if strcmp(barstyle(1:3),'col')
       cfac=repmat(nspds(1:end-1),1,Na);
       cfac=cfac(:)';
       patch(X+real(zz),Y+imag(zz),ones(size(zz))+1, repmat(cfac,5,1),'edgecolor','none','clip','off' );
    else
       % alternate the colours from white to black with the variour rows
       cfac=repmat(rem([0:Ns-1]',2),1,Na);
       cfac=repmat(cfac(:)',1,1,3);
       patch(X+real(zz),Y+imag(zz),ones(size(zz)),cfac,'edgecolor','k','clip','off');
    end


    % The central circle

    if Ncalm>0 
       patch(X+real(N(1,:).*scfacS/nrings(end)),Y+imag(N(1,:).*scfacS/nrings(end)),'w','edgecolor','k');
       text(X,Y,3,sprintf('%.0f%%',Ncalm/length(spd)*100),'horizontal','center','vertical','middle',...
         'fontsize',10);

    end

end;


