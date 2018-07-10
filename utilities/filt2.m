function Zf = filt2(Z,res,lambda,filtertype) 
% filt2 performs a highpass, lowpass, bandpass, or bandstop 2D gaussian filter on gridded data such as 
% topographic, atmospheric, oceanographic, or any kind of geospatial data. This function is designed to
% make it easy to remove features longer or shorter than a given characteristic wavelength. The
% input grid can contain NaNs!  
% 
%% Syntax
% 
%  Zf = filt2(Z,res,lambda,filtertype)
% 
%% Description 
% 
% Zf = filt2(Z,res,lambda,filtertype) filters 2D dataset Z that has resolution res, 
% to an approximate wavelength lambda.  If the filtertype is 'lp' or 'hp' for lowpass
% or highpass, lambda must be a scalar value.  If the filtertype is 'bp' or 'bs' for 
% bandpass or bandstop, lambda must be a two-element array of the two cutoff wavelengths. 
% 
%% Explanation of this type of filter 
% For an explanation of this filter with pictures, type 
% 
%   showdemo filt2_documentation 
% 
% into your command window. 
% 
%% Example
% Consider a 100 km by 100 km elevation dataset that has a resolution of 200 m.  It has 
% some long 25 km wavelength features aligned with north/south direction, some short 
% ~5 km features oriented diagonally, and a considerable amount of random noise.  There's
% also a block of missing data. Here's what your datset looks like: 
% 
% res = 0.2; % 200 m resolution
% x = 0:res:100; % eastings from 0 to 100 km
% y = 0:res:100; % northings from 0 to 100 km
% [X,Y] = meshgrid(x,y);
% 
% % Z contains 25 km features, ~5 km diagonal features, and noise: 
% Z = cos(2*pi*X/25)+cos(2*pi*(X+Y)/7)+randn(size(X)); 
% 
% % Z also has some missing data: 
% Z(100:120,100:120) = nan; 
% 
% subplot(1,3,1) 
% imagesc(x,y,Z); 
% axis xy image 
% caxis([-1 1])
% title ' original Z matrix '
% xlabel ' eastings (km) '
% ylabel ' northings (km) '
% 
% % Get the lowpass-filtered version of Z: 
% Zlow = filt2(Z,res,15,'lp'); 
% 
% subplot(1,3,2) 
% imagesc(x,y,Zlow); 
% axis xy image 
% caxis([-1 1])
% title ' 15 km lowpass filtered Z matrix ' 
% xlabel ' eastings (km) '
% ylabel ' northings (km) '
% 
% % Get the highpass-filtered version of Z: 
% Zhi = filt2(Z,res,15,'hp'); 
% 
% subplot(1,3,3) 
% imagesc(x,y,Zhi); 
% axis xy image 
% caxis([-1 1])
% title ' 15 km highpass filtered Z matrix ' 
% xlabel ' eastings (km) '
% ylabel ' northings (km) '
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas Institute for 
% Geophysics in November 2016; however, all I did was repackage Carlos Adrian Vargas Aguilera's 
% superb ndnanfilter function, which can be found here: http://www.mathworks.com/matlabcentral/fileexchange/20417. 
% Many thanks to Carlos for his well-thought-out code and clear documentation.  
% 
% See also conv2, imgaussfilt, and imfilter.

%% Input checks

narginchk(4,4) 
% assert(license('test','image_toolbox')==1,'Error: I''m sorry, the filt2 function requires the Image Processing Toolbox.') 
assert(ismatrix(Z)==1,'Input error: Z must be a 2d matrix.')
assert(isscalar(res)==1,'Input error: res must be a scalar value.') 
assert(ismember(lower(filtertype),{'lp','hp','bp','bs'}),'Input error: filtertype must be ''hp'', ''lp'', ''bp'', or ''bs''.') 
if lambda<=(2*res) 
   warning('Nyquist says the wavelength should exceed two times the resolution of the dataset, which is an unmet condition based on these inputs. I''ll give you some numbers, but I would''t trust ''em if I were you.') 
end

if ismember(lower(filtertype),{'bp','bs'})
   assert(numel(lambda)==2,'Input error: Wavelength lambda must be a two-element array for a bandpass filter.') 
else
   assert(isscalar(lambda)==1,'Input error: Wavelength lambda must be a scalar for lowpass or bandpass filters.') 
end

%% Design filter: 

% 2*pi*sigma is the wavelength at which the amplitude is multiplied by a factor of about 0.6 (more exactly, exp(-0.5))
sigma = (lambda(1)/res) /(2*pi); 
%f = fspecial('gaussian',2*ceil(2.6*sigma)+1,sigma);
% WJP implementation without image processing toolbox, subfunc at bottom of
% this function
f = gaussian2D([2*ceil(2.6*sigma)+1 2*ceil(2.6*sigma)+1], sigma);

%% Now filter the data

switch lower(filtertype)
   %WJP incorrect usage originally but happened to work if imfilter is available
   case 'lp'
      Zf = ndnanfilter(Z,f,size(f),[],[],{'replicate'}); % ndnanfilter is Carlos Adrian Vargas Aguilera's excellent function, which is included as a subfunction below. 
      
   case 'hp'
      Zf = Z - ndnanfilter(Z,f,size(f),[],[],{'replicate'}); 
      
   case 'bp' 
      Zf =  filt2(filt2(Z,res,max(lambda),'hp'),res,min(lambda),'lp'); 
      
   case 'bs' 
      Zf = filt2(Z,res,max(lambda),'lp') - filt2(Z,res,min(lambda),'hp'); 
      
   otherwise 
      error('No such filter type.') 
end


end


function [Y,W] = ndnanfilter(X,HWIN,F,DIM,WINOPT,PADOPT,WNAN)
% NDNANFILTER   N-dimensional zero-phase digital filter, ignoring NaNs.
%
%   Syntax:
%         Y = ndnanfilter(X,HWIN,F);
%         Y = ndnanfilter(X,HWIN,F,DIM);
%         Y = ndnanfilter(X,HWIN,F,DIM,WINOPT);
%         Y = ndnanfilter(X,HWIN,F,DIM,WINOPT,PADOPT);
%         Y = ndnanfilter(X,HWIN,F,DIM,WINOPT,PADOPT,WNAN);
%     [Y,W] = ndnanfilter(...);
%
%   Input:
%     X      - Data to be filtered with/without NaNs.
%     HWIN   - Window function handle (or name) or numeric multidimensional
%              window to be used (without NaNs). See WINDOW for details.   
%              Default:   @rectwin  or 'rectwin' (moving average).
%     F      - A vector specifying the semi-width of the window for each
%              dimension. The final window's width will be 2*F+1. 
%              Default: 3 (i.e. a 1-dimensional window of width 6).
%     DIM    - If F is a single scalar, the window will be applied through
%              this dimension; otherwise, this will be ignored. 
%              Default: columns (or the first non-singleton dimension).
%     WINOPT - Cell array specifying optional arguments for the window
%              function HWIN (in addition to the width). 
%              Default: {} (window's defaults).
%     PADOPT - Cell array specifying the optional arguments for the
%              PADARRAY MATLAB's function (in addition to the array X and
%              the padsize: 2*F+1). If the function is not found, data is
%              padded with zeros or the specified value: try {mean(X(:))}
%              for example.
%              Default: {'replicate'} (repeats border elements of X).
%              Default: {0} (pads with zeros if PADARRAY not found).
%     WNAN   - Integer indicating NaNs treatment and program behaviour!:
%              0: Filters data and interpolates NaNs         (default). 
%              1: Filters data but do not interpolates NaNs 
%              2: "Do not filters data" but interpolates NaNs!
%              See the NOTEs below
%
%   Output:
%     Y      - Filtered X data (same size as X!).
%     W      - N-dimensional window with central symmetry generated by a
%              special subfunction called NDWIND. See the description below
%              for details.
%
%   Description:
%     This function applies a N-dimensional convolution of X with W, using
%     the MATLAB's IMFILTER or CONVN function. One important aspect of the
%     function is the generation of the N-dimensional window (W) from the
%     specified function and width, which cannot be done with MATLAB's
%     functions. Besides, unlike MATLAB's FILTER, FILTER2 and IMFILTER,
%     NaNs elements are taken into account (ignored).
%
%     The N-dimensional window is generated from rotating the 1-dimensional
%     output of the HWIN function, through each of the N-dimensions, and
%     then shrinking it through each of its axes in order to fit the
%     specified semi-widths (F). This is done in the included subfunction
%     named NDWIND. In this way, the window has central symmetry and do not
%     produce a phase shift on X data.
%
%     By default, the edges are padded with the values of X at the borders
%     with the PADARRAY MATLAB's function. In this way, the edges are
%     treated smoothly. When PADARRAY is not found, the program performs
%     zero-padding.
%
%   Notes: 
%     * The use of semi-widths F's is to force the generated window to be
%       even and, therefore, the change of phase is null.  
%     * The window function HWIN should output an even function, otherwise,
%       it won't generate an error but the user should be aware that this
%       program will consider only the last half of it.
%     * The function window should return a monotonically decreasing
%       result, this restriction is because I try to avoid the use of FZERO
%       function, for example, to find the expanding/shrinking factors.
%     * If the user has an already generated window, it can be used in HWIN
%       instead of a function handle or name. 
%     * Accepts empty value for any input. When X is empty, the program can
%       be used as a N-dimensional window generator.
%     * NaNs elements surrounded by no-NaNs elements (which will depend on
%       window width) are the ones that will be interpolated. The others
%       are leaved untouched.
%     * When WNAN=2, the programs acts like an NAN-interpolat/GAP-filling,
%       leaving untouched the no-NaNs elements but the filtering is
%       perfomed anyway. I recomend the default behaviour (WNAN=0) in order
%       to keep the filtered data in the workspace, and then use the code
%       at the end of this function to get/remove the interpolated NaNs
%       (see the example).
%     * The program looks for the IMFILTER and PADARRAY functions from the
%       Image Processing Toolbox. If not found, then CONVN is used instead
%       (slower) and pads with zeros or the given value. In this latter
%       case, if border elements are NaNs, the window won't work properly.
%
%   Example:
%     FWIN = 'hamming';
%     F = [13 8];
%     N = 100;
%     Pnoise = 0.30;
%     PNaNs  = 0.20;
%     X = peaks(N);                                     % original
%     Y = X + ((rand(size(X))-0.5)*2)*max(X(:))*Pnoise; % add noise
%     Y(round(1 + (N^2-1).*rand(N^2*PNaNs,1))) = NaN;   % add NaNs
%     [Z0,W] = ndnanfilter(Y,FWIN,F);                   % filters
%     Z1 = Z0; Z2 = Y; inan = isnan(Y);
%     Z1(inan) = NaN;
%     Z2(inan) = Z0(inan);  
%     subplot(231), imagesc(X), clim = caxis; axis equal tight
%                   title('Original data')
%     subplot(232), imagesc(Y),  caxis(clim), axis equal tight 
%                   title('Data + NOISE + NaNs')
%     subplot(234), imagesc(Z0), caxis(clim), axis equal tight 
%                   title('FILTERS + NaNs interpolation')
%     subplot(235), imagesc(Z1), caxis(clim), axis equal tight 
%                   title('FILTERS ignoring NaNs')
%     subplot(236), imagesc(Z2), caxis(clim), axis equal tight 
%                   title('GAP-filling with interpolated NaNs')
%     subplot(233), imagesc(-F(1):F(1),-F(2):F(2),W), axis equal tight, 
%                    title([upper(FWIN) ' 2D window']), view(2) 
%
%   See also: FILTER, FILTER2 and CONVN; WINDOW from the Signal Processing
%   Toolbox; and FWIND1, FWIND2, FSPECIAL, IMFILTER and PADARRAY from the
%   Image Processing Toolbox. 

%   Copyright 2008 Carlos Adrian Vargas Aguilera
%   $Revision: 1.2 $  $Date: 2008/06/30 18:00:00 $

%   Written by
%   M.S. Carlos Adrian Vargas Aguilera
%   Physical Oceanography PhD candidate
%   CICESE 
%   Mexico, 2008
%   nubeobscura@hotmail.com
%
%   Download from:
%   http://www.mathworks.com/matlabcentral/fileexchange/loadAuthor.do?objec
%   tType=author&objectId=1093874

%   1.0    Release (2008/06/23 10:30:00)
%   1.1    Fixed Bug adding an extra dimension of unitary width. 
%   1.2    Fixed Bug with ynan.

% Use the IMFILTER function? (faster than CONVN):
yimfilter = (exist('imfilter','file')==2);

% Use the PADARRAY function (or zero padding): 
ypadarray = (exist('padarray','file')==2);

% Check inputs and sets defaults of principal arguments:
if nargin<3 || nargin>7
 error('Filtern:IncorrectNumberOfInputs',...
  'At least three inputs are needed and less than 7.')
end
if isempty(HWIN)
 HWIN = 'rectwin';
end
if isempty(F)
 F = 3;
end
N = length(F);
S = size(X);
% Secondary arguments:
if N && (nargin<4 || isempty(DIM))
 DIM = find(S~=1,1);   %  DIM = min(find(S~=1));
 if isempty(DIM), DIM = 1; end
end
if nargin<5 || isempty(WINOPT)
 WINOPT = {};
end
if nargin<6 || isempty(PADOPT)
 if ypadarray
  PADOPT = {'replicate'};
 else
  PADOPT = {0};
 end
elseif ~ypadarray && ~isnumeric(PADOPT{1})
 PADOPT = {0};
end
if nargin<7 || isempty(WNAN)
 WNAN = 0;
end

% Selects the 1-dimensional filter or set a row vector: 
if N==1
 a = zeros(1,DIM);
 a(DIM) = F;
 F = a;
 clear a
end

% Checks if the window input is a function or an array:
if ~isa(HWIN,'function_handle') && ~ischar(HWIN)
 W = HWIN;
else
 W = [];
end

% If no input data but two outputs then generates the window only:
if isempty(X)
 Y = [];
 if nargout==2 && ~isempty(W)
  W = ndwind(HWIN,F,WINOPT{:});
 end
 return
end

% Generates the window:
if isempty(W)
 W = ndwind(HWIN,F,WINOPT{:});
end

% Check for NaN's:
inan = isnan(X);
ynan = any(inan(:));                       % Bug fixed 30/jun/2008
if ynan
 X(inan) = 0;
else
 factor = sum(W(:));
end

% Filtering:
if yimfilter                                % Use IMFILTER (faster)
 if ~isa(X,'double')
  X = double(X);
 end
 if ~isa(W,'double')
  W = double(W);
 end
 if ynan
  Y = imfilter(X,W       ,PADOPT{:},'conv');
 else
  Y = imfilter(X,W/double(factor),PADOPT{:},'conv');
 end
else                                        % Use CONVN
 % Sets F and S of equal sizes.
 F = reshape(F,1,N);
 Nx = numel(S);
 if N<Nx
  F(N+1:Nx) = 0;
 elseif N>Nx
  S(Nx+1:N) = 1;
 end
 F2 = 2*F;
 % Pads the borders:
 if ypadarray
  ind    = padarray(false(S),F2,true     );    % Index of the padding.
  Y      = padarray(X       ,F2,PADOPT{:},'both');
 elseif length(PADOPT{1})==1
  ind2 = cell(N,1);
  for n = 1:N
   ind2{n} = F2(n) + (1:S(n)).';
  end
  ind          = true(2*F2+S);
  Y            = repmat(PADOPT{1},2*F2+S);
  ind(ind2{:}) = false;
  Y(ind2{:})   = X;
 else % No padding at all
  Y    = X;
  ind  = false(S); 
  warning('Ndnanfilter:PaddingOption','Do not perfom any padding.')
 end
 % Convolutes both arrays:
 if ynan
  Y = convn(Y,W       ,'same');
 else
  Y = convn(Y,W/factor,'same');
 end
 %  Eliminates the padding:
 Y(ind) = [];
 Y      = reshape(Y,S);   
end

% Estimates the averages when NaNs are present:
if ynan
 if yimfilter
  factor       = imfilter(double(~inan),W,PADOPT{:},'conv');
 else
  if ypadarray
   factor      = padarray(~inan,F2,PADOPT{:});
  elseif length(PADOPT{1})==1 % (won't work properly with NaNs at borders)
   factor          = ind;
   factor(ind2{:}) = ~inan;
  else
   factor = ~inan;
  end
  factor      = convn(factor,W,'same');
  factor(ind) = [];
  factor      = reshape(factor,S);
 end
 Y = Y./factor;
end

% What about NaNs?:
if     WNAN == 1       % Leave NaNs elements untouched!
 Y(inan) = NaN;
elseif WNAN == 2       % Leave no-NaNs elements untouched!!!
 X(inan) = Y(inan);
 Y = X;
end  
end

function W = ndwind(HWIN,F,varargin)
% NDWIND Generate a N-Dimensional zero-phase window.
%
%   Syntax:
%     W = ndwind(HWIN,F);
%     W = ndwind(HWIN,F,OPT);
%
%   Input:
%     HWIN - Window function handle. See WINDOW for details. By default
%            uses: @rectwin (a rectangular window).
%     F    - A vector specifying the semiwidth of the window for each
%            dimension. The window's width will be 2*F+1. By default uses:
%            3 (i.e. a window of width 6). 
%     OPT  - Cell array specifying optional arguments for the window
%            function. By default uses: {[]} (window's defaults).
%
%   Output:
%     W    - N-Dimensional window with central symmetry.
%
%   Description:
%     In the axes of each dimension, W has a 1-D window defined as
%              feval(HWIN,2*F(n)+1), n = 1,...,N.
%     That is, they are defined by the same window function but have
%     different widths. So, this program creates another widther window (at
%     least 201 points), with the same definition, and finds how much the
%     former windows should be expanded in order to fit the latter one. 
%
%     Afterwards, the coordinates of every point are expanded accordingly
%     and the value of the window in those points are found by linear
%     interpolation with the bigger window. 
%
%     In resume, it is like rotating this big window through every
%     dimension and then shrinking it through each of its axes to fix the
%     specified widths.
%
%   Notes: 
%     * Because of the use of the semi-widths F's, all the generated
%       windows are even. Therefore the change of phase is null. 
%     * The window function HWIN should output an even function, otherwise,
%       it won't generate an error but this program will consider only the
%       last half of it.
%     * The window should be monotonically decreasing.
%     * Instead of the handle window, it can be given as a string:
%       'hamming' instead of @hamming, for example.
%     * Uses the MATLAB's function FUNC2STR.
%
%   Example:
%     W = ndwind(@hamming,[3 2])
%     % Results:
%     W =
%     
%              0         0    0.0800         0         0
%              0    0.1417    0.3100    0.1417         0
%              0    0.3966    0.7700    0.3966         0
%         0.0800    0.5400    1.0000    0.5400    0.0800
%              0    0.3966    0.7700    0.3966         0
%              0    0.1417    0.3100    0.1417         0
%              0         0    0.0800         0         0
%
%
%   See also: WINDOW from the Signal Processing Toolbox; and FWIND1,
%   FWIND2, and FSPECIAL from the Image Processing Toolbox.

%   Copyright 2008 Carlos Adrian Vargas Aguilera
%   $Revision: 1.1 $  $Date: 2008/06/26 19:30:00 $

%   Written by
%   M.S. Carlos Adrian Vargas Aguilera
%   Physical Oceanography PhD candidate
%   CICESE 
%   Mexico, 2008
%   nubeobscura@hotmail.com
%
%   Download from:
%   http://www.mathworks.com/matlabcentral/fileexchange/loadAuthor.do?objec
%   tType=author&objectId=1093874

%   1.0    Release (2008/06/23 10:30:00)
%   1.1    Fixed Bug adding an extra dimension of unitary width.  

% Check inputs:
if nargin<1 || isempty(HWIN)
 HWIN = 'rectwin';
end
if nargin<2 || isempty(F)
 F = 3;
end

% Rectangular wind?:
if isa(HWIN,'function_handle')
 HWIN = func2str(HWIN);
end
if strcmpi(HWIN,'rectwin')
 W = ones([2*F(:).'+1 1]);
 return
end

% Generate the BIG window (only the last half):
FBIG         = max([100; F(:)]);
BIGw         = feval(HWIN,2*FBIG+1,varargin{:});
BIGw(1:FBIG) = [];       % Deletes the first half.
rBIGw        = 0:FBIG;   % Window argument (distance).

% Axial windows widths:
N  = numel(F);
F  = reshape(F,1,N); 
F  = [F 0];             % BUG fixed by adding an extra dimension.
N  = N+1;
F2 = 2*F+1;


% Pre-allocates the final window and the expanded axis:
W  = zeros(F2);
An = cell(N,1);
Ae = An;

% Generates the index and expanded axes:
for n = 1:N
 
 % Generate temporally the window in the n-axis:
 wn = feval(HWIN,F2(n),varargin{:});
 
 % Finds the expansion factors (Note: the window should tends to zero):
 if F(n)
  piv = wn(end);
  ind = (BIGw == piv);
  if ~any(ind)
   ind1 = (BIGw >= piv); ind1 = length(ind1(ind1));
   ind2 = (BIGw <= piv); ind2 = length(ind2(~ind2))+1;
   if ind2>FBIG+1
    r = rBIGw(ind1);
   else
    r = interp1(BIGw([ind1 ind2]), rBIGw([ind1 ind2]),piv);
   end
  else
   r = rBIGw(ind);
  end
  Ef = r/F(n);
 else
  Ef = 1;
 end
 
 % Reversed index and expanded n-axis (for the following grid):
 An{n} = (F(n):-1:0);
 Ae{n} = An{n}*Ef;
 
end

% Estimates the expanded distances outside the axes (only at the 1st
% quarter):
% Note: In a 2-Dimensional matrix, by the 1st quarter of a matrix I mean
% the first 1/4 piece of the matrix after you divided it throuh the middle
% row and column. In N-dimensions it would be the 1st 1/2^N part.
gride4      = cell(N,1);
[gride4{:}] = ndgrid(Ae{:});
R4          = sqrt(sum(reshape([gride4{:}],prod(F+1),N).^2,2));

% Generates the window and linear index in the 1st quarter:
grid4     = cell(N,1);
[grid4{:}]= ndgrid(An{:});
in        = (R4<=rBIGw(end));           % Looks for elements inside window.
W4        = zeros(F+1);                 % 1st quarter of the window.
W4(in)    = interp1(rBIGw,BIGw,R4(in)); % Interpolates the window values.
for n=1:N                               % Linear index on the 1st quarter.
 grid4{n} = flip(grid4{n}+1,n);
end
ind4      = sub2ind(F2,grid4{:});

% Index of permutations to fill the N-D window:
np = 2^N-1;
ip = zeros(1,np);
for n = 1:N
 ini  = 2^(n-1);
 step = ini*2;
 ip(ini:step:np) = n;
end

% Fills the N-D window by flipping W4 and the index: 
ones4       = false(F2);    % Avoids using new FALSE function
ones4(ind4) = true;
W(ones4)    = W4;
for kp = ip
 W4         = flip(W4,kp);
 ones4      = flip(ones4,kp);
 W(ones4)   = W4;
end
end

function h = gaussian2D(siz, std)

    % create the grid of (x,y) values
    siz = (siz-1)./2;
    [x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));

    % analytic function
    h = exp(-(x.*x + y.*y)/(2*std*std));

    % truncate very small values to zero
    h(h<eps*max(h(:))) = 0;

    % normalize filter to unit L1 energy 
    sumh = sum(h(:));
    if sumh ~= 0
        h = h/sumh;
    end
end
