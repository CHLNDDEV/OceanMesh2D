function Processed=ProcessSignal(x,y,DerivativeMode,w,type,ends,Sharpen,factor1,factor2,SlewRate,MedianWidth)
% Command line function that performs smoothing and differentiation on the
% time-series data set x,y. Returns the processed signal as a vector that
% has the same shape as x, regardless of the shape of y.
% Version 1.0, April, 2013.  Tom O'Haver (toh@umd.edu)
%
% 'DerivativeMode' determines the derivative order (O, 1, 2, 3, 4, 5); 
%    See http://terpconnect.umd.edu/~toh/spectrum/Differentiation.html
%
% 'w' is the smooth width; 
% 'type' determines the smooth mode:
%        If type=0, the signal is not smoothed.
%        If type=1, rectangular (sliding-average or boxcar) 
%        If type=2, triangular (2 passes of sliding-average)
%        If type=3, pseudo-Gaussian (3 passes of sliding-average)
%        If type=4, Savitzky-Golay smooth 
% 'ends' controls how the "ends" of the signal (the first w/2 points and the last w/2 points) are handled.
%        If ends=0, the ends are zeroed
%        If ends=1, the ends are smoothed with progressively 
%        smaller smooths the closer to the end.
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html
%
% 'Sharpen' controls peak sharpening (resolution enhancement). 0=off; 1=on.
%   The sharpening str5ength is controled by the factor1 and factor2 The
%   optimum values depend on the peak shape and width. For details, see
%   http://terpconnect.umd.edu/~toh/spectrum/InteractiveResEnhance.htm). 
%  
% 'SlewRate' sets the maximum slew rate (maximum rate of change of the
%    signal when it is non-zero (useful for reducing the effect of large 
%    steps in the background).
% 
% 'MedianWidth' implements a median filter when it is non-zero, which
%    replaces each point in the signal with the median (rather than the
%    average) of MedianWidth adjacent points, cuseful for dealing with narrow
%    spike artifacts (set MedianWidth to the number of points in the spikes)
%
% EXAMPLE 1:
% >> ProcessSignal([1 3 5 7 9],[0 0 1 0 0],0,3,1,0,0,0,0,0,0)
% ans =
%             0      0.33333      0.33333      0.33333            0
% EXAMPLE 2:
% Use of smoothed second derivative to reduce influence of background.
% on the weak peak at x=5.
% >> x=[0:.01:10]';
% >> y=10./(1+(x/4).^2)+exp(-(x-5).^2)+.1.*randn(size(x));
% >> subplot(2,1,1);plot(x,y);
% >> subplot(2,1,2);plot(x,ProcessSignal(x,y,2,60,3,0,0,0,0,0,0))
%
% EXAMPLE 3: Single peak with random spikes. Use of median filter
% (MedianWidth=1) to remove single-point spikes.
% >> x=-5:.01:5;
% >> y=exp(-(x).^2);for n=1:1000,if randn()>2,y(n)=rand()+y(n);,end,end;
% >> plot(x,ProcessSignal(x,y,0,0,0,0,0,0,0,0,1))
%
if SlewRate,
    for n=1:length(y)-1,
        if y(n+1)-y(n)>SlewRate,y(n+1)=y(n)+SlewRate;end
        if y(n+1)-y(n)<-SlewRate,y(n+1)=y(n)-SlewRate;end
    end
end % SlewRate
if MedianWidth,
    mY=y;
    for n=1:length(x)-(1+MedianWidth*2),
        mY(n+MedianWidth)=median(y((n):(n+1+MedianWidth*2)));
        y=mY;
    end
end  % MedianWidth
if type==0,w=1;end
if type==4,
    if w<2,w=3;end
    if DerivativeMode>4,
        if w<5,w=5;end
    end 
    % The polynomial order, 2+DerivativeMode, must be less than the
    % frame size, 2*w+1, and 2*w+1 must be odd.  
        Processed=savitzkyGolayFilt(y,2+DerivativeMode,DerivativeMode,2*w+1);
        if DerivativeMode==1,Processed=-Processed;end
        if DerivativeMode==3,Processed=-Processed;end;
else
    switch DerivativeMode
        case 0
            Processed=fastsmooth(y,w,type,ends);
        case 1
            Processed=fastsmooth(deriv(x,y),w,type,ends);
        case 2
            D2=fastsmooth(secderiv(x,y),w,1,ends);
            Processed=fastsmooth(D2,w,type,ends);
        case 3
            D2=fastsmooth(secderiv(x,y),w,1,ends);
            D3=fastsmooth(deriv(x,D2),w,1,ends);
            Processed=fastsmooth(D3,w,type,ends);
        case 4
            D2=fastsmooth(secderiv(x,y),w,1,ends);
            D4=fastsmooth(secderiv(x,D2),w,1,ends);
            D4=fastsmooth(D4,w,1,ends);
            
            Processed=fastsmooth(D4,w,type,ends);
        case 5
            D2=fastsmooth(secderiv(x,y),w,1,ends);
            D4=fastsmooth(secderiv(x,D2),w,1,ends);
            D4=fastsmooth(D4,w,1,ends);
            Processed=fastsmooth(deriv(x,D4),w,type,ends);
    end
end
if Sharpen,
    type=4; 
    if w<3;w=3;end
    Processed=enhance(x,Processed,factor1,factor2,w,type);
end
Processed=reshape(Processed,size(x));
% ----------------------------------------------------------------------
function Enhancedsignal=enhance(x,signal,factor1,factor2,SmoothWidth,type)
% Resolution enhancement function by even derivative method. The
% arguments factor1 and factor 2 are 2nd and 4th derivative weighting
% factors. Larger values of factor1 and factor2 will reduce the 
% peak width but will cause artifacts in the baseline near 
% the peak.  Adjust the factors for the the best compromise. 
% Use minimum smooth width needed to reduce excess noise. 
datasize=size(signal);
if datasize(1)>datasize(2),signal=signal';end
if type==4,
    Enhancedsignal = signal-factor1.*savitzkyGolayFilt(signal,4,2,2*SmoothWidth+1)+...
        factor2.*savitzkyGolayFilt(signal,6,4,2*SmoothWidth+1);
else
    d2=secderiv(x,signal);  % Computes second derivative
    d4=secderiv(x,d2);   % Computes fourth derivative
    Enhancedsignal = signal-factor1.*fastsmooth(d2,SmoothWidth,type)+...
        factor2.*fastsmooth(fastsmooth(fastsmooth(d4,SmoothWidth,2),SmoothWidth,2),SmoothWidth,2);
end
Enhancedsignal=Enhancedsignal';
% ----------------------------------------------------------------------
function d=secderiv(x,a)
% Second derivative of y with respect to x using 3-point central difference.
%  T. C. O'Haver, 2011.
n=length(a);
d=zeros(size(a));
for j = 2:n-2;
  x1=x(j-1);x2=x(j);x3=x(j+1);
  d(j)=((a(j+1)-a(j))./(x3-x2) - (a(j)-a(j-1))./(x2-x1))./((x3-x1)/2);
end
d(1)=d(2);
d(n)=d(n-1);
% ----------------------------------------------------------------------
function d=deriv(x,y)
% First derivative of y with respect to x using 2-point central difference.
%  T. C. O'Haver, 2011.
n=length(y);
d=zeros(size(y));
d(1)=(y(2)-y(1))./(x(2)-x(1));
d(n)=(y(n)-y(n-1))./(x(n)-x(n-1));
for j = 2:n-1;
  d(j)=(y(j+1)-y(j-1)) ./ (1.*(x(j+1)-x(j-1)));
end
% ----------------------------------------------------------------------
function SmoothY=fastsmooth(Y,w,type,ends)
% fastbsmooth(Y,w,type,ends) smooths vector Y with smooth 
%  of width w. Version 2.0, May 2008.
% The argument "type" determines the smooth type:
%   If type=1, rectangular (sliding-average or boxcar) 
%   If type=2, triangular (2 passes of sliding-average)
%   If type=3, pseudo-Gaussian (3 passes of sliding-average)
% The argument "ends" controls how the "ends" of the signal 
% (the first w/2 points and the last w/2 points) are handled.
%   If ends=0, the ends are zero.  (In this mode the elapsed 
%     time is independent of the smooth width). The fastest.
%   If ends=1, the ends are smoothed with progressively 
%     smaller smooths the closer to the end. (In this mode the  
%     elapsed time increases with increasing smooth widths).
% fastsmooth(Y,w,type) smooths with ends=0.
% fastsmooth(Y,w) smooths with type=1 and ends=0.
% Example:
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3)= [0 1 4 7 10 7 4 1 1 0]
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3,1,1)= [1 1 4 7 10 7 4 1 1 1]
%  T. C. O'Haver, May, 2008.
if nargin==2, ends=0; type=1; end
if nargin==3, ends=0; end
  switch type
    case 0
       SmoothY=sa(Y,w,ends);  
    case 1
       SmoothY=sa(Y,w,ends);
    case 2   
       SmoothY=sa(sa(Y,w,ends),w,ends);
    case 3
       SmoothY=sa(sa(sa(Y,w,ends),w,ends),w,ends);
  end
function SmoothY=sa(Y,smoothwidth,ends)
w=round(smoothwidth);
SumPoints=sum(Y(1:w));
s=zeros(size(Y));
halfw=round(w/2);
L=length(Y);
for k=1:L-w,
   s(k+halfw-1)=SumPoints;
   SumPoints=SumPoints-Y(k);
   SumPoints=SumPoints+Y(k+w);
end
s(k+halfw)=sum(Y(L-w+1:L));
SmoothY=s./w;
% Taper the ends of the signal if ends=1.
  if ends==1,
    startpoint=(smoothwidth + 1)/2;
    SmoothY(1)=(Y(1)+Y(2))./2;
    for k=2:startpoint,
       SmoothY(k)=mean(Y(1:(2*k-1)));
       SmoothY(L-k+1)=mean(Y(L-2*k+2:L));
    end
    SmoothY(L)=(Y(L)+Y(L-1))./2;
  end% ----------------------------------------------------------------------
function y=savitzkyGolayFilt(x,N,DN,F,W,DIM)
%savitzkyGolayFilt Savitzky-Golay Filtering.
%   savitzkyGolayFilt(X,N,DN,F) filters the signal X using a Savitzky-Golay 
%   (polynomial) filter.  The polynomial order, N, must be less than the
%   frame size, F, and F must be odd.  DN specifies the differentiation
%   order (DN=0 is smoothing). For a DN higher than zero, you'll have to
%   scale the output by 1/T^DN to acquire the DNth smoothed derivative of
%   input X, where T is the sampling interval. The length of the input X
%   must be >= F.  If X is a matrix, the filtering is done on the columns
%   of X.
%
%   Note that if the polynomial order N equals F-1, no smoothing
%   will occur.
%
%   savitzkyGolayFilt(X,N,DN,F,W) specifies a weighting vector W with
%   length F containing real, positive valued weights employed during the
%   least-squares minimization. If not specified, or if specified as
%   empty, W defaults to an identity matrix.
%
%   savitzkyGolayFilt(X,N,DN,F,[],DIM) or savitzkyGolayFilt(X,N,DN,F,W,DIM)
%   operates along the dimension DIM.
%   Copyright (c) 2011, Diederick
%   See also savitzkyGolay, FILTER, sgolayfilt

%   References:
%     [1] Sophocles J. Orfanidis, INTRODUCTION TO SIGNAL PROCESSING,
%              Prentice-Hall, 1995, Chapter 8.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.11.4.4 $  $Date: 2009/08/11 15:47:54 $

error(nargchk(4,6,nargin,'struct'));

% Check if the input arguments are valid
if round(F) ~= F, error(generatemsgid('MustBeInteger'),'Frame length must be an integer.'), end
if rem(F,2) ~= 1, error(generatemsgid('SignalErr'),'Frame length must be odd.'), end
if round(N) ~= N, error(generatemsgid('MustBeInteger'),'Polynomial order must be an integer.'), end
if N > F-1, error(generatemsgid('InvalidRange'),'The Polynomial order must be less than the frame length.'), end
if DN > N, error(generatemsgid('InvalidRange'),'The Differentiation order must be less than or equal to the Polynomial order.'), end

if nargin < 5 || isempty(W)
   % No weighting matrix, make W an identity
   W = ones(F,1);
else
   % Check for right length of W
   if length(W) ~= F, error(generatemsgid('InvalidDimensions'),'The weight vector must be of the same length as the frame length.'),end
   % Check to see if all elements are positive
   if min(W) <= 0, error(generatemsgid('InvalidRange'),'All the elements of the weight vector must be greater than zero.'), end
end

if nargin < 6, DIM = []; end

% Compute the projection matrix B
pp = fix(-F./2):fix(F./2);
B = savitzkyGolay(pp,N,DN,pp,W);

if ~isempty(DIM) && DIM > ndims(x)
	error(generatemsgid('InvalidDimensions'),'Dimension specified exceeds the dimensions of X.')
end

% Reshape X into the right dimension.
if isempty(DIM)
	% Work along the first non-singleton dimension
	[x, nshifts] = shiftdim(x);
else
	% Put DIM in the first dimension (this matches the order 
	% that the built-in filter function uses)
	perm = [DIM,1:DIM-1,DIM+1:ndims(x)];
	x = permute(x,perm);
end

if size(x,1) < F, error(generatemsgid('InvalidDimensions'),'The length of the input must be >= frame length.'), end

% Preallocate output
y = zeros(size(x));

% Compute the transient on (note, this is different than in sgolayfilt,
% they had an optimization leaving out some transposes that is only valid
% for DN==0)
y(1:(F+1)/2-1,:) = fliplr(B(:,(F-1)/2+2:end)).'*flipud(x(1:F,:));

% Compute the steady state output
ytemp = filter(B(:,(F-1)./2+1),1,x);
y((F+1)/2:end-(F+1)/2+1,:) = ytemp(F:end,:);

% Compute the transient off
y(end-(F+1)/2+2:end,:) = fliplr(B(:,1:(F-1)/2)).'*flipud(x(end-(F-1):end,:));

% Convert Y to the original shape of X
if isempty(DIM)
	y = shiftdim(y, -nshifts);
else
	y = ipermute(y,perm);
end
% ----------------------------------------------------------------------
function [fc, df] = savitzkyGolay(x,n,dn,x0,W,flag)
% Function:
%       Savitzky-Golay Smoothing and Differentiation Filter
%       Copyright (c) 2011, Diederick
%       The Savitzky-Golay smoothing/differentiation filter (i.e., the
%       polynomial smoothing/differentiation filter, or  the least-squares
%       smoothing/differentiation filters) optimally fit a set of data
%       points to polynomials of different degrees. 
%       See for details in Matlab Documents (help sgolay). The sgolay
%       function in Matlab can deal with only symmetrical and uniformly
%       spaced data of even number.
%       This function presented here is a general implement of the sgolay
%       function in Matlab. The Savitzky-Golay filter coefficients for even
%       number, nonsymmetrical and nonuniformly spaced data can be
%       obtained. And the filter coefficients for the initial point or the
%       end point can be obtained too. In addition, either numerical
%       results or symbolical results can be obtained. Lastly, this
%       function is faster than MATLAB's sgolay.
%
% Usage:
%       [fc,df] = savitzkyGolay(x,n,dn,x0,flag)
%   input:
%       x    = the original data point, e.g., -5:5 
%       n    = polynomial order
%       dn   = differentation order (0=smoothing),  default=0
%       x0   = estimation point, can be a vector    default=0
%       W    = weight vector, can be empty          
%              must have same length as x0          default=identity
%       flag = numerical(0) or symbolical(1),       default=0
%
%   output:
%       fc   = filter coefficients obtained (B output of sgolay).
%       df   = differentiation filters (G output of sgolay).
% Notes:
% 1.    x can be arbitrary, e.g., odd number or even number, symmetrical or
%       nonsymmetrical, uniformly spaced or nonuniformly spaced, etc.       
% 2.    x0 can be arbitrary, e.g., the initial point, the end point, etc.
% 3.    Either numerical results or symbolical results can be obtained.
% Example:
%       sgsdf([-3:3],2,0,0,[],0)
%       sgsdf([-3:3],2,0,0,[],1)
%       sgsdf([-3:3],2,0,-3,[],1)
%       sgsdf([-3:3],2,1,2,[],1)
%       sgsdf([-2:3],2,1,1/2,[],1)
%       sgsdf([-5:2:5],2,1,0,[],1)     
%       sgsdf([-1:1 2:2:8],2,0,0,[],1)
% Author:
%       Diederick C. Niehorster <dcniehorster@hku.hk> 2011-02-05
%       Department of Psychology, The University of Hong Kong
%
%       Originally based on
%       http://www.mathworks.in/matlabcentral/fileexchange/4038-savitzky-golay-smoothing-and-differentiation-filter
%       Allthough I have replaced almost all the code (partially based on
%       the comments on the FEX submission), increasing its compatibility
%       with MATLABs sgolay (now supports a weight matrix), its numerical
%       stability and it speed. Now, the help is pretty much all that
%       remains.
%       Jianwen Luo <luojw@bme.tsinghua.edu.cn, luojw@ieee.org> 2003-10-05
%       Department of Biomedical Engineering, Department of Electrical Engineering
%       Tsinghua University, Beijing 100084, P. R. China  
% Reference
%[1]A. Savitzky and M. J. E. Golay, "Smoothing and Differentiation of Data
%   by Simplified Least Squares Procedures," Analytical Chemistry, vol. 36,
%   pp. 1627-1639, 1964.
%[2]J. Steinier, Y. Termonia, and J. Deltour, "Comments on Smoothing and
%   Differentiation of Data by Simplified Least Square Procedures,"
%   Analytical Chemistry, vol. 44, pp. 1906-1909, 1972.
%[3]H. H. Madden, "Comments on Savitzky-Golay Convolution Method for
%   Least-Squares Fit Smoothing and Differentiation of Digital Data,"
%   Analytical Chemistry, vol. 50, pp. 1383-1386, 1978.
%[4]R. A. Leach, C. A. Carter, and J. M. Harris, "Least-Squares Polynomial
%   Filters for Initial Point and Slope Estimation," Analytical Chemistry,
%   vol. 56, pp. 2304-2307, 1984.
%[5]P. A. Baedecker, "Comments on Least-Square Polynomial Filters for
%   Initial Point and Slope Estimation," Analytical Chemistry, vol. 57, pp.
%   1477-1479, 1985.
%[6]P. A. Gorry, "General Least-Squares Smoothing and Differentiation by
%   the Convolution (Savitzky-Golay) Method," Analytical Chemistry, vol.
%   62, pp. 570-573, 1990.
%[7]Luo J W, Ying K, He P, Bai J. Properties of Savitzky-Golay Digital
%   Differentiators, Digital Signal Processing, 2005, 15(2): 122-136.
%
%See also:
%       sgolay, savitzkyGolayFilt

% Check if the input arguments are valid and apply defaults
error(nargchk(2,6,nargin,'struct'));

if round(n) ~= n, error(generatemsgid('MustBeInteger'),'Polynomial order (n) must be an integer.'), end
if round(dn) ~= dn, error(generatemsgid('MustBeInteger'),'Differentiation order (dn) must be an integer.'), end
if n > length(x)-1, error(generatemsgid('InvalidRange'),'The Polynomial Order must be less than the frame length.'), end
if dn > n, error(generatemsgid('InvalidRange'),'The Differentiation order must be less than or equal to the Polynomial order.'), end

% set defaults if needed
if nargin<6
    flag=false;
end
if nargin < 5 || isempty(W)
   % No weighting matrix, make W an identity
   W = eye(length(x0));
else
   % Check W is real.
   if ~isreal(W), error(generatemsgid('NotReal'),'The weight vector must be real.'),end
   % Check for right length of W
   if length(W) ~= length(x0), error(generatemsgid('InvalidDimensions'),'The weight vector must be of the same length as the frame length.'),end
   % Check to see if all elements are positive
   if min(W) <= 0, error(generatemsgid('InvalidRange'),'All the elements of the weight vector must be greater than zero.'), end
   % Diagonalize the vector to form the weighting matrix
   W = diag(W);
end
if nargin<4
    x0=0;
end
if nargin<3
    dn=0;
end

% prepare for symbolic output
if flag
    x=sym(x);
    x0=sym(x0);
end

Nx  = length(x);
x=x(:);
Nx0 = length(x0);
x0=x0(:);

if flag
    A=ones(length(x),1);
    for k=1:n
        A=[A x.^k];
    end
    df = inv(A'*A)*A';                          % backslash operator doesn't work as expected with symbolic inputs, but the "slowness and inaccuracy" of this method doesn't matter when doing the symbolic version
else
    df = cumprod([ones(Nx,1) x*ones(1,n)],2) \ eye(Nx);
end
df = df.';

hx = [(zeros(Nx0,dn)) ones(Nx0,1)*prod(1:dn)];  % order=0:dn-1,& dn,respectively
for k=1:n-dn                                    % order=dn+1:n=dn+k
    hx = [hx x0.^k*prod(dn+k:-1:k+1)];
end

% filter coeffs
fc = df*hx'*W;