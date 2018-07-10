function [pY,PowerSpectrum,maxy,miny,area,stdev]=isignal(datamatrix,xcenter,xrange,sm,sw,em,dm,rm,s1,s2,sr,mw,spm)
% Y=isignal(DataMatrix,xcenter,xrange,SmoothMode,SmoothWidth,ends,...
% DerivativeMode,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth,SpectrumMode)
% An interactive function that performs smoothing, differentiation, and
% peak sharpening of a time-series signal in the form of a 2-column 
% matrix with the independent variable (x-vfixedparametersalues) in the first
% column and dependent variable (y values) in the second column, or as
% separate x and y vectors. Returns the processed independent axis (Y)
% vector as the output argument. The lower half of the figure window shows
% a plot of the entire signal, and the upper half shows a selected portion
% controlled by the pan and zoom keystrokes or by optional input  
% arguments 'xcenter' and 'xrange', respeSctively. Other keystrokes
% also allow you to control the smooth type, width, and ends
% treatment, the derivative order (0th through 5th), and peak 
% sharpening. (Alternatively, the initial values of these parameters 
% can be passed to the function via the optional input arguments.) 
%
% Version 5.71  Adds Shift-L to replace signal with processed version
% Shift-V for Fourier convolution/Deconvolution function menu
% 
% By T. C. O'Haver (toh@umd.edu); Savitzky-Golay smooth code by Diederick.  
% See http://terpconnect.umd.edu/~toh/spectrum/iSignal.html
% The S key (or optional argument "sm") determines the smooth mode:
%     If sm=0, the signal is not smoothed.
%     If sm=1, rectangular (sliding-average or boxcar) 
%     If sm=2, triangular (2 passes of sliding-average)
%     If sm=3, pseudo-Gaussian (3 passes of sliding-average)
%     If sm=4, Savitzky-Golay smooth 
% The A and Z keys (or optional argument sw) control the smooth width.
% The Z key (or argument "em") controls how the "ends" of the signal 
%   (the first w/2 points and the last w/2 points) are handled.
%     If ends=0, the ends are zeroed
%     If ends=1, the ends are smoothed with progressively 
%     smaller smooths the closer to the end.
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html
%
% The D key (or optional input argument "dm") determines the derivative 
%   order (O, 1, 2, 3, 4, 5, and back to 0). See
%   http://terpconnect.umd.edu/~toh/spectrum/Differentiation.html
%
% The E key (or optional argument "rm") turns off and on peak 
%  sharpening (resolution enhancement). The sharpening strength is
%  controled by the F and V keys (optional argument "s1") and B and G 
%  keys (optional argument "s2"). The optimum values depend on the 
%  peak shape and width; For details, see
%  http://terpconnect.umd.edu/~toh/spectrum/InteractiveResEnhance.htm). 
%
% Shift-S key toggles on and off Frequency Spectrum mode, which computes
%  frequency spectrum of the segment of the signal displayed in the upper
%  window and displays it in the lower window (in red). Use the pan and zoom
%  keys to adjust the region of the signal to be viewed. (Press Ctrl-A to
%  select the entire signal). In the frequency spectrum mode, you can press
%  Shift-A to cycle through four plot modes (linear, semilog X, semilog Y,
%  or log-log) and press Shift-X to toggle between a frequency on the x axis
%  and time on the x-axis. Press Shift-S again to return to the normal mode.
%  Shift-Z toggles on and off peak detection and labeling on the 
%  frequency/time spectrum. Adjust peak detection in lines 2196-2199; see
%  http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
%
% The P key toggles off and on the peak measure mode, which measures and
%  displays the peak position, height, width, and area of the one peak
%  at a time if it is centered and zoomed in; a red "cap" on the peak
%  indicates that portion of the signal that is taken for the measurement.
% The 'R' key prints out the peak measures in the command window.
% The L key toggles off and on the Overlay mode, which overlays the
%  selected portion in the upper plot with the original signal as a dotted
%  line, for comparison. 
% The H key switches between linear and log y-axis on the lower plot. The 0
%  (zero) key set minimun signal to zero.
% The ; key (semicolon) sets the entire selected region to zero (use to
%  remove stray data points). 
% The Tab key resets smooth, derivative, and sharpen effects to zero. 
% The O (letter O) key saves the X,Y processed signal as a "mat" file, in a
%  location and with a file name that you specify. 
% The C key condenses the signal by the specified factor N, replacing each
%  group of N points with their average; 
% The I key replaces the signal with a linearily interploated version
%  containing N data points.
% The M key implements a median filter for removing spikes. The ~ key
%  limits the maximum rate of change of the signal.
% The + key computes the absolute value of the entire signal.
%
% Press K to see all keyboard commands.
%
% EXAMPLE 1: Data in two columns of a matrix: [x y].
%             >> load data.mat
%             >> isignal(DataMatrix);
% 
% EXAMPLE 2: Data in separate x,y vectors or single y vector
%             >> isignal(x,y);  or
%             >> isignal(y);  
%
% EXAMPLE 3: As above, but specifies initial values of pan (xcenter) and 
%            zoom (xrange) in the last two input arguments. 
%             >> isignal(DataMatrix,180,40); or
%             >> isignal(x,y,180,40);
%
% EXAMPLE 4: As above, but additionally specifies initial values of 
%            SmoothMode, SmoothWidth, ends, and DerivativeMode. 
%             >> isignal(DataMatrix,180,40,2,9,0,1);
% 
% EXAMPLE 5: As above, but additionally specifies initial values of the  
%            peak sharpening parameters Sharpen, Sharp1, and Sharp2.
%             >> isignal(DataMatrix,180,40,4,19,0,0,1,51,6000);
%                (Press 'E' key to toggle sharpening ON/OFF)
%
% EXAMPLE 6:   >> x=[0:.005:2];y=humps(x);Data=[x;y];
%             4th derivative of the peak at x=0.9:
%              >> isignal(Data,0.9,0.5,1,3,1,4);
%             Peak sharpening applied to the peak at x=0.3:
%              >> isignal(Data,0.3,0.5,4,3,1,0,1,220,5400);
%                 (Press 'E' key to toggle sharpening ON/OFF)
%
% EXAMPLE 7: Measurement of peak area.  This example generates four 
% Gaussian peaks, all with the exact same peak height (1.00) and area 
% (1.77). The first peak (at x=4) is isolated, the second peak (x=9) 
% is slightly overlapped with the third one, and the last two peaks 
% (at x= 13 and 15) are strongly overlapped.  To measure the area under 
% a peak using the perpendicular drop method, position the dotted red
% marker lines at the minimum between the overlapped peaks.  
% 
% >> x=[0:.01:20];y=exp(-(x-4).^2)+exp(-(x-9).^2)+exp(-(x-13).^2)+exp(-(x-15).^2);
% >> isignal(x,y);
%
% EXAMPLE 8: Single peak with random spikes. Compare smoothing vs spike
%            filter (M key) and slew rate limit (~ key) to remove spikes.
% >> x=-5:.01:5;
% >> y=exp(-(x).^2);for n=1:1000,if randn()>2,y(n)=rand()+y(n);,end,end;
% >> isignal(x,y);
%
% Example 9: Weak peak at x=128 on a smooth, curved background. 
% Try second derivative + smoothing
% >> x=1:.1:256;
% >> y=gaussian(x,-100,300)+.02.*gaussian(x,128,30)+0.001.*randn(size(x));
% >> isignal(x,y);
%  
% Example 10: Spectrum mode
% >> x=0:.1:60; y=sin(x)+sin(10.*x);
% >> [pY,PowerSpectrum]=isignal([x;y],30,30,4,3,1,0,0,1,0,0,0,1);
% >> plot(PowerSpectrum)
%
% Example 11: Noisy 4th derivative signal. Adjust smoothing to reveal peak
% at x=150, height=1e-4; SNR=90.
% isignal(x,deriv4(100000.*gaussian(x,150,PeakWidth)+.1*randn(size(x))));
%
% KEYBOARD CONTROLS, version 5.5:
%  Pan signal left and right...Coarse pan: < and >
%                              Fine pan: left and right cursor arrows
%                              Nudge: [ and ]
%  Zoom in and out.............Coarse zoom: / and "  
%                              Fine zoom: up and down cursor arrows
%  Resets pan and zoom.........ESC
%  Select entire signal........Ctrl-A
%  Display Grid (on/off).......Shift-G  Temporarily displays grid on plots
%  Adjust smooth width.........A,Z  (A=>more, Z=>less) 
%  Adjust smooth type..........S (No, Rectanglular, Triangle, Gaussian, Savitzky-Golay)
%  Toggle smooth ends..........X (0=ends zeroed  1=ends smoothed (slower)
%  Cycle derivative orders.....D/Shift-D Increase/Decrease derivative order
%  Toggle peak sharpening......E (0=OFF 1=ON)
%  Sharpening for Gaussian.....Y  Set sharpen settings for Gaussian
%  Sharpening for Lorentzian...U  Set sharpen settings for Lorentzian
%  Adjust sharp1...............F,V  F=>sharper, V=>less sharpening
%  Adjust sharp2...............G,B  G=>sharper, B=>less sharpening
%  Slew rate limit (0=OFF).....`  Largest allowed change between points
%  Spike filter width (O=OFF)..m  Spike filter eliminates sharp spikes
%  Toggle peak parabola........P  fits parabola to center, labels vertex
%  Fit polynomial to segment...Shift-o  Asks for polynomial order
%  Fits peak in upper window...Shift-F (Asks for shape, number of peaks, etc)
%  Spectrum mode on/off........Shift-S (Shift-A and Shift-X to change axes)
%  Peak labels on spectrum.....Shift-Z in spectrum/time mode
%  Display Waterfall spectrum..Shift-W  Allows choice of mesh, surf, contour, or pcolor
%  Lock in current processing..Shift-L  Replace signal with processed version
%  ConVolution/DeconVolution...Shift-V  Convolution/Deconvolution menu
%  Click on graph..............C  Prints out x and y coordinates
%  Print peak report...........R  prints position, height, width, area
%  Toggle overlay mode.........L  Overlays original signal as dotted line
%  Toggle log y mode...........H  semilog plot in lower window
%  Cycles baseline mode........T  none, linear, quadratic, or flat baseline mode
%  Restores original signal....Tab or Ctrl-Z key resets to previous signal
%  Baseline subtraction........Backspace, then click baseline at multiple points
%  Restore background..........\  to cancel previous background subtraction
%  Invert signal...............Shift-N  Invert (negate) the signal (flip + and -)
%  Remove offset...............0  (zero) set minimun signal to zero 
%  Sets region to zero.........;  sets selected region to zero.
%  Absolute value..............+  Computes absolute value of entire signal')
%  Condense signal.............C  Condense oversampled signal by factor N
%  Interpolate signal..........i  Interpolate (resample) to N points
%  Print keyboard commands.....K  prints this list
%  Print signal report.........Q  prints signal info and current settings
%  Print isignal arguments.....W  prints isignal (current arguments)
%  Save output to disk.........O as .mat file with processed signal matrix
%  Play signal as sound........Spacebar or Shift-P Play selected section
%  Play signal as sound........Shift-R Change sampling rate for playing sound
%  Switch to ipf.m.............Shift-Ctrl-F transfer current signal to ipf.m
%  Switch to iPeak.............Shift-Ctrl-P transfer current signal to iPeak.m

% Copyright (c) 2016, Thomas C. O'Haver
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
% 
global X Y xo dx DerivativeMode Sharpen Sharp1 Sharp2 SmoothWidth SlewRate MedianWidth
global SmoothType ends SmoothMode SavedSignal SavedXvalues PeakLabels Report autozero
global plotmode xmode SpectrumMode samplerate LabelPeaks NumSigs DataMatrix
format short g
format compact
warning off all
DataMatrix=datamatrix;
NumSigs=1;
switch nargin % Process arguments
    % 'nargin' is the number of arguments
    case 1  % One argument only
        % Might be isignal(DataMatrix) ot isignal(Y-vector)
        % If data is in the wrong transposition, fix it.
        datasize=size(DataMatrix);  
        if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
        datasize=size(DataMatrix);
        if datasize(2)==1, %  Must be isignal(Y-vector)
            X=1:length(DataMatrix); % Create an independent variable vector
            Y=DataMatrix;
        else
            % Must be isignal(DataMatrix)
            X=DataMatrix(:,1); % Split matrix argument
            Y=DataMatrix(:,2);
        end
        SmoothMode=0; % Initial SmoothMode = Rect
        SmoothWidth=3; % Initial smooth width
        SmoothType='Rect.';  % Label for initial SmoothMode
        DerivativeMode=0; % Derivative mode off initially
        ends=0;  % Initial smooth ends setting zero (no tapering)
        Sharpen=0;  % Initially no sharpening
        Sharp1=10; % Initial factor1 for resolution enhancement
        Sharp2=1000; % Initial factor2 for resolution enhancement
        SlewRate=0;
        MedianWidth=0;
        xo=length(Y)/2; % Initial Pan setting
        dx=length(Y); % Initial Zoom setting
        SpectrumMode=0; % Frequency spectrum initially off. 
    case 2
        % Two arguments, might be separate x and y data vectors,
        % or one data matrix and a peak density estimate.
        if isscalar(xcenter) % if second argument is scalar
            % Must be isignal(DataMatrix,xcenter)
            % If DataMatrix is in the wrong transposition, fix it.
            datasize=size(DataMatrix);
            if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
            X=DataMatrix(:,1); % Split matrix argument
            Y=DataMatrix(:,2);
            xo=val2ind(X,xcenter);
        else % if second argument is not scalar
            % Must be isignal(x,y)
            xdatasize=size(DataMatrix);
            if xdatasize(1)<xdatasize(2),DataMatrix=DataMatrix';end
            X=DataMatrix;  % First argument is X
            xdatasize=size(xcenter);
            if xdatasize(1)<xdatasize(2),xcenter=xcenter';end
            Y=xcenter; % Second argument is Y
            xo=length(Y)/2; %  % Default initial zoom setting
        end  % if isscalar
        SmoothMode=0; % Initial SmoothMode = No
        SmoothWidth=3; % Initial smooth width
        SmoothType='Rect.';  % Label for initial SmoothMode
        DerivativeMode=0; % Derivative mode off initially
        ends=0;  % Initial smooth ends setting zero (no tapering)
        Sharpen=0;  % Initially no sharpening
        Sharp1=10; % factor1 for resolution enhancement
        Sharp2=1000; % factor2 for resolution enhancement
        SlewRate=0;
        MedianWidth=0;
        dx=length(Y); %  % Default initial zoom setting
        SpectrumMode=0; % Frequency spectrum initially off. 
    case 3
        % Might be isignal(DataMatrix,xcenter,xrange) or isignal(x,y,xcenter)
        if isscalar(xcenter) % if second argument is scalar
            % Must be isignal(DataMatrix,xcenter,xrange)
            % If DataMatrix is in the wrong transposition, fix it.
            datasize=size(DataMatrix);
            if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
            X=DataMatrix(:,1); % Split matrix argument
            Y=DataMatrix(:,2);
            [xo,dx]=panandzoom(X,xcenter,xrange);
        else % if second argument is not isscalar
            % Must be isignal(x,y,xcenter)
            xdatasize=size(DataMatrix);
            if xdatasize(1)<xdatasize(2),DataMatrix=DataMatrix';end
            X=DataMatrix;  % First argument is X
            xdatasize=size(xcenter);
            if xdatasize(1)<xdatasize(2),xcenter=xcenter';end
            Y=xcenter; % Second argument is Y
            xo=xrange; % third argument is xcenter
            dx=length(Y); % Default initial zoom setting
        end  % if isscalar
        SmoothMode=0; % Initial SmoothMode = No
        SmoothWidth=3; % Initial smooth width
        SmoothType='Rect.';  % Label for initial SmoothMode
        DerivativeMode=0; % Derivative mode off initially
        ends=0;  % Initial smooth ends setting zero (no tapering)
        Sharpen=0;  % Initially no sharpening
        Sharp1=10; % factor1 for resolution enhancement
        Sharp2=1000; % factor2 for resolution enhancement
        SlewRate=0;
        MedianWidth=0;
        SpectrumMode=0; % Frequency spectrum initially off. 
    case 4   % Must be isignal(x,y,xcenter,xrange)
        xdatasize=size(DataMatrix);
        if xdatasize(1)<xdatasize(2),DataMatrix=DataMatrix';end
        X=DataMatrix;  % First argument is X
        xdatasize=size(xcenter);
        if xdatasize(1)<xdatasize(2),xcenter=xcenter';end
        Y=xcenter; % Second argument is Y
        SmoothMode=0; % Initial SmoothMode = No
        SmoothWidth=3; % Initial smooth width
        SmoothType='Rect.';  % Label for initial SmoothMode
        DerivativeMode=0; % Derivative mode off initially
        ends=0;  % Initial smooth ends setting zero (no tapering)
        Sharpen=0;  % Initially no sharpening
        Sharp1=10; % factor1 for resolution enhancement
        Sharp2=1000; % factor2 for resolution enhancement
        SlewRate=0;
        MedianWidth=0;
        [xo,dx]=panandzoom(X,xrange,sm);
        SpectrumMode=0; % Frequency spectrum initially off. 
    case 7
        % One data matrix, all smoothing and derivative parameters specified
        % in arguments, default values for resolution enhancement.
        % If DataMatrix is in the wrong transposition, fix it.
        datasize=size(DataMatrix);
        if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
        X=DataMatrix(:,1); % Split matrix argument
        Y=DataMatrix(:,2);
        [xo,dx]=panandzoom(X,xcenter,xrange);
        SmoothMode=sm; % SmoothMode (O, 1, 2, 3 or 4)
        SmoothWidth=sw; % Smooth width
        DerivativeMode=dm; % Derivative mode (0, 1, 2, 3, 4)
        ends=em;  % Smooth ends setting (0 or 1)
        Sharpen=0;  % Initially no sharpening
        Sharp1=10; % factor1 for resolution enhancement
        Sharp2=1000; % factor2 for resolution enhancement
        SlewRate=0;
        MedianWidth=0;
        SpectrumMode=0; % Frequency spectrum initially off. 
    case 10
        % One data matrix, all signal processing parameters specified
        % in arguments, including  resolution enhancement.
        % If DataMatrix is in the wrong transposition, fix it.
        datasize=size(DataMatrix);
        if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
        X=DataMatrix(:,1); % Split matrix argument
        Y=DataMatrix(:,2);
        [xo,dx]=panandzoom(X,xcenter,xrange);
        SmoothMode=sm; % SmoothMode sm
        SmoothWidth=sw; % Smooth width
        DerivativeMode=dm; % Derivative mode (0, 1, 2, 3, 4)
        ends=em;  % Smooth ends setting (0 or 1)
        Sharpen=rm;  % Sharpen mode
        Sharp1=s1; % factor1 for resolution enhancement
        Sharp2=s2; % factor2 for resolution enhancement
        SlewRate=0;
        MedianWidth=0;
        SpectrumMode=0; % Frequency spectrum initially off. 
    case 11
        % One data matrix, all signal processing parameters specified
        % in arguments, except MedianWidth.
        % If DataMatrix is in the wrong transposition, fix it.
        datasize=size(DataMatrix);
        if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
        X=DataMatrix(:,1); % Split matrix argument
        Y=DataMatrix(:,2);
        [xo,dx]=panandzoom(X,xcenter,xrange);
        SmoothMode=sm; % SmoothMode sm
        SmoothWidth=sw; % Smooth width
        DerivativeMode=dm; % Derivative mode (0, 1, 2, 3, 4)
        ends=em;  % Smooth ends setting (0 or 1)
        Sharpen=rm;  % Sharpen mode
        Sharp1=s1; % factor1 for resolution enhancement
        Sharp2=s2; % factor2 for resolution enhancement
        SlewRate=sr;
        MedianWidth=0;
        SpectrumMode=0; % Frequency spectrum initially off. 
    case 12
        % One data matrix, all signal processing parameters specified
        % in arguments, except SpectrumMode.
        % If DataMatrix is in the wrong transposition, fix it.
        datasize=size(DataMatrix);
        if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
        X=DataMatrix(:,1); % Split matrix argument
        Y=DataMatrix(:,2);
        [xo,dx]=panandzoom(X,xcenter,xrange);
        SmoothMode=sm; % SmoothMode sm
        SmoothWidth=sw; % Smooth width
        DerivativeMode=dm; % Derivative mode (0, 1, 2, 3, 4)
        ends=em;  % Smooth ends setting (0 or 1)
        Sharpen=rm;  % Sharpen mode
        Sharp1=s1; % factor1 for resolution enhancement
        Sharp2=s2; % factor2 for resolution enhancement
        SlewRate=sr;
        MedianWidth=mw;
        SpectrumMode=0; % Frequency spectrum initially off. 
     case 13
        % One data matrix, all signal processing parameters specified
        % in arguments, except SpectrumMode.
        % If DataMatrix is in the wrong transposition, fix it.
        datasize=size(DataMatrix);
        if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
        X=DataMatrix(:,1); % Split matrix argument
        Y=DataMatrix(:,2);
        [xo,dx]=panandzoom(X,xcenter,xrange);
        SmoothMode=sm; % SmoothMode sm
        SmoothWidth=sw; % Smooth width
        DerivativeMode=dm; % Derivative mode (0, 1, 2, 3, 4)
        ends=em;  % Smooth ends setting (0 or 1)
        Sharpen=rm;  % Sharpen mode
        Sharp1=s1; % factor1 for resolution enhancement
        Sharp2=s2; % factor2 for resolution enhancement
        SlewRate=sr;
        MedianWidth=mw;
        SpectrumMode=spm;
    otherwise
        disp('Invalid number of arguments')
        disp('Expected forms are:')
        disp('isignal(y);  % Data in single y vector')
        disp('isignal(x,y);  % Data in separate x and y vectors')
        disp('isignal(DataMatrix); % Data in two columns of DataMatrix')
        disp('isignal(x,y,xcenter,xrange); ')
        disp('isignal(DataMatrix,xcenter,xrange;) ')
        disp('isignal(DataMatrix,xcenter,xrange,SmoothMode,SmoothWidth,ends,DerivativeMode); ' )
        disp('isignal(DataMatrix,xcenter,xrange,SmoothMode,SmoothWidth,ends,DerivativeMode,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth); ')
        beep
        return
end % switch nargin
% Define smooth type string for xlabel
switch SmoothMode
     case 0
          SmoothType='No';
     case 1
          SmoothType='Rect.';
     case 2
          SmoothType='Tri.';
     case 3
          SmoothType='Gauss';
     case 4
           SmoothType='Savitzky-Golay';
end
PeakLabels=0;  % Start with peak label turned off
% Save original signal in SavedSignal for undo function 
SavedSignal=Y;
SavedXvalues=X;
Overlay=0;  % Start with overlay turned off
Report=0;
autozero=0;
logymode=0; % Start in linear y mode.
plotmode=2; % Frequency spectrum initially in semilog y mode.
% NumPeaksUW=1;
xmode=0;% Frequency spectrum initially in frequency mode.
PowerSpectrum=0;
samplerate=44100;
LabelPeaks=0;
extra=0;
NumTrials=1;
% Plot the signal
[xx,yy]=RedrawSignal(X,Y,xo,dx);
pY=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
Y=pY;
% xo=length(Y)/2; % Initial Pan setting           
% dx=length(Y); % Initial Zoom setting
[xx,yy]=RedrawSignal(X,Y,xo,dx);
maxy=max(yy);
miny=min(yy);
area=trapz(xx,yy);
stdev=std(yy);
if SpectrumMode==1;
    % Plot the power spectrum  in the lower half of the window.
     [f realsy PowerSpectrum]=PlotFrequencySpectrum(X,Y,xo,dx,plotmode,xmode,LabelPeaks);
    subplot(2,1,1)
    title('iSignal 5   Frequency Spectrum Mode (Press Shift-S again to cancel')
end

% Attaches KeyPress test function to the figure.
set(gcf,'KeyPressFcn',@ReadKey)
uicontrol('Style','text')
% end of outer function
% ----------------------------SUBFUNCTIONS--------------------------------
function ReadKey(obj,eventdata)
% Interprets key presses from the Figure window. When a key is pressed,
% executes the code in the corresponding section in the SWITCH statement.
% Note: If you don't like my key assignments, you can change the numbers
% in the case statements here to re-assign that function to any other key.
% If you press a key that has not yet beZen assigned to a function, it
% displays the key code number in the command window so you can easily
% add that to the SWITCH statement for your own custom functions.
global X Y xx yy xo dx SmoothMode SmoothWidth DerivativeMode SlewRate extra SpectrumMode
global Sharpen Sharp1 Sharp2 SavedSignal SavedXvalues SavedBackground plotmode
global SmoothType ends PeakLabels Overlay Report GaussEstimate MedianWidth xmode
global LorentzEstimate logymode autozero Shape NumTrials NumPeaksUW fixedparameters
global PowerSpectrum samplerate LabelPeaks SS NumSigs DataMatrix AllMode SavedMatrix
key=get(gcf,'CurrentCharacter');
if isscalar(key),
    ly=length(Y);
    switch double(key),
        case 28
            % Pans down when right arrow pressed.
            xo=xo-dx/20;
            if xo<1,xo=1;end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 29
            % Pans up when left arrow pressed.
            xo=xo+dx/20;
            if xo>ly,xo=ly;end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 91
            % Nudge down 1 point when [ pressed.
            xo=xo-1;
            if xo<1,xo=1;end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 93
            % Nudge up 1 point when ] pressed.
            xo=xo+1;
            if xo>ly,xo=ly;end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 44
            % Pans down when > key pressed.
            xo=xo-dx;
            if xo<1,xo=1;end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 46
            % Pans up when < key pressed.
            xo=xo+dx;
            if xo>ly,xo=ly;end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 31
            % Zooms up when up arrow pressed.
            dx=dx+dx/10;
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 30
            % Zooms down when down arrow pressed.
            dx=dx-dx/10;
            if dx>ly,dx=ly;end
            if dx<2,dx=2;end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 47
            % Zooms x 2 up when / pressed.
            dx=dx*2;
            if dx>ly,dx=ly;end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 39
            % Zooms x 1/2 down when ' pressed.
            dx=round(dx/2);
            if dx<2,dx=2;end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 1 % Ctrl-A selects entire signal
            xo=length(Y)/2; % Initial Pan setting           
            dx=length(Y); % Initial Zoom setting
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 71 % Shift-G temporatiy displays grid on both upper and lower panels
            subplot(211);grid
            subplot(212);grid
%         case {49,50,51,52,53,54,55,56,57}
%             % When a number key is pressed, sets the signal
%             % number
%             % NumSigs=NumSigs; % Testing
%             AllMode=0;
%             SS=key-48;
%             if SS>NumSigs,
%                 disp('The signal matrix is not that large.')
%                 SS=NumSigs;
%             end
%             datasize=size(DataMatrix);
%             % if SS>sizey(2),SS=sizey(2);end
%             % SignalSelected=[SS NumSigs]
%             %             sizeX=size(X)
%             %             sizeY=size(Y)
%             switch SmoothMode(SS)
%                 case 0
%                     SmoothType='No';
%                 case 1
%                     SmoothType='Rect.';
%                 case 2
%                     SmoothType='Tri.';
%                 case 3
%                     SmoothType='Gauss';
%                 case 4
%                     SmoothType='Savitzky-Golay';
%             end
%             Y=ProcessSignal(X,DataMatrix(:,SS),DerivativeMode(SS),SmoothWidth(SS),SmoothMode(SS),ends(SS),Sharpen(SS),Sharp1(SS),Sharp2(SS),SlewRate(SS),MedianWidth(SS));
%             [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 41  % When Shift-0 ')' key is pressed, asks for the signal number
            SS=input('Enter desired signal number and press Enter:');
            if isempty(SS),SS=1;end
            if isnan(SS),SS=1;end
            if SS>NumSigs,
                disp('The signal matrix is not that large.')
                SS=NumSigs;
            end
            switch SmoothMode(SS)
                case 0
                    SmoothType='No';
                case 1
                    SmoothType='Rect.';
                case 2
                    SmoothType='Tri.';
                case 3
                    SmoothType='Gauss';
                case 4
                    SmoothType='Savitzky-Golay';
            end
            Y=DataMatrix(:,SS);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 63 % Shift-?
            disp(['SizeDataMatrix= ' num2str(size(DataMatrix)) ] )
            MaxY=max(max(Y));disp(['MaxY = ' num2str(MaxY) ] )
            MinY=min(min(Y));disp(['MinY = ' num2str(MinY) ] )
            disp(['DerivativeMode = ' num2str(DerivativeMode) ] )
            disp(['SmoothWidth = ' num2str(SmoothWidth) ] )
            disp(['SmoothMode = ' num2str(SmoothMode) ] )
            disp(['Sharpen = ' num2str(Sharpen) ] )
            disp(['Sharp1 = ' num2str(Sharp1) ] )
            disp(['Sharp2 = ' num2str(Sharp2) ] )
            disp(['MedianWidth = ' num2str(MedianWidth) ] )
            disp(['SlewRate = ' num2str(SlewRate) ] )
            disp(['Peak Shape = ' num2str(Shape) ] )
            disp(['NumPeaks = ' num2str(NumPeaksUW) ] )
            disp(['NumTrials = ' num2str(NumTrials) ] )  
            disp(['Baseline mode = ' num2str(autozero) ] )  
        case 77 % Shift-M key applies current processing variable to all signals
            for sig=1:NumSigs
                DerivativeMode(sig)=DerivativeMode(SS);
                SmoothWidth(sig)=SmoothWidth(SS);
                SmoothMode(sig)=SmoothMode(SS);
                Sharpen(sig)=Sharpen(SS);
                ends(sig)=ends(SS);
                Sharp1(sig)=Sharp1(SS);
                Sharp2(sig)=Sharp2(SS);
                MedianWidth(sig)=MedianWidth(SS);
                Y=ProcessSignal(X,SavedMatrix(:,sig),DerivativeMode(SS),SmoothWidth(SS),SmoothMode(SS),ends(SS),Sharpen(SS),Sharp1(SS),Sharp2(SS),SlewRate(SS),MedianWidth(SS));
                DataMatrix(:,sig)=Y;
            end
            disp('Current signal processing applied to all signals')
            [xx,yy]=RedrawSignal(X,Y,xo,dx);   
        case 8
            % When 'Backspace' key is pressed, user clicks the graph
            % along the presumed background points, then the program
            % subtracts the interploated background between the points.
            SavedBackground=Y;
            disp('Multi-point baseline subtraction')
            BaselinePoints=input('Number of baseline points to click): ');
            if isempty(BaselinePoints),BaselinePoints=8;end
            % Acquire background points from user mouse clicks
            subplot(2,1,2)
            title(['Click on ' num2str(BaselinePoints) ' points on the baseline between the peaks.'])
            bX=[];bY=[];
            for g=1:BaselinePoints;
                [clickX,clickY] = ginput(1);
                bX(g)=clickX;
                bY(g)=clickY;
                xlabel(['Baseline point '  num2str(g) ' / ' num2str(BaselinePoints) ])
            end
            tempY=Y;
            for k=1:length(bX)-1,
                fp=val2ind(X,bX(k)); % First point in segment
                lp=val2ind(X,bX(k+1));  % Last point in segment
                % Subtract piecewise linear background from Y
                tempY(fp:lp)=Y(fp:lp)-((bY(k+1)-bY(k))/(bX(k+1)-bX(k))*(X(fp:lp)-bX(k))+bY(k));
            end
            Y=tempY;
            SavedSignal=Y;
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 92
            % When '\' key is pressed, restoreed original signal
            SavedSignal=SavedBackground;
            X=SavedXvalues;
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 122
            % When 'z' key is pressed, DEcreases "SmoothWidth" by 1 or 10%
            if SmoothMode==0,
                SmoothMode=1;
                SmoothType='Rect.';
            end
            if SmoothWidth>20,
                SmoothWidth=round(SmoothWidth-.1.*SmoothWidth);
                SmoothWidth=2*round(SmoothWidth/2)-1;
            else
                SmoothWidth=SmoothWidth-1;
            end
            if SmoothWidth<1, SmoothWidth=1;end
            if SmoothMode==0,
                Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            else
                Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 97
            % When 'a' key is pressed, INcreases "SmoothWidth" by 1 or 10%
            
            if SmoothMode==0,
                SmoothMode=1;
                SmoothType='Rect.';
            end
            if SmoothWidth>20,
                SmoothWidth=round(SmoothWidth+.1.*SmoothWidth);
                SmoothWidth=2*round(SmoothWidth/2)+1;
            else
                SmoothWidth=SmoothWidth+1;
            end
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 115 % When 's' key is pressed, steps through SmoothModes
            SmoothMode=SmoothMode+1;
            if SmoothMode==5,SmoothMode=0; end
            switch SmoothMode
                case 0
                    SmoothType='No';
                case 1
                    SmoothType='Rect.';
                case 2
                    SmoothType='Tri.';
                case 3
                    SmoothType='Gauss';
                case 4
                    SmoothType='Savitzky-Golay';
            end
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 120 % When 'x' key is pressed, toggles between ends 0 and 1
            if ends==0,
                ends=1;
            else
                ends=0;
            end
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 100
            % When 'd' key is pressed, cycles through DerivativeModes 0,1,2,3,4,5->0
            % if length(Y)>10000,disp('Warning: Derivatives can be slow for for signal lengths above 10,000 points'),end
            DerivativeMode=DerivativeMode+1;
            if DerivativeMode==6,DerivativeMode=0; end
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 68
            % When 'Shift-d' key is pressed, cycles through DerivativeModes 5,4,3,2,1->0
            % if length(Y)>10000,disp('Warning: Derivatives can be slow for for signal lengths above 10,000 points'),end
            DerivativeMode=DerivativeMode-1;
            if DerivativeMode==-1,DerivativeMode=5; end
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 101 % When 'e' key is pressed, toggles between Sharpen 0 and 1
            % if length(Y)>10000,disp('Warning: Sharpening can be slow for for signal lengths above 10,000 points'),end
            if Sharpen==0,
                Sharpen=1;
                SmoothMode=4;
                if SmoothWidth<3;SmoothWidth=3;end
                SmoothType='Savitzky-Golay';
            else
                Sharpen=0;
            end
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 121  % When 'y' key is pressed, sets Sharp1 and 2 for Gaussian
            GaussEstimate=1;
            % if length(Y)>10000,disp('Warning: Sharpening can be slow for
            % for signal lengths above 10,000 points'),end
            PeakLabels=1;
            SmoothMode=4;
            if SmoothWidth<3;SmoothWidth=3;end
            SmoothType='Savitzky-Golay';
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
            Sharpen=1;
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 117  % When 'u' key is pressed, sets Sharp1 and 2 for Lorentzian
            LorentzEstimate=1;
            % if length(Y)>10000,disp('Warning: Sharpening can be slow for for signal lengths above 10,000 points'),end
            SmoothMode=4;
            if SmoothWidth<3;SmoothWidth=3;end
            SmoothType='Savitzky-Golay';
            PeakLabels=1;
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
            Sharpen=1;
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 102
            % When 'f' key is pressed, increases Sharp1
            if Sharpen==0,Sharpen=1;end
            Sharp1=Sharp1+.1*Sharp1;
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 118
            % When 'v' key is pressed, decreases Sharp1
            if Sharpen==0,Sharpen=1;end
            Sharp1=Sharp1-.1*Sharp1;
            if Sharp1<0, Sharp1=0;end
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 103
            % When 'g' key is pressed, increases Sharp2
            if Sharpen==0,Sharpen=1;end
            Sharp2=Sharp2+.1*Sharp2;
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 98
            % When 'b' key is pressed, decreases Sharp2
            if Sharpen==0,Sharpen=1;end
            Sharp2=Sharp2-.1*Sharp2;
            if Sharp2<0, Sharp2=0;end
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 112
            % When 'p' key is pressed, toggles on/off peak labels in upper panel
            if PeakLabels==0,
                PeakLabels=1;
            else
                PeakLabels=0;
            end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 108 % When 'L' key is pressed, toggles between ends 0 and 1
            if Overlay==0,
                Overlay=1;
            else
                Overlay=0;
            end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 104 % When 'H' key is pressed, toggles between normal and logy plot
            if logymode==0,
                logymode=1;
            else
                logymode=0;
            end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 114 % When 'R' key is pressed, toggles between Report 0 and 1
            if Report==0,
                Report=1;
                switch autozero,
                    case 0
                        disp('No baseline correction')
                    case 1
                        disp('Linear baseline subtraction')
                    case 2
                        disp('Quadratic subtraction baseline')
                    case 3
                        disp('Flat baseline correction')
                end
                disp('Position     Height       Width     Gauss. Area   Total Area     SNR       FWHM')
            else
                Report=0;
            end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 116 % When 'T' key is pressed, toggles between normal and autozero plot
            % When 't' key is pressed, steps through AUTOZERO modes
            autozero=autozero+1;
            if autozero==4,autozero=0;end
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case {9,26} % When 'Tab' or Ctrl-Z key is pressed, resets to original signal and modes
            Y=SavedSignal;
            X=SavedXvalues;
            SmoothMode=0;
            SmoothWidth=1;
            SmoothType='No';
            DerivativeMode=0;
            Sharpen=0;
            SlewRate=0;
            MedianWidth=0;
            Y=ProcessSignal(X,Y,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 27 % When 'ESC' key is pressed, resets pan and zoom
            xo=length(Y)/2; % Initial Pan setting
            dx=length(Y)/4; % Initial Zoom setting
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 59
            % When ';' key is pressed, replaces selected segment with zeros
            doit=input('Replace selected region with zeros?','s');
            if doit=='y',
                startpoint=round(xo-dx/2);
                if startpoint<1;startpoint=1;end
                endpoint=round(xo+dx/2);
                if endpoint>length(Y);endpoint=length(Y);end
                Y(startpoint:endpoint)=0;
            end
             [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 48
            % When '0' key is pressed, minimum value of Y is set to zero
            Y=Y-min(Y);
            SavedSignal=Y;
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
            disp('Signal re-zeroed')
        case 78
            % When 'Shift-N' key is pressed, invert the signal
            SavedSignal=-SavedSignal;
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 61
            % When '+' plus key is pressed, compute absolute value of the signal
            SavedSignal=abs(SavedSignal);
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 96
            % When '~' (tilde) key is pressed, enforce maximum slew rate
            SavedSignal=Y;
            disp(['Current slew rate limit =' num2str(SlewRate)])
            SlewRate=input('Enter desired slew rate:');
            if SlewRate=='',SlewRate=0;end
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 109
            % When 'm' key is pressed, performs median filter
            SavedSignal=Y;
            disp(['Current spike width =' num2str(MedianWidth)])
            MedianWidth=input('Enter spike width (1,2,3,...):');
            if MedianWidth=='',MedianWidth=0;end
            MedianWidth=round(MedianWidth);
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 127
            % When 'Delete' key is pressed, sets the single point under the
            % green cursor to zero
            Y(round(xo))=0;
            SavedSignal=Y;
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
%       case 95  % FUTURE ADDITION?
%             % When '_' key (Shift '-') is pressed, replaces selected region
%             startpoint=round(xo-dx/2);
%             if startpoint<1;startpoint=1;end
%             endpoint=round(xo+dx/2)-1;
%             if endpoint>length(Y);endpoint=length(Y);end
%             lxx=length(xx);
%             bkgsize=2;
%             X1=xx(1:round(lxx/bkgsize));
%             X2=xx((lxx-round(lxx/bkgsize)):lxx);
%             MX=[X1;X2];
%             Y1=yy(1:round(length(xx)/bkgsize));
%             Y2=yy((lxx-round(lxx/bkgsize)):lxx);
%             MY=[Y1;Y2];
%             bkgcoef=polyfit(MX,MY,1);  % Fit straight line to sub-group of points
%             bkg=polyval(bkgcoef,xx);
%             Y(startpoint:endpoint)=bkg;
%             Y=ProcessSignal(X,Y,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
%             [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 73
            % When 'I' key (upper-case i or Shift-i) is pressed, integrates the signal
            sum=0;
            for n=1:length(X),
                sum=sum+Y(n);
                Y(n)=sum;
            end
            [xx,yy]=RedrawSignal(X,Y.*(X(2)-X(1)),xo,dx);
        case 6 % Shift-Ctrl-F transfers current signal to Interactive Curve Fitter
            ipf(X,Y);
        case 16 % Shift-Ctrl-P transfers current signal to Interactive Peak Detector
            ipeak(X,Y);
        case 21 % Shift-Ctrl-U transfers current signal to Interactive Fourier Filter
            ifilter(X,Y);
        case 105
            % When 'i' key (lower-case i) is pressed, interpolates the signal
            % to find XI,YI, the values of the underlying function Y at the points
            % linearly interpolated between the points of X, using interp1.
            disp(['X,Y size before interpolation = ' num2str(size(X)) ' , '  num2str(size(Y)) ] )
            InterPoints=input('Number of points in interpolated signal: ');
            if InterPoints>1,
                Xi=linspace(min(X),max(X),InterPoints);
                Y=interp1(X,Y,Xi)';
                X=Xi';
                disp(['X,Y size after interpolation = ' num2str(size(X)) ' , '  num2str(size(Y)) ] )
                xo=length(Y)/2; % Initial Pan setting
                dx=length(Y)/4; % Initial Zoom setting
                SavedSignal=Y;
                SavedXvalues=X;
                pY=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
                Y=pY;
                RedrawSignal(X,Y,xo,dx);
            end
        case 99
            % When C key is pressed, condenses signal by specified factor
            CondenseFactor=input('Condense oversampled signal by factor of (e.g. 2, 3, 4...): ');
            if CondenseFactor>1,
                disp([ 'X,Y size before condensation = ' num2str(size(X)) ' , '  num2str(size(Y)) ] )
                X=condense(X,CondenseFactor)';
                Y=condense(Y,CondenseFactor)';
                xo=length(Y)/2; % Initial Pan setting
                dx=length(Y)/4; % Initial Zoom setting
                SavedSignal=Y;
                SavedXvalues=X;
                disp([ 'X,Y size after condensation = ' num2str(size(X)) ' , '  num2str(size(Y)) ] )
                pY=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
                Y=pY;
                RedrawSignal(X,Y,xo,dx);
            end
        case 113
            % When 'Q' is pressed, prints a report listing signal
            % characteristics and current settings.
            disp('--------------------------------------------------------')
            % SizeX=size(X)
            % SizeY=size(Y)
            disp(['X,Y size = ' num2str(size(X)) ' , '  num2str(size(Y)) ] )
            disp([ num2str(length(Y)) ' total points from x= ' num2str(X(1)) ' to '  num2str(X(length(X))) ] )
            disp(['Interval between points = ' num2str(X(2)-X(1)) ' to ' num2str(X(length(X))-X(length(X)-1)) ] )
            disp(sprintf('min/max Y = %0.3g / %0.4g', min(Y), max(Y)))
                switch autozero,
                    case 0
                        disp('No baseline correction')
                    case 1
                        disp('Linear baseline subtraction')
                    case 2
                        disp('Quadratic subtraction baseline')
                    case 3
                        disp('Flat baseline correction')
                end
            if SlewRate,
                disp(['Maximum slew rate = ' num2str(SlewRate) ] ),
            end
            if MedianWidth,
                disp(['spike filter width = ' num2str(MedianWidth) ] ),
            end
            disp(['Smooth: = ' num2str(SmoothWidth) ' point ' SmoothType ', Ends = ' num2str(ends) ] )
            if DerivativeMode,
                disp(['Derivative order = ' num2str(DerivativeMode) ] ),
            end
            if Sharpen
                disp(['Sharpen factor 1 = ' num2str(Sharp1) '  Sharpen factor 2 = ' num2str(Sharp2) ] )
            end
            disp([ 'Selected range: ' num2str(length(xx)) ' points from x=' num2str(min(xx)) ' to ' num2str(max(xx)) ])
            disp(sprintf('  Peak-to-peak Y: %0.4g  \r  Standard deviation: %0.3g ', max(yy)-min(yy), std(yy)))
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
            disp(sprintf('  Area: %0.4g', trapz(xx,yy)))
            disp(sprintf('  Percent of total area: %0.4g',100*trapz(xx,yy)./trapz(X,Y)))
        case 107
            % When 'k' key is pressed, prints out table of keyboard commands
            disp('iSignal KEYBOARD CONTROLS, version 5.5:')
            disp(' Pan signal left and right...Coarse pan: < and >')
            disp('                             Fine pan: left and right cursor arrows')
            disp('                             Nudge: [ and ] ')
            disp(' Zoom in and out.............Coarse zoom: / and "  ')
            disp('                             Fine zoom: up and down cursor arrows')
            disp(' Resets pan and zoom.........ESC')
            disp(' Select entire signal........Ctrl-A')
            disp(' Display grid................Shift-G  temporarily display grid on both panels')
            disp(' Adjust smooth width.........A,Z (A=>more, Z=>less) ')
            disp(' Cycle smooth types..........S (No, Rectanglular, Triangle, Gaussian, Savitzky-Golay)')
            disp(' Toggle smooth ends..........X (0=ends zeroed  1=ends smoothed (slower)')
            disp(' Cycle derivative orders.....D/Shift-D Increase/Decrease derivative order')
            disp(' Toggle peak sharpening......E (0=OFF 1=ON)')
            disp(' Sharpening for Gaussian.....Y  Set sharpen settings for Gaussian')
            disp(' Sharpening for Lorentzian...U  Set sharpen settings for Lorentzian')
            disp(' Adjust sharp1...............F,V  F=>sharper, V=>less sharpening')
            disp(' Adjust sharp2   ............G,B  G=>sharper, B=>less sharpening')
            disp(' Slew rate limit (0=OFF).....~  Largest allowed y change between points')
            disp(' Spike filter width (O=OFF)..m  spike filter eliminates sharp spikes')
            disp(' Toggle peak parabola........P  fits parabola to center, labels vertex')
            disp(' Fit polynomial to segment...Shift-o  Asks for polynomial order')
            disp(' Fits peak in upper window...Shift-F (Asks for shape, number of peaks, etc)')
            disp(' Spectrum mode on/off........Shift-S (Shift-A and Shift-X to change axes)')            
            disp(' Peak labels on spectrum.....Shift-Z in spectrum mode ')
            disp(' Click graph to print x,y....Shift-C  Click graph to print coordinates')
            disp(' Display Waterfall spectrum..Shift-W  Allows choice of mesh, surf, contour, or pcolor')
            disp(' Transfer power spectrum.....Shift-T  Replaces signal with power spectrum')
            disp(' Lock in current processing..Shift-L  Replace signal with processed version')
            disp(' ConVolution/DeconVolution...Shift-V  Convolution/Deconvolution menu')
            disp(' Print peak report...........R  prints position, height, width, area')
            disp(' Toggle log y mode...........H  semilog plot in lower window')
            disp(' Cycles baseline mode........T  none, linear, quadratic, or flat baseline mode')
            disp(' Restores original signal....Tab or Ctrl-Z key resets to original signal and modes')
            disp(' Toggle overlay mode.........L  Overlays original signal as dotted line')
            disp(' Baseline subtraction........Backspace, then click baseline at multiple points')
            disp(' Restore background..........\  to cancel previous background subtraction')
            disp(' Invert signal...............Shift-N  Invert (negate) the signal (flip + and -)')
            disp(' Remove offset...............0  (zero) set minimun signal to zero ')
            disp(' Sets region to zero.........;  sets selected region to zero')
            disp(' Absolute value..............+  Computes absolute value of entire signal')
            disp(' Condense signal.............C  Condense oversampled signal by factor of N')
            disp(' Interpolate signal..........i  Interpolate (resample) to N points')
            disp(' Print report................Q  prints signal info and current settings')
            disp(' Print keyboard commands.....K  prints this list')
            disp(' Print isignal arguments.....W  prints isignal function with all current arguments')
            disp(' Save output to disk.........O  Save .mat file with processed signal matrix')
            disp(' Play signal as sound........Spacebar or Shift-P  Play selected signal through speaker')
            disp(' Play signal as sound........Shift-R Change sampling rate for playing sound')
            disp(' Switch to ipf.m.............Shift-Ctrl-F  transfer current signal to Interactive Curve Fitter')
            disp(' Switch to iPeak.............Shift-Ctrl-P  transfer current signal to Interactive Peak Detector')
        case 119
            % When 'W' is pressed, prints 'isignal(current arguments)'
            firstpoint=xo+dx/2;
            if firstpoint>length(X),firstpoint=length(X);end
            lastpoint=xo-dx/2;
            if lastpoint<1,lastpoint=1;end
            disp(['isignal(DataMatrix,'  num2str(X(round(xo))) ',' num2str(X(round(firstpoint))-X(round(lastpoint)))  ',' num2str(SmoothMode)  ',' num2str(SmoothWidth) ',' num2str(ends) ',' num2str(DerivativeMode)  ',' num2str(Sharpen)  ',' num2str(Sharp1)  ',' num2str(Sharp2)  ',' num2str(SlewRate) ',' num2str(MedianWidth)  ',' num2str(SpectrumMode) ');' ] )
            disp(['peakfit(DataMatrix,'  num2str(X(round(xo))) ',' num2str(X(round(firstpoint))-X(round(lastpoint))) ')' ] )        
        case 111
            % When 'o' key is pressed, processed signal X,Y matrix is saved as in
            % mat file as the variable 'Output"
            Output=[X Y];
            uisave('Output');
            if SpectrumMode==1;
                [f realsy PowerSpectrum]=PlotFrequencySpectrum(X,Y,xo,dx,plotmode,xmode,LabelPeaks);
               uisave('PowerSpectrum');
            end 
        case 70
            % When 'Shift-F' is pressed, applies peakfit function only to
            %  peaks in the upper window.
            Startx=round(xo-(dx/2));
            Endx=abs(round(xo+(dx/2)-1));
            if Endx>length(Y),Endx=length(Y);end
            if Startx<1,Startx=1;end
            PlotRange=Startx:Endx;
            if (Endx-Startx)<2, PlotRange=xo:xo+2;end
            xx=X(PlotRange);
            yy=Y(PlotRange);

%             if Shape==11||12,
%                 fixedstart=[];
%                 for pk=1:NumPeaks,
%                     fixedstart(pk)=start(2*pk-1);
%                 end
%                 peakfit([xx;yy],0,0,NumPeaks,shapesvector,extra,1,start,AUTOZERO,FIXEDPARAMETERS,1,BIPOLAR,MINWIDTH,delta);
%             end
%             [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
            disp('Gaussians: y=exp(-((x-pos)./(0.6005615.*width)) .^2)')
            disp('  Gaussians with independent positions and widths...................1 (default)')
            disp('  Exponentially-broadened Gaussian (equal time constants)...........5 ')
            disp('  Exponentially-broadened equal-width Gaussian......................8 ')
            disp('  Fixed-width exponentionally-broadened Gaussian...................36 ')
            disp('  Exponentially-broadened Gaussian (independent time constants)....31 ')
            disp('  Gaussians with the same widths....................................6 ')
            disp('  Gaussians with preset fixed widths...............................11 ')
            disp('  Fixed-position Gaussians.........................................16 ')
            disp('  Asymmetrical Gaussians with unequal half-widths on both sides....14 ')
            disp('Lorentzians: y=ones(size(x))./(1+((x-pos)./(0.5.*width)).^2)')
            disp('  Lorentzians with independent positions and widths.................2 ')
            disp('  Exponentially-broadened Lorentzian...............................18 ')
            disp('  Equal-width Lorentzians...........................................7')
            disp('  Fixed-width Lorentzian...........................................12')
            disp('  Fixed-position Lorentzian........................................17')
            disp('Gaussian/Lorentzian blend (equal blends)...........................13')
            disp('  Fixed-width Gaussian/Lorentzian blend............................35')
            disp('  Gaussian/Lorentzian blend with independent blends)...............33')
            disp('Voigt profile with equal alphas)...................................20')
            disp('  Fixed-width Voigt profile with equal alphas......................34')
            disp('  Voigt profile with independent alphas............................30')
            disp('Logistic: n=exp(-((x-pos)/(.477.*wid)).^2); y=(2.*n)./(1+n).........3 ')
            disp('Pearson: y=ones(size(x))./(1+((x-pos)./((0.5.^(2/m)).*wid)).^2).^m..4')
            disp('  Fixed-width Pearson..............................................37')
            disp('  Pearson with independent shape factors, m........................32')
            disp('Breit-Wigner-Fano..................................................15')
            disp('Exponential pulse: y=(x-tau2)./tau1.*exp(1-(x-tau2)./tau1)..........9')
            disp('Alpha function: y=(x-spoint)./pos.*exp(1-(x-spoint)./pos);.........19')
            disp('Up Sigmoid (logistic function): y=.5+.5*erf((x-tau1)/sqrt(2*tau2)).10')
            disp('Down Sigmoid y=.5-.5*erf((x-tau1)/sqrt(2*tau2))....................23')
            disp('Triangular.........................................................21')
            disp(' ')
            disp(['Select the peak shape of the model from the table above (type 1-37 and press Enter key):'])
            disp(['Current shape number is ' num2str(Shape) '. Press Enter to keep.' ]) 
            Shapeinput=input('Peak shape number (1-37): ');
            if isempty(Shapeinput),
            else
                Shape=Shapeinput;
            end
            if Shape>37, Shape=37;end
            if Shape<1, Shape=1;end
            % disp(' ')
            switch Shape
                case 1
                    ShapeString='Gaussian';
                case 2
                    ShapeString='Lorentzian';
                case 3
                    ShapeString='logistic';
                case 4
                    ShapeString='Pearson';
                    disp(['Current shape number is ' num2str(extra) '. Press Enter to keep.' ])
                    inputextra=input('Shape factor: ');
                    if isempty(inputextra),
                    else
                        extra=inputextra;
                    end
                case 5
                    ShapeString='ExpGaussian';
                    disp(['Current shape number is ' num2str(extra) '. Press Enter to keep.' ]) 
                    inputextra=input('Exponentional factor: ');
                    if isempty(inputextra),
                    else
                        extra=inputextra;
                    end
                case 6
                    ShapeString='Equal-width Gaussian';
                case 7
                    ShapeString='Equal-width Lorentzian';
                case 8
                    ShapeString='Equal-width ExpGauss.';
                    disp(['Current shape number is ' num2str(extra) '. Press Enter to keep.' ]) 
                    inputextra=input('Exponentional factor: ');
                    if isempty(inputextra),
                    else
                        extra=inputextra;
                    end
                case 9
                    ShapeString='Exponental pulse';
                case 10
                    ShapeString='Sigmoid';
                case 11
                    ShapeString='Fixed-width Gaussian';
                    disp(['Current peak width vector is [' num2str(fixedparameters) ']. Press Enter to keep.' ])
                    inputwidth=input('Peak width vector, in brackets: ');
                    if isempty(inputwidth),
                    else
                        fixedparameters=inputwidth;
                    end
                case 12
                    ShapeString='Fixed-width Lorentzian';
                    disp(['Current peak width vector is [' num2str(fixedparameters) ']. Press Enter to keep.' ])
                    inputwidth=input('Peak width: ');
                    if isempty(inputwidth),
                    else
                        fixedparameters=inputwidth;
                    end
                case 13
                    ShapeString='Gauss/Lorentz blend';
                    disp(['Current percent Gaussian is ' num2str(extra) '. Press Enter to keep.' ])
                    inputextra=input('Percent Gaussian: ');
                    if isempty(inputextra),
                    else
                        extra=inputextra;
                    end
                case 14
                    ShapeString='bifurcated Gaussian';
                    disp(['Current asymmetry factor is ' num2str(extra) '. Press Enter to keep.' ])
                    inputextra=input('Asymmetry factor: ');
                    if isempty(inputextra),
                    else
                        extra=inputextra;
                    end
                case 15
                    ShapeString='Breit-Wigner-Fano';
                    disp(['Current asymmetry factor is ' num2str(extra) '. Press Enter to keep.' ])
                    inputextra=input('Asymmetry factor: ');
                    if isempty(inputextra),
                    else
                        extra=inputextra;
                    end
                case 16
                    ShapeString='Fixed-position Gaussians';
                    inputpositions=input('Peak positions as a vector, e.g. [200 400 600]: ');
                    if isempty(inputpositions),
                    else
                        fixedparameters=inputpositions;
                    end
                case 17
                    ShapeString='Fixed-position Lorentzians';
                    inputpositions=input('Peak positions as a vector, e.g. [200 400 600]: ');
                    if isempty(inputpositions),
                    else
                        fixedparameters=inputpositions;
                    end
                case 18
                    ShapeString='ExpLorentzian';
                    disp(['Current shape number is ' num2str(extra) '. Press Enter to keep.' ])
                    inputextra=input('Exponentional factor: ');
                    if isempty(inputextra),
                    else
                        extra=inputextra;
                    end
                case 19
                    ShapeString='Alpha function';
                case 20
                    ShapeString='Voigt profile';
                    disp(['Current alpha value is ' num2str(extra) '. Press Enter to keep.' ])
                    inputextra=input('Alpha: ');
                    if isempty(inputextra),
                    else
                        extra=inputextra;
                    end
                case 21
                    ShapeString='triangular';
                case 22
                    ShapeString=num2str(shapesvector);
                case 24
                    ShapeString='Negative Binomial Distribution';
                case 25
                    ShapeString='Lognormal Distribution';
                case 26
                    ShapeString='slope';
                case 27
                    ShapeString='Gaussian First derivative';
                case 28
                    ShapeString='Polynomial';
                case 29
                    ShapeString='Segmented linear';
                case 30
                    ShapeString='Voigt (variable alphas)';
                case 31
                    ShapeString='ExpGaussian (var. time constant)';
                case 32
                    ShapeString='Pearson (var. shape constant)';
                case 33
                    ShapeString='Variable Gaussian/Lorentzian';
                case 34
                    ShapeString='Fixed-width Voigt';
                    disp(['Current peak width vector is ' num2str(fixedparameters) '. Press Enter to keep.' ])
                    inputwidth=input('Peak width: ');
                    if isempty(inputwidth),
                    else
                        fixedparameters=inputwidth;
                    end
                case 35
                    ShapeString='Fixed-width G/L blend';
                    disp(['Current peak width vector is ' num2str(fixedparameters) '. Press Enter to keep.' ])
                    inputwidth=input('Peak width: ');
                    if isempty(inputwidth),
                    else
                        fixedparameters=inputwidth;
                    end
                case 36
                    ShapeString='Fixed-width ExpGaussian';
                    disp(['Current peak width vector is ' num2str(fixedparameters) '. Press Enter to keep.' ])
                    inputwidth=input('Peak width: ');
                    if isempty(inputwidth),
                    else
                        fixedparameters=inputwidth;
                    end
                case 37
                    ShapeString='Fixed-width Pearson';
                    disp(['Current peak width vector is ' num2str(fixedparameters) '. Press Enter to keep.' ])
                    inputwidth=input('Peak width: ');
                    if isempty(inputwidth),
                    else
                        fixedparameters=inputwidth;
                    end
                otherwise
                    ShapeString='';
            end % switch Shape''
                disp(' ')
                disp(['Current number of peaks is ' num2str(NumPeaksUW) '. Press Enter to keep.' ])
                inputNumPeaks=input('Number of peaks: ');
                if isempty(inputNumPeaks),
                else
                    NumPeaksUW=inputNumPeaks;
                end

            disp(' ')
            disp(['Current number of trials is ' num2str(NumTrials) '. Press Enter to keep.' ])
            inputNumTrials=input('Number of Trial fits: ');
            if isempty(inputNumTrials),
            else
                if isnumeric(inputNumTrials),
                    NumTrials=inputNumTrials;
                else
                    NumTrials=1;
                end
            end
            X1=min(xx);
            X2=max(xx);
            lyy=min(yy);
            uyy=max(yy)+(max(yy)-min(yy))/10;
            n=X2-X1;
            width=n/(5);
            % When 'c' key is pressed, user clicks graph to enter start positons,
            % then fit is computed and graph rte-drawn.
            % Acquire first-guess peak positions from user mouse pointer
            disp('Click on the estimated positions of each proposed component peak.')
            % figure(1);
            subplot(2,1,1);xlabel('Click on the estimated positions of each proposed component peak.')
            [clickX,clickY] = ginput(NumPeaksUW);
            % Create a "start" vector using these peak positions, with peak
            % widths
            n=max(xx)-min(xx);
            width=n/(5*NumPeaksUW);
            start=[];
            for k=1:NumPeaksUW,
                start=[start clickX(k) width];
            end
            disp(' ')
            disp(['Least-squares fit of selected peaks to ' ShapeString ' peak model using the peakfit function:' ])
            figure(2)
            % peakfit called
            [FitResults,MeanFitError]=peakfit([xx,yy],0,0,NumPeaksUW,Shape,extra,NumTrials,start,autozero,fixedparameters);
            disp(['Fitting Error = ' num2str(MeanFitError(1)) '%     R2 = ' num2str(MeanFitError(2)) ] )
            disp('          Peak#     Position     Height      Width         Area  ')
            % for peak=1:NumPeaksUW,FitResults(peak,1)=PUW(peak,1);end
            disp(FitResults(:,1:5))
            disp('Peakfit plot shown in Figure 2')
            figure(1)
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        case 83 % If Shift-S is pressed, plots frequency spectrum in the lower window
            if SpectrumMode==1, 
                SpectrumMode=0;
                [xx,yy]=RedrawSignal(X,Y,xo,dx);
            else
                SpectrumMode=1;
                % Plot the power spectrum  in the lower half of
                % the window.
                 [f,realsy,PowerSpectrum]=PlotFrequencySpectrum(X,Y,xo,dx,plotmode,xmode,LabelPeaks);
                subplot(2,1,1)
                title('iSignal 5   Frequency Spectrum Mode (Press Shift-S again to cancel')
            end  
        case 84 % If Shift-T is pressed, transfers spectrum to signal in upper window
            if SpectrumMode==1,
                [f,realsy,PowerSpectrum]=PlotFrequencySpectrum(X,Y,xo,dx,plotmode,xmode,LabelPeaks);
                X=f';
                Y=realsy;
                SavedSignal=Y;
            end
            SpectrumMode=0;
            xo=length(Y)/2; % Initial Pan setting
            dx=length(Y); % Initial Zoom setting
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
            
            disp('Signal replaced with its power spectrum. Press Tab to cancel.')
        case 90 % If Shift-Z is pressed, prints frequency spectrum peaks in the lower window
             if LabelPeaks==1, 
                 LabelPeaks=0;
                 [xx,yy]=RedrawSignal(X,Y,xo,dx);
             else
                LabelPeaks=1;
                % Plot the power spectrum  in the lower half of
                % the window.
                [f,realsy,PowerSpectrum]=PlotFrequencySpectrum(X,Y,xo,dx,plotmode,xmode,LabelPeaks);
             end  
        case 87 % Shift-W computes 3D matrix 
            n=input('Number of segments: ');
            disp('Type of 3D display:');
            disp('1 mesh');
            disp('2 contour');
            disp('3 pcolor');
            disp('4 surf');
            disp('5 Waterfall');
            plotype=input('Type 1 - 5: ');
            lxx=round(length(X)./n);
            for k=1:n-1,
                [f,realsy]=CompFrequencySpectrum(X,Y,k.*lxx,lxx,0,xmode);
                M(k,:)=log10(realsy);
            end
            figure(2);
            switch plotype
                case 1
                    mesh(M');
                case 2
                    contour(M');
                case 3
                    pcolor(M');
                case 4
                    surf(M');
                case 5
                    waterfall(M');
                otherwise
                   waterfall(M');
            end
        case 67  % Shift-C Click on plot points
                [clickX,clickY] = ginput(1);
                f=CompFrequencySpectrum(X,Y,xo,dx,0,xmode);
                disp('          x        y')
                disp([clickX,clickY])
        case 65 % If Shift-A is pressed, changes plot mode for Spectrum
            plotmode=plotmode+1;
            if plotmode==5;plotmode=1;end
             [f realsy PowerSpectrum]=PlotFrequencySpectrum(X,Y,xo,dx,plotmode,xmode,LabelPeaks);
            subplot(2,1,1)
            title('iSignal 5   Frequency Spectrum Mode (Press Shift-S again to cancel')
         case 88 % If Shift-X is pressed,changes xmode for Spectrum
            xmode=xmode+1;
            if xmode==2;xmode=0;end
             [f realsy PowerSpectrum]=PlotFrequencySpectrum(X,Y,xo,dx,plotmode,xmode,LabelPeaks);
            subplot(2,1,1)
           title('iSignal 5   Frequency Spectrum Mode (Press Shift-S again to cancel')
        case {80,32}  % Shift-P or Space bar; Play signal as sound through computer sound system
            sound(yy./max(yy),samplerate)
            figure(1)
        case 82   % Shift-R  Enter signal sampling rate
            samplerate=input('Sampling rate in Hz (Press Enter for 8000 Hz): ');
            if isempty(samplerate),samplerate=8000;end
            sound(yy./max(yy),samplerate)
            figure(1)
        case 79 % Shift-o
            disp(' ')
            polyorder=input('Polynomial order (1=linear, 2=quadratic, etc): ');
            if isempty(polyorder),polyorder=0;end
            if polyorder>0,
                clf
                figure(1);subplot(2,1,1); % Select upper window
                [coef, RSquared]=plotit(xx,yy,polyorder);
                residual=yy-polyval(coef,xx);
                subplot(2,1,2);plot(xx,residual,'r.')
                xlabel('Residuals')
                disp(' ')
                disp([ 'Selected x range: ' num2str(length(xx)) ' points from x = ' num2str(min(xx)) ' to ' num2str(max(xx)) ])
                disp([ 'Polynomial order: ' num2str(polyorder)])
                disp([ 'Polynomial coefficients (slope, intercept): ' num2str(coef)])
                disp([ 'Coefficient of determination (R2): ' num2str(RSquared)])
            end
        case 76 % Shift-L  Lock in current processing
            SavedSignal=Y;
            SmoothMode=0;
            SmoothWidth=1;
            SmoothType='No';
            DerivativeMode=0;
            Sharpen=0;
            SlewRate=0;
            MedianWidth=0;
            Y=ProcessSignal(X,Y,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
            subplot(2,1,1);
            title('Signal replaced with processed version; all settings reset.')
        case 86  % Shift-V Fourier convolution/deconvolution menu
            disp('Fourier convolution/deconvolution menu ')
            disp('  1. Convolution')
            disp('  2. Deconvolution')
            vmode=input('Select mode 1 or 2: ');
            disp(' ')
            disp('Shape of convolution/deconvolution function:')
            disp('  1. Gaussian')
            disp('  2. Exponential')
            smode=input('Select shape 1 or 2: ');
            disp(' ')
            switch vmode
                case 1 % Convolution
                    switch smode
                        case 1 % '1. Gaussian
                            vwidth=input('Enter the Gaussian width: ');
                            gc=gaussian(X,max(X)/2,vwidth); % convolution function
                            Y=conv(Y,gc,'same')./addup(gc);
                        case 2 % '2. Exponential'
                            vwidth=input('Enter the exponential time constant: ');
                            ec=exp(-X./vwidth);
                            Y=ifft(fft(Y).*fft(ec))./addup(ec);
                    end
                case 2 % Deconvolution
                    switch smode
                        case 1 % '1. Gaussian
                            vwidth=input('Enter the Gaussian width: ');
                            Y=deconvgauss(X,Y,vwidth);
                        case 2 % '2. Exponential'
                            vwidth=input('Enter the exponential time constant: ');
                            Y=deconvexp(X,Y,vwidth);
                    end
            end
            figure(1)
            SavedSignal=Y;
            Y=ProcessSignal(X,SavedSignal,DerivativeMode,SmoothWidth,SmoothMode,ends,Sharpen,Sharp1,Sharp2,SlewRate,MedianWidth);
            [xx,yy]=RedrawSignal(X,Y,xo,dx);
        otherwise
            UnassignedKey=double(key)
            disp('Press k to print out list of keyboard commands')
    end % switch double(key),
end % if  isscalar(key),
% ----------------------------------------------------------------------
function [xo,dx]=panandzoom(X,xcenter,xrange)
 xo=val2ind(X,xcenter);
 hirange=val2ind(X,xcenter+xrange);
 lorange=val2ind(X,xcenter-xrange);
 dx=(hirange-lorange)./2;
 if xcenter<min(X),
      disp(['Lowest X value is ' num2str(min(X)) ]),
      xcenter=min(X)+xrange;
 end
 if xcenter>max(X),
       disp(['Highest X value is ' num2str(max(X)) ]),
       xcenter=max(X)-xrange;
  end
% ----------------------------------------------------------------------    
function [xx,yy]=RedrawSignal(x,y,xo,dx)
% Plots the entire signal (x,y) in the lower half of the plot window and an
% isolated segment (xx,yy) in the upper half, controlled by Pan and Zoom
% keys.
global X SmoothType SmoothWidth Sharpen Sharp1 Sharp2 ends SpectrumMode LabelPeaks
global DerivativeMode PeakLabels Overlay Report autozero MedianWidth xmode
global SavedSignal GaussEstimate LorentzEstimate logymode SlewRate plotmode 
Startx=round(xo-(dx/2));
Endx=abs(round(xo+(dx/2)-1));
if Endx>length(y),Endx=length(y);end
if Startx<1,Startx=1;end
PlotRange=Startx:Endx;
if (Endx-Startx)<2, PlotRange=xo:xo+2;end
xx=x(PlotRange);
yy=y(PlotRange); 
datasize=size(yy);if datasize(1)>datasize(2),yy=yy';end
datasize=size(xx);if datasize(1)>datasize(2),xx=xx';end
% Remove local baseline from data segment
bkgsize=round(length(xx)/10);
if bkgsize<2,bkgsize=2;end
lxx=length(xx);
bkgcoef=0;
if autozero==1, % linear autozero operation  
    XX1=xx(1:round(lxx/bkgsize));
    XX2=xx((lxx-round(lxx/bkgsize)):lxx);
    Y1=yy(1:(round(length(xx)/bkgsize)));
    Y2=yy((lxx-round(lxx/bkgsize)):lxx);
    bkgcoef=polyfit([XX1,XX2],[Y1,Y2],1);  % Fit straight line to sub-group of points
    bkg=polyval(bkgcoef,xx);
    yy=yy-bkg;
end % if
if autozero==2, % Quadratic autozero operation  
    XX1=xx(1:round(lxx/bkgsize));
    XX2=xx((lxx-round(lxx/bkgsize)):lxx);
    Y1=yy(1:round(length(xx)/bkgsize));
    Y2=yy((lxx-round(lxx/bkgsize)):lxx);
    bkgcoef=polyfit([XX1,XX2],[Y1,Y2],2);  % Fit parabola to sub-group of points
    bkg=polyval(bkgcoef,xx);
    yy=yy-bkg;
end % if autozero
clf
% Plots isolated segment (xx,yy) in the upper half
% figure(1);
subplot(2,1,1);
if Overlay,
    hold on
    plot(xx,SavedSignal(PlotRange),'b:'); 
end % Overlay
plot(xx,yy,'b')
    switch autozero,
        case 0
            title('iSignal 5. No baseline correction.  Press K for keyboard commands')
        case 1
            title('iSignal 5. Linear baseline subtraction.  Press K for keyboard commands')
        case 2
            title('iSignal 5. Quadratic subtraction baseline.  Press K for keyboard commands')
        case 3
            title('iSignal 5. Flat baseline correction.  Press K for keyboard commands')
    end
lyy=min(yy);
uyy=max(yy)+(max(yy)-min(yy))/10;
if lyy<uyy;
   axis([x(Startx) x(Endx) lyy uyy ]);
end
center=x(round(xo));
hold on;plot([center center],[lyy uyy],'g-')
yyrange=max(yy)-min(yy);
xlabel(sprintf('%s    y: %0.3g at %0.5g      P/P: %0.3g      Area: %0.3g    Std. Dev.: %0.2g ',' ',y(round(xo)), center, yyrange, trapz(xx,yy), std(yy) ))

% Bottom half of the figure shows full signal
subplot(2,1,2);cla
if SpectrumMode,
    PlotFrequencySpectrum(x,y,xo,dx,plotmode,xmode,LabelPeaks);
else
    if Overlay,
        hold on
        if logymode,
            semilogy(x,SavedSignal,'b:');
        else
            plot(x,SavedSignal,'b:');
        end % if logymode
    end % Overlay
    if logymode,
        semilogy(x,y,'b')  % Graph the signal with linear y axis on lower half
    else
        plot(x,y,'b')  % Graph the signal with linear y axis on lower half
    end % if logymode
    axis([x(1) x(length(x)) min(min(y)) max(max(y))]); % Update plot
    title('Smooth: S, A/Z   Derivatives: D   Peak Meas: P   Spectrum: Shift-S   lin/log: H   Baseline mode: T   ')
    if Sharpen,
        xlabel([ num2str(SmoothWidth) ' point '  SmoothType ' smooth. Ends: ' num2str(ends) '   Der: ' num2str(DerivativeMode)  '   S1: ' num2str(round(10*Sharp1)/10)  '   S2: ' num2str(round(100*Sharp2)/100) '   Slew: ' num2str(SlewRate) '   Median: ' num2str(MedianWidth) ])
    else
        xlabel([ num2str(SmoothWidth) ' point '  SmoothType ' smooth. Ends: ' num2str(ends) '    Der: ' num2str(DerivativeMode) '    Slew: ' num2str(SlewRate) '    Median: ' num2str(MedianWidth) ])
    end
    if logymode,
        ylabel('Log y')
    else
        ylabel('Linear y')
    end
    hold on
    % Mark the zoom range on the full signal with two magenta dotted vertical lines
    checkzero=abs(y);
    checkzero(~checkzero)=NaN; % Find smallest non-zero value
    plot([min(xx) min(xx)],[min(y) max(y)],'m--')
    plot([max(xx) max(xx)],[min(y) max(y)],'m--')
    plot([center center],[min(y) max(y)],'g-')
    if PeakLabels,
        % Compute limited range around center of zoom region
        FitLow=round(xo-(dx/10));
        if FitLow<1,FitLow=1;end
        FitHigh=abs(round(xo+(dx/10)-1));
        if FitHigh>length(X),FitHigh=length(X);end
        FitRange=FitLow:FitHigh;
        xxx=x(FitRange);
        yyy=y(FitRange);
        datasize=size(yyy);
        if datasize(1)<datasize(2),yyy=yyy';end
        datasize=size(xxx);
        if datasize(1)<datasize(2),xxx=xxx';end
        if autozero==1||2,
            yyy=yyy-polyval(bkgcoef,xxx);
        end
        [coef,S]=polyfit(xxx,yyy,2);  % Fit parabola to sub-group of points
        c1=coef(3);c2=coef(2);c3=coef(1);
        % Compute the correlation coefficient and R-Squared
        cc=corrcoef(polyval(coef,xxx),yyy);
        if size(cc)>1;
            RSquared=cc(2).^2;
        else
            RSquared=0;
        end
        % disp([c1 c2 c3 RSquared])  % for testing
        PeakX=-c2/(2*c3); % x-value of vertex
        PeakY=(c1-(c2*c2/(4*c3))); % y-value of vertex
        MeasuredWidth=norm(2.35482/(sqrt(2)*sqrt(-1*c3))); % Estimated Gaussian half-width
        % Label the peaks on the upper graph with number, position, height, and
        % width
        residual=yyy-polyval(coef,xxx);
        SNR=abs(PeakY./std(residual));
        if c2>0,
            % Fit parabola to log10 of sub-group
            [coef,S,MU]=polyfit(xxx,rmnan(log(yyy)),2);
            d1=coef(3);d2=coef(2);d3=coef(1);
            % Compute peak position and height of fitted parabola
            PeakX=-((MU(2).*d2/(2*d3))-MU(1));
            PeakY=exp(d1-d3*(d2/(2*d3))^2);
            MeasuredWidth=norm(MU(2).*2.35482/(sqrt(2)*sqrt(-1*d3)));
            % cc=corrcoef(polyval(coef,xxx),log(abs(yyy)));
            % RSquared=cc(2).^2;
            residual=yyy-PeakY*gaussian(xxx,PeakX,MeasuredWidth);
            SNR=abs(PeakY./std(residual));
        end
        subplot(2,1,1);
        hold on
        topaxis=axis;
        hpos=min(xx);
        yrange=topaxis(4)-topaxis(3);
        pos1=.1*yrange;
        pos2=.2*yrange;
        pos3=.3*yrange;
        pos4=.4*yrange;
        pos5=.5*yrange;
        pos6=.6*yrange;
        pos7=.7*yrange;
        pos8=.8*yrange;
        pos9=.9*yrange;
        text(hpos,topaxis(4)-pos1,[' Position=' num2str(PeakX)])
        text(hpos,topaxis(4)-pos2,[' Height=' num2str(PeakY)])
        text(hpos,topaxis(4)-pos3,[' Gaussian Width=' num2str(MeasuredWidth)]) 
        text(hpos,topaxis(4)-pos5,[' Gaussian area=' num2str(1.0645*PeakY*MeasuredWidth) ])
        area=trapz(xx,yy); % Compute the area of displayed segment
        text(hpos,topaxis(4)-pos6,[' Displayed area=' num2str(area) ])
        text(hpos,topaxis(4)-pos7,[' Percent of total area=' num2str(100*area./trapz(x,y)) ])
        text(hpos,topaxis(4)-pos8,[' R2=' num2str(RSquared)])
        text(hpos,topaxis(4)-pos9,[' SNR=' num2str(round(10.*SNR)/10) ])
        try
            FWHM=halfwidth(xx,yy);
            text(hpos,topaxis(4)-pos4,[' FWHM=' num2str(FWHM)])
        catch  
        end
        if Report,
        try
            FWHM=halfwidth(xx,yy);
        catch  
        end
            % disp([PeakX PeakY MeasuredWidth 1.0645*PeakY*MeasuredWidth area SNR  FWHM]);
            disp(sprintf('%0.5g       %0.5g       %0.5g       %0.5g       %0.5g       %0.3g       %0.5g',PeakX, PeakY, MeasuredWidth, 1.0645*PeakY*MeasuredWidth, area, SNR, FWHM));
            Report=0;
        end % Report
        plotspace=linspace(min(xxx),max(xxx));
        if c2>0,
            plot(plotspace,PeakY.*exp(-((plotspace-PeakX)./(0.6005615.*MeasuredWidth)).^2),'r')
        else
            plot(plotspace,c3.*plotspace.^2+c2.*plotspace+c1,'r')
        end
        hold off
        if xo<2;xo=2;end
        xinterval=X(round(xo))-X(round(xo-1));
        if GaussEstimate,
            SmoothWidth=round(0.4*MeasuredWidth./xinterval);
            Sharp1=((MeasuredWidth/xinterval)^2)/25;
            Sharp2=((MeasuredWidth/xinterval)^4)/800;
            GaussEstimate=0;
        end % if GaussEstimate
        if LorentzEstimate,
            SmoothWidth=round(0.3*MeasuredWidth./xinterval);
            Sharp1=((MeasuredWidth/xinterval)^2)/8;
            Sharp2=((MeasuredWidth/xinterval)^4)/700;
            LorentzEstimate=0;
        end % if LorentzEstimate
    end  % PeakLabels
end % if SpectrumMode
hold off
% ----------------------------------------------------------------------
function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Tom O'Haver (toh@umd.edu) October 2006
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);
% ----------------------------------------------------------------------
function Processed=ProcessSignal(x,y,DerivativeMode,w,type,ends,Sharpen,factor1,factor2,SlewRate,MedianWidth)
global X
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
% following lines for testing
% disp(['X,Y size = ' num2str(size(x)) ' , '  num2str(size(Processed)) ] )
% sizeProcessed=size(Processed)
% sizeX=size(X)
Processed=reshape(Processed,size(X));
% datasize=size(Processed);if datasize(1)>datasize(2),Processed=Processed';end
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
    case 4
       SmoothY=sa(sa(sa(sa(Y,w,ends),w,ends),w,ends),w,ends);
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
  end
% ----------------------------------------------------------------------
function sy=condense(y,n)
% Condense y by a factor of n, where n is a non-zero positive integer.
% Produces a shorter, approximate version of vector y, with each group 
% of n adjacent points in y replaced by its average. Use for reducing the 
% length and processing time of over-sampled signals or for preliminary 
% and exploratory analysis of very large signals to locate the interesting 
% bits, which can then be selected out of the full-length signal for 
% more precise analysis. For x,y data sets, use this function on both 
% independent variable x AND dependent variable y so that the features 
% of y will appear at the same x values.
% Example: condense([1 2 3 4 5 6 7 8 9 10 11 12],3) yields [2 5 8 11]
% condense([.9 1.1 .9 1 .9 1.1 .9 1 .9 1.1 .9 1],3) = [0.9667 1 0.9333 1]
% condense([0 0 0 0 0 0 0 1 1 1 1 1 1 1],2) = [0 0 0 .5 1 1 1]
n=round(n);
m=floor(length(y)/n);
if n > 1
    sy=mean(reshape(y(1:n*m),n,m));
else
    sy=y;
end
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.60056120439323.*wid)) .^2);
% ----------------------------------------------------------------------
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
% ----------------------------------------------------------------------
function [f realsy PowerSpectrum]=PlotFrequencySpectrum(X,Y,xo,dx,plotmode,XMODE,LabelPeaks)
global SmoothType SmoothWidth Sharpen Sharp1 Sharp2 ends 
global DerivativeMode MedianWidth SlewRate 
% xodx=[xo dx]
Startx=round(xo-(dx/2));
Endx=abs(round(xo+(dx/2)-1));
if Endx>length(Y),Endx=length(Y);end
if Startx<1,Startx=1;end
PlotRange=Startx:Endx;
if (Endx-Startx)<2, PlotRange=xo:xo+2;end
xx=X(PlotRange);
yy=Y(PlotRange);
subplot(2,1,1)
title('iSignal 5 Frequency Spectrum Mode (Press Shift-S again to cancel)')
xlabel('Press Shift-A to cycle through spectrum log/linear plot modes ')
subplot(2,1,2)
fy=fft(yy);
sy=fy .* conj(fy); % Compute power spectrum
plotrange=2:length(fy)/2;
if XMODE,
    f=range(xx)./(plotrange-1);
else
    f=((plotrange-1)./range(xx));
end
realsy=real(sy(plotrange));
maxpower=max(realsy);
maxx=val2ind(realsy,maxpower);
maxf=f(maxx);
hold off
switch plotmode,
    case 1
        plot(f,realsy,'r.-')
        ylabel('Linear y')
    case 2
        semilogx(f,realsy,'r.-')
        ylabel('Linear y')
    case 3
        semilogy(f,realsy,'r.-')
        ylabel('Log y')
    case 4
        loglog(f,realsy,'r.-')
        ylabel('Log y')
    otherwise,
end
% spectrumaxis=axis;
% hpos=min(realsy);
% text(spectrumaxis(1),0.5.*spectrumaxis(4),['Peak ' num2str(maxf) ' at harmonic #' num2str(maxx) ])
 if Sharpen,
        xlabel([ num2str(SmoothWidth) ' point '  SmoothType ' smooth. Ends: ' num2str(ends) '   Der: ' num2str(DerivativeMode)  '   S1: ' num2str(round(10*Sharp1)/10)  '   S2: ' num2str(round(100*Sharp2)/100) '   Slew: ' num2str(SlewRate) '   Median: ' num2str(MedianWidth) ])
    else
        xlabel([ num2str(SmoothWidth) ' point '  SmoothType ' smooth. Ends: ' num2str(ends) '    Der: ' num2str(DerivativeMode) '    Slew: ' num2str(SlewRate) '    Median: ' num2str(MedianWidth) ])
 end
if XMODE,
    title('x=Time. Press Shift-X to change to frequency.')
else
    title('x=Frequency (e.g. 1/time). Press Shift-X to change to time.')
end
PowerSpectrum=([f' realsy]);
if LabelPeaks % If Shift-Z has been pressed to set LabelPeaks=1
    % Change the values in the next 4 lines to adjust peak detection
    AmpT=0.03.*max(PowerSpectrum(:,2));
    SlopeT=0.000001;
    SmoothW=2;
    FitW=3;
    xaxis=((plotrange-1)./range(xx));
    P=findpeaksG(xaxis,realsy,SlopeT,AmpT,SmoothW,FitW);
    if XMODE,
       subplot(2,1,2)
       text(1./P(:,2),P(:,3),num2str(1./P(:,2)))
       subplot(2,1,1)
       title('iSignal 5 Frequency Spectrum Mode (Press Shift-S again to cancel)')
    else
       subplot(2,1,2)
       text(P(:,2),P(:,3),num2str(P(:,2)))
       subplot(2,1,1)
       title('iSignal 5 Frequency Spectrum Mode (Press Shift-S again to cancel)')
    end  
end
% ----------------------------------------------------------------------
function [f,realsy]=CompFrequencySpectrum(X,Y,xo,dx,plotmode,XMODE)
Startx=round(xo-(dx/2));
Endx=abs(round(xo+(dx/2)-1));
if Endx>length(Y),Endx=length(Y);end
if Startx<1,Startx=1;end
PlotRange=Startx:Endx;
if (Endx-Startx)<2, PlotRange=xo:xo+2;end
xx=X(PlotRange);
yy=Y(PlotRange);
fy=fft(yy);
sy=fy .* conj(fy); % Compute power spectrum
plotrange=2:length(fy)/2;
if XMODE,
    f=range(xx)./(plotrange-1);
else
    f=((plotrange-1)./range(xx));
end
realsy=real(sy(plotrange));
% ----------------------------------------------------------------------
function [FitResults,GOF,baseline,coeff,residual,xi,yi,BootResults]=peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero,fixedparameters,plots,bipolar,minwidth,DELTA,clipheight)

% A command-line peak fitting program for time-series signals.
% Version 7: March, 2015, adds peak shapes with three unconstrained
% iterated variables: 30=voigt (variable alpha), 31=ExpGaussian (variable
% time constant), 32=Pearson (variable shape factor), 34=Gaussian/
% Lorentzian blend (variable percent). See Examples 25-28 below.
%
% For more details, see
% http://terpconnect.umd.edu/~toh/spectrum/CurveFittingC.html and
% http://terpconnect.umd.edu/~toh/spectrum/InteractivePeakFitter.htm
% 
global AA xxx PEAKHEIGHTS FIXEDPARAMETERS AUTOZERO delta BIPOLAR CLIPHEIGHT
format short g
format compact
warning off all
NumArgOut=nargout;
datasize=size(signal);
if datasize(1)<datasize(2),signal=signal';end
datasize=size(signal);
if datasize(2)==1, %  Must be isignal(Y-vector)
    X=1:length(signal); % Create an independent variable vector
    Y=signal;
else
    % Must be isignal(DataMatrix)
    X=signal(:,1); % Split matrix argument 
    Y=signal(:,2);
end
X=reshape(X,1,length(X)); % Adjust X and Y vector shape to 1 x n (rather than n x 1)
Y=reshape(Y,1,length(Y));
% If necessary, flip the data vectors so that X increases
if X(1)>X(length(X)),
    disp('X-axis flipped.')
    X=fliplr(X);
    Y=fliplr(Y);
end

% Isolate desired segment from data set for curve fitting
if nargin==1 || nargin==2,center=(max(X)-min(X))/2;window=max(X)-min(X);end
% Y=Y-min(Y);
xoffset=center-window;  % <<<<<<<<<<<<
n1=val2ind(X,center-window/2);
n2=val2ind(X,center+window/2);
if window==0,n1=1;n2=length(X);end
xx=X(n1:n2)-xoffset;
% minxx=xx(1)

yy=Y(n1:n2);
ShapeString='Gaussian';
coeff=0;
CLIPHEIGHT=max(Y);
LOGPLOT=0;
% Define values of any missing arguments
switch nargin
    case 1  % Only data specified
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        xx=X;yy=Y;
        start=0;
        AUTOZERO=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=xx(2)-xx(1);
        delta=1;
        CLIPHEIGHT=max(Y);
    case 2
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        xx=signal;yy=center;
        start=0;
        AUTOZERO=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=xx(2)-xx(1);
        delta=1;
        CLIPHEIGHT=max(Y);
    case 3
        NumPeaks=1;
        peakshape=1;
        extra=0;
        NumTrials=1;
        start=0;
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=xx(2)-xx(1);
        delta=1;
        CLIPHEIGHT=max(Y);
    case 4 % NumPeaks specified in arguments
        peakshape=1;
        extra=0;
        NumTrials=1;
        start=0;
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=xx(2)-xx(1);
        delta=1;
        CLIPHEIGHT=max(Y);
    case 5 % NumPeaks, peakshape specified in arguments
        extra=zeros(1,NumPeaks);
        NumTrials=1;
        start=0;
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
        CLIPHEIGHT=max(Y);
    case 6  % extra, Numpeaks, peakshape specified in arguments
        NumTrials=1;
        start=0;
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
    case 7 % extra, Numpeaks, peakshape specified in arguments 
        start=0;
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
        CLIPHEIGHT=max(Y);
    case 8 % start, extra, Numpeaks, peakshape included in input arguments
        AUTOZERO=0;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
        CLIPHEIGHT=max(Y);
    case 9
        AUTOZERO=autozero;
        FIXEDPARAMETERS=0;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
    case 10
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        plots=1;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
    case 11
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        BIPOLAR=0;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
        CLIPHEIGHT=max(Y);
    case 12
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        BIPOLAR=bipolar;
        MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
        delta=1;
        CLIPHEIGHT=max(Y);
    case 13
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        BIPOLAR=bipolar;
        MINWIDTH=minwidth;
        delta=1;
    case 14
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        BIPOLAR=bipolar;
        MINWIDTH=minwidth;
        delta=DELTA;
        CLIPHEIGHT=max(Y);
    case 15
        AUTOZERO=autozero;
        FIXEDPARAMETERS=fixedparameters;
        BIPOLAR=bipolar;
        MINWIDTH=minwidth;
        delta=DELTA;
        CLIPHEIGHT=clipheight;
    otherwise
end % switch nargin

% Saturation Code, skips points greater than set maximum
if CLIPHEIGHT<max(Y),
    apnt=1;
    for pnt=1:length(xx),
        if yy(pnt)<CLIPHEIGHT,
            axx(apnt)=xx(pnt);
            ayy(apnt)=yy(pnt);
            apnt=apnt+1;
        end
    end
    xx=axx;yy=ayy;
end
% Default values for placeholder zeros1
if NumTrials==0;NumTrials=1;end
shapesvector=peakshape;
if isscalar(peakshape),
else
    % disp('peakshape is vector');
    shapesvector=peakshape;
    NumPeaks=length(peakshape);
    peakshape=22;
end
if peakshape==0;peakshape=1;end
if NumPeaks==0;NumPeaks=1;end

firststart=start; % <<<<<<<<<<<
if start==0;
    start=calcstart(xx,NumPeaks,xoffset);
else
     for Peak=1:NumPeaks,
         newstart(2*Peak-1)=start(2*Peak-1)-xoffset;
         newstart(2*Peak)=start(2*Peak);
     end
     start=newstart;
end
newstart=start; % <<<<<<<<<<<
if FIXEDPARAMETERS==0, FIXEDPARAMETERS=length(xx)/10;end
if peakshape==16;FIXEDPOSITIONS=fixedparameters;end
if peakshape==17;FIXEDPOSITIONS=fixedparameters;end
if AUTOZERO>3,AUTOZERO=3;disp('AUTOZERO must be between 0 and 3');end
if AUTOZERO<0,AUTOZERO=0;disp('AUTOZERO must be between 0 and 3');end
Heights=zeros(1,NumPeaks);
FitResults=zeros(NumPeaks,6);

% % Remove linear baseline from data segment if AUTOZERO==1
baseline=0;
bkgcoef=0;
bkgsize=round(length(xx)/10);
if bkgsize<2,bkgsize=2;end
lxx=length(xx);
if AUTOZERO==1, % linear autozero operation  
    XX1=xx(1:round(lxx/bkgsize));
    XX2=xx((lxx-round(lxx/bkgsize)):lxx);
    Y1=yy(1:(round(length(xx)/bkgsize)));
    Y2=yy((lxx-round(lxx/bkgsize)):lxx);
    bkgcoef=polyfit([XX1,XX2],[Y1,Y2],1);  % Fit straight line to sub-group of points
    bkg=polyval(bkgcoef,xx);
    yy=yy-bkg;
end % if
if AUTOZERO==2, % Quadratic autozero operation  
    XX1=xx(1:round(lxx/bkgsize));
    XX2=xx((lxx-round(lxx/bkgsize)):lxx);
    Y1=yy(1:round(length(xx)/bkgsize));
    Y2=yy((lxx-round(lxx/bkgsize)):lxx);
    bkgcoef=polyfit([XX1,XX2],[Y1,Y2],2);  % Fit parabola to sub-group of points
    bkg=polyval(bkgcoef,xx);
    yy=yy-bkg;
end % if autozero

PEAKHEIGHTS=zeros(1,NumPeaks);
n=length(xx);
% newstart=start;
% Assign ShapStrings
switch peakshape(1)
    case 1
        ShapeString='Gaussian';
    case 2
        ShapeString='Lorentzian';
    case 3
        ShapeString='Logistic';
    case 4
        ShapeString='Pearson';
    case 5
        ShapeString='ExpGaussian';
    case 6
        ShapeString='Equal width Gaussians';
    case 7
        ShapeString='Equal width Lorentzians';
    case 8
        ShapeString='Exp. equal width Gaussians';
    case 9
        ShapeString='Exponential Pulse';
    case 10
        ShapeString='Up Sigmoid (logistic function)';
    case 23
        ShapeString='Down Sigmoid (logistic function)';  
    case 11
        ShapeString='Fixed-width Gaussian';
    case 12
        ShapeString='Fixed-width Lorentzian';
    case 13
        ShapeString='Gaussian/Lorentzian blend';
    case 14
        ShapeString='BiGaussian';    
    case 15
        ShapeString='Breit-Wigner-Fano';   
    case 16
        ShapeString='Fixed-position Gaussians';
    case 17
        ShapeString='Fixed-position Lorentzians';
    case 18
        ShapeString='Exp. Lorentzian';
    case 19
        ShapeString='Alpha function';
    case 20
        ShapeString='Voigt (equal alphas)';
    case 21
        ShapeString='triangular';
    case 22
        ShapeString=num2str(shapesvector);
    case 24
        ShapeString='Negative Binomial Distribution';
    case 25
        ShapeString='Lognormal Distribution';
    case 26
        ShapeString='slope';
    case 27
        ShapeString='First derivative';
    case 28
        ShapeString='Polynomial';
    case 29
        ShapeString='Segmented linear';
    case 30
        ShapeString='Voigt (variable alphas)';
    case 31
        ShapeString='ExpGaussian (var. time constant)';
    case 32
        ShapeString='Pearson (var. shape constant)';
    case 33
        ShapeString='Variable Gaussian/Lorentzian';
    case 34
        ShapeString='Fixed-width Voigt';
    case 35
        ShapeString='Fixed-width G/L blend';
    case 36
        ShapeString='Fixed-width ExpGaussian';
    case 37
        ShapeString='Fixed-width Pearson';
    otherwise
end % switch peakshape
  
% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.00001,'TolFun',.00001,'Display','off','MaxFunEvals',1000 );
LowestError=1000; % or any big number greater than largest error expected
FitParameters=zeros(1,NumPeaks.*2); 
BestStart=zeros(1,NumPeaks.*2); 
height=zeros(1,NumPeaks); 
bestmodel=zeros(size(yy));

for k=1:NumTrials, 
    % StartMatrix(k,:)=newstart;
    % disp(['Trial number ' num2str(k) ] ) % optionally prints the current trial number as progress indicator
    switch peakshape(1)
        case 1
            TrialParameters=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 2
            TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 3
            TrialParameters=fminsearch(@(lambda)(fitlogistic(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 4
            TrialParameters=fminsearch(@(lambda)(fitpearson(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 5
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexpgaussian(lambda,zxx,zyy,-extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 6
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitewgaussian(lambda,xx,yy)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 7
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitewlorentzian(lambda,xx,yy)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 8
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 9
            TrialParameters=fminsearch(@(lambda)(fitexppulse(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 10
            TrialParameters=fminsearch(@(lambda)(fitupsigmoid(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 23
            TrialParameters=fminsearch(@(lambda)(fitdownsigmoid(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 11
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFWGaussian(lambda,xx,yy)),fixedstart,options);
        case 12
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFWLorentzian(lambda,xx,yy)),fixedstart,options);
        case 13
            TrialParameters=fminsearch(@(lambda)(fitGL(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 14
            TrialParameters=fminsearch(@(lambda)(fitBiGaussian(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 15
            TrialParameters=fminsearch(@(lambda)(fitBWF(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 16
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
                fixedstart(pc)=fixedstart(pc)+.1*(rand-.5).*fixedstart(pc);
            end
            TrialParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(Peak)<MINWIDTH,
                    TrialParameters(Peak)=MINWIDTH;
                end
            end
        case 17
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
                fixedstart(pc)=fixedstart(pc)+.1*(rand-.5).*fixedstart(pc);
            end
            TrialParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(Peak)<MINWIDTH,
                    TrialParameters(Peak)=MINWIDTH;
                end
            end
        case 18
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[ones(size(yy)).*yy(1) yy zeros(size(yy)).*yy(length(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexplorentzian(lambda,zxx,zyy,-extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 19
            TrialParameters=fminsearch(@(lambda)(fitalphafunction(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 20
            TrialParameters=fminsearch(@(lambda)(fitvoigt(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 21
            TrialParameters=fminsearch(@(lambda)(fittriangular(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 22
            TrialParameters=fminsearch(@(lambda)(fitmultiple(lambda,xx,yy,shapesvector,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH(Peak),
                    TrialParameters(2*Peak)=MINWIDTH(Peak);
                end
            end
        case 24
            TrialParameters=fminsearch(@(lambda)(fitnbinpdf(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 25
            TrialParameters=fminsearch(@(lambda)(fitlognpdf(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 26
            TrialParameters=fminsearch(@(lambda)(fitlinslope(lambda,xx,yy)),polyfit(xx,yy,1),options);
             coeff=TrialParameters;
        case 27
            TrialParameters=fminsearch(@(lambda)(fitd1gauss(lambda,xx,yy)),newstart,options);
        case 28
            coeff=fitpolynomial(xx,yy,extra);
            TrialParameters=coeff;
        case 29
            cnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cnewstart(pc)=newstart(2.*pc-1)+(delta*(rand-.5)/50);
            end
            TrialParameters=fminsearch(@(lambda)(fitsegmented(lambda,xx,yy)),cnewstart,options);
        case 30
            nn=max(xx)-min(xx);
            start=[];
            startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker)+ xoffset;
                start=[start markx nn/5 extra];
            end % for marker
             newstart=start;
            for parameter=1:3:3*NumPeaks,
                newstart(parameter)=newstart(parameter)*(1+randn/100);
                newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
            end
            TrialParameters=fminsearch(@(lambda)(fitvoigtv(lambda,xx,yy)),newstart);
         case 31
            nn=max(xx)-min(xx);
            start=[];
            startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker)+ xoffset;
                start=[start markx nn/5 extra];
            end % for marker
             newstart=start;
            for parameter=1:3:3*NumPeaks,
                newstart(parameter)=newstart(parameter)*(1+randn/100);
                newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
            end
            % newstart=newstart
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[ones(size(yy)).*yy(1) yy zeros(size(yy)).*yy(length(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexpgaussianv(lambda,zxx,zyy)),newstart);
        case 32
            nn=max(xx)-min(xx);
            start=[];
            startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker)+ xoffset;
                start=[start markx nn/5 extra];
            end % for marker
             newstart=start;
            for parameter=1:3:3*NumPeaks,
                newstart(parameter)=newstart(parameter)*(1+randn/100);
                newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
            end
            % newstart=newstart
            TrialParameters=fminsearch(@(lambda)(fitpearsonv(lambda,xx,yy)),newstart);
        case 33
             nn=max(xx)-min(xx);
            start=[];
            startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker)+ xoffset;
                start=[start markx nn/5 extra];
            end % for marker
             newstart=start;
            for parameter=1:3:3*NumPeaks,
                newstart(parameter)=newstart(parameter)*(1+randn/100);
                newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
            end
            % newstart=newstart
            TrialParameters=fminsearch(@(lambda)(fitGLv(lambda,xx,yy)),newstart);
        case 34
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end            
            TrialParameters=fminsearch(@(lambda)(fitFWVoigt(lambda,xx,yy,extra)),fixedstart,options);
        case 35
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end            
            TrialParameters=fminsearch(@(lambda)(fitFWGL(lambda,xx,yy,extra)),fixedstart,options);
        case 36
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end            
            TrialParameters=fminsearch(@(lambda)(fitFWExpGaussian(lambda,xx,yy,extra)),fixedstart,options);
        case 37
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end            
            TrialParameters=fminsearch(@(lambda)(fitFWPearson(lambda,xx,yy,extra)),fixedstart,options);
        otherwise
    end % switch peakshape

% Construct model from Trial parameters
A=zeros(NumPeaks,n);
for m=1:NumPeaks,
    switch peakshape(1)
        case 1
            A(m,:)=gaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 2
            A(m,:)=lorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 3
            A(m,:)=logistic(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 4
            A(m,:)=pearson(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 5
            A(m,:)=expgaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
        case 6
            A(m,:)=gaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
        case 7
            A(m,:)=lorentzian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
        case 8
            A(m,:)=expgaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1),-extra)';
        case 9
            A(m,:)=exppulse(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 10
            A(m,:)=upsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 11
            A(m,:)=gaussian(xx,TrialParameters(m),FIXEDPARAMETERS(m));
        case 12
            A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDPARAMETERS(m));
        case 13
            A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 14
            A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
        case 15
            A(m,:)=BWF(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);        
        case 16
            A(m,:)=gaussian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
        case 17
            A(m,:)=lorentzian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
        case 18
            A(m,:)=explorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
        case 19
            A(m,:)=alphafunction(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 20
            A(m,:)=voigt(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);        
        case 21
            A(m,:)=triangular(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 22
            A(m,:)=peakfunction(shapesvector(m),xx,TrialParameters(2*m-1),TrialParameters(2*m),extra(m));        
        case 23
            A(m,:)=downsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));        
        case 24
            A(m,:)=nbinpdf(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 25
            A(m,:)=lognormal(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 26
            A(m,:)=linslope(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 27
            A(m,:)=d1gauss(xx,TrialParameters(2*m-1),TrialParameters(2*m));
        case 28
            A(m,:)=polynomial(xx,coeff);
        case 29
            A(m,:)=segmented(xx,yy,PEAKHEIGHTS);
        case 30
            A(m,:)=voigt(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
        case 31
            A(m,:)=expgaussian(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),-TrialParameters(3*m));        
        case 32
            A(m,:)=pearson(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
        case 33
            A(m,:)=GL(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
        case 34
             width(m)=abs(FIXEDPARAMETERS(m));
%                 gD(m)=width(m);
%                 gL(m)=extra.*gD(m);
%                 width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2))
            A(m,:)=Voigt(xx,TrialParameters(m), width(m),extra);
        case 35
            A(m,:)=GL(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);    
        case 36
            A(m,:)=expgaussian(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);    
        case 37
            A(m,:)=pearson(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);    
    end % switch
    for parameter=1:2:2*NumPeaks,
        newstart(parameter)=newstart(parameter)*(1+delta*(rand-.5)/50);
        newstart(parameter+1)=newstart(parameter+1)*(1+delta*(rand-.5)/100);
    end
end % for NumPeaks
% newstart=newstart
% Multiplies each row by the corresponding amplitude and adds them up
if peakshape(1)==29, % Segmented linear
    model=segmented(xx,yy,PEAKHEIGHTS);
    TrialParameters=PEAKHEIGHTS;
    Heights=ones(size(PEAKHEIGHTS));
else
    if AUTOZERO==3,
        baseline=PEAKHEIGHTS(1);
        Heights=PEAKHEIGHTS(2:1+NumPeaks);
        model=Heights'*A+baseline;
    else
%          size(PEAKHEIGHTS) % error check
%          size(A)
        model=PEAKHEIGHTS'*A;
        Heights=PEAKHEIGHTS;
        baseline=0;
    end
end
if peakshape(1)==28, % polynomial;
    model=polynomial(xx,coeff);
    TrialParameters=PEAKHEIGHTS;
    Heights=ones(size(PEAKHEIGHTS));
end
% Compare trial model to data segment and compute the fit error
    MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy));
  % Take only the single fit that has the lowest MeanFitError
  if MeanFitError<LowestError, 
      if min(Heights)>=-BIPOLAR*10^100,  % Consider only fits with positive peak heights
        LowestError=MeanFitError;  % Assign LowestError to the lowest MeanFitError
        FitParameters=TrialParameters;  % Assign FitParameters to the fit with the lowest MeanFitError
        BestStart=newstart; % Assign BestStart to the start with the lowest MeanFitError
        height=Heights; % Assign height to the PEAKHEIGHTS with the lowest MeanFitError
        bestmodel=model; % Assign bestmodel to the model with the lowest MeanFitError
      end % if min(PEAKHEIGHTS)>0
  end % if MeanFitError<LowestError
%  ErrorVector(k)=MeanFitError;
end % for k (NumTrials)
    Rsquared=1-(norm(yy-bestmodel)./norm(yy-mean(yy)));
    SStot=sum((yy-mean(yy)).^2);
    SSres=sum((yy-bestmodel).^2);
    Rsquared=1-(SSres./SStot);
    GOF=[LowestError Rsquared];
% Uncomment following 4 lines to monitor trail fit starts and errors.
% StartMatrix=StartMatrix;
% ErrorVector=ErrorVector;
% matrix=[StartMatrix ErrorVector']
% std(StartMatrix)
% Construct model from best-fit parameters
AA=zeros(NumPeaks,600);
xxx=linspace(min(xx),max(xx),600);
% minxxx=min(xxx)
% xxx=linspace(min(xx)-length(xx),max(xx)+length(xx),200);
for m=1:NumPeaks,
   switch peakshape(1)
    case 1
        AA(m,:)=gaussian(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 2
        AA(m,:)=lorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 3
        AA(m,:)=logistic(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 4
        AA(m,:)=pearson(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 5
        AA(m,:)=expgaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra*length(xxx)./length(xx))';
    case 6
        AA(m,:)=gaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1));
    case 7
        AA(m,:)=lorentzian(xxx,FitParameters(m),FitParameters(NumPeaks+1));
    case 8
        AA(m,:)=expgaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1),-extra*length(xxx)./length(xx))';
    case 9
        AA(m,:)=exppulse(xxx,FitParameters(2*m-1),FitParameters(2*m));  
    case 10
        AA(m,:)=upsigmoid(xxx,FitParameters(2*m-1),FitParameters(2*m));   
    case 11
        AA(m,:)=gaussian(xxx,FitParameters(m),FIXEDPARAMETERS(m));
    case 12
        AA(m,:)=lorentzian(xxx,FitParameters(m),FIXEDPARAMETERS(m));
    case 13
        AA(m,:)=GL(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
    case 14
        AA(m,:)=BiGaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 15
        AA(m,:)=BWF(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 16
        AA(m,:)=gaussian(xxx,FIXEDPOSITIONS(m),FitParameters(m));
    case 17
        AA(m,:)=lorentzian(xxx,FIXEDPOSITIONS(m),FitParameters(m));
    case 18
        AA(m,:)=explorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra*length(xxx)./length(xx))';
    case 19
        AA(m,:)=alphafunction(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 20
        AA(m,:)=voigt(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
    case 21
        AA(m,:)=triangular(xxx,FitParameters(2*m-1),FitParameters(2*m));
    case 22
        AA(m,:)=peakfunction(shapesvector(m),xxx,FitParameters(2*m-1),FitParameters(2*m),extra(m));        
    case 23
        AA(m,:)=downsigmoid(xxx,FitParameters(2*m-1),FitParameters(2*m));  
    case 24
        AA(m,:)=nbinpdf(xxx,FitParameters(2*m-1),FitParameters(2*m));    
    case 25
        AA(m,:)=lognormal(xxx,FitParameters(2*m-1),FitParameters(2*m));    
    case 26
        AA(m,:)=linslope(xxx,FitParameters(2*m-1),FitParameters(2*m));   
    case 27
        AA(m,:)=d1gauss(xxx,FitParameters(2*m-1),FitParameters(2*m));  
    case 28
        AA(m,:)=polynomial(xxx,coeff);
    case 29
    case 30
        AA(m,:)=voigt(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m));        
    case 31
        AA(m,:)=expgaussian(xxx,FitParameters(3*m-2),FitParameters(3*m-1),-FitParameters(3*m)*length(xxx)./length(xx));        
    case 32
        AA(m,:)=pearson(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m));        
    case 33
        AA(m,:)=GL(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m)); 
    case 34
                  width(m)=abs(FIXEDPARAMETERS(m));
%                 gD(m)=width(m);
%                 gL(m)=extra.*gD(m);
%                 width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 +
%                 gD(m).^2));
        AA(m,:)=Voigt(xxx,FitParameters(m),width(m),extra);
    case 35
        AA(m,:)=GL(xxx,FitParameters(m),FIXEDPARAMETERS(m),extra);
    case 36
        AA(m,:)=expgaussian(xxx,FitParameters(m),FIXEDPARAMETERS(m),extra);    
    case 37
        AA(m,:)=pearson(xxx,FitParameters(m),FIXEDPARAMETERS(m),extra);
    otherwise
  end % switch
end % for NumPeaks

% Multiplies each row by the corresponding amplitude and adds them up
if peakshape(1)==29, % Segmented linear
    mmodel=segmented(xx,yy,PEAKHEIGHTS);
    baseline=0;
else
    heightsize=size(height');
    AAsize=size(AA);
    if heightsize(2)==AAsize(1),
        mmodel=height'*AA+baseline;
    else
        mmodel=height*AA+baseline;
    end
end
% Top half of the figure shows original signal and the fitted model.
if plots,
    subplot(2,1,1);plot(xx+xoffset,yy,'b.'); % Plot the original signal in blue dots
    hold on
end
if peakshape(1)==28, % Polynomial
     yi=polynomial(xxx,coeff);
else
    for m=1:NumPeaks,
        if plots, plot(xxx+xoffset,height(m)*AA(m,:)+baseline,'g'),end  % Plot the individual component peaks in green lines
        area(m)=trapz(xxx+xoffset,height(m)*AA(m,:)); % Compute the area of each component peak using trapezoidal method
        yi(m,:)=height(m)*AA(m,:); % Place y values of individual model peaks into matrix yi
    end
end
xi=xxx+xoffset; % Place the x-values of the individual model peaks into xi

if plots,
    % Mark starting peak positions with vertical dashed magenta lines
    if peakshape(1)==16||peakshape(1)==17
    else
        if peakshape(1)==29, % Segmented linear
            subplot(2,1,1);plot([PEAKHEIGHTS' PEAKHEIGHTS'],[0 max(yy)],'m--')
        else
            for marker=1:NumPeaks,
                markx=BestStart((2*marker)-1);
                subplot(2,1,1);plot([markx+xoffset markx+xoffset],[0 max(yy)],'m--')
            end % for
        end
    end % if peakshape

    % Plot the total model (sum of component peaks) in red lines
    if peakshape(1)==29, % Segmented linear
        mmodel=segmented(xx,yy,PEAKHEIGHTS);
       plot(xx+xoffset,mmodel,'r');  
    else
       plot(xxx+xoffset,mmodel,'r');  
    end
    hold off;
    lyy=min(yy);
    uyy=max(yy)+(max(yy)-min(yy))/10;
    if BIPOLAR,
        axis([min(xx+xoffset) max(xx+xoffset) lyy uyy]);
        ylabel('+ - mode')
    else
        axis([min(xx+xoffset) max(xx+xoffset) lyy uyy]);
        ylabel('+ mode')
    end
    switch AUTOZERO,
        case 0
            title(['peakfit.m Version 7.7   No baseline correction'])
        case 1
            title(['peakfit.m Version 7.7   Linear baseline subtraction'])
        case 2
            title(['peakfit.m Version 7.7   Quadratic subtraction baseline'])
        case 3
            title(['peakfit.m Version 7.7  Flat baseline correction'])
    end
 
    switch peakshape(1)
        case {4,20,34,37}
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      Shape Constant = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000) ] )
        case {5,8,18,36}
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      Time Constant = ' num2str(extra)   '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000)  ] )
        case 13
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      % Gaussian = ' num2str(extra)   '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000)  ] )
        case {14,15,22,35}
            xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      extra = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000) ] )
        case 28
            xlabel(['Shape = ' ShapeString '      Order = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '%  R2 = ' num2str(round(1000*LowestError)/1000) ] )
        otherwise
            if peakshape(1)==29, % Segmented linear
                xlabel(['Breakpoints = ' num2str(NumPeaks) '     Shape = ' ShapeString  '     Error = ' num2str(round(1000*LowestError)/1000) '%  R2 = ' num2str(round(100000*Rsquared)/100000) ] )
            else
                xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH)  '     Error = ' num2str(round(1000*LowestError)/1000) '%  R2 = ' num2str(round(100000*Rsquared)/100000) ] )
            end % if peakshape(1)==29
    end % switch peakshape(1)

    % Bottom half of the figure shows the residuals and displays RMS error
    % between original signal and model
    residual=yy-bestmodel;
    subplot(2,1,2);plot(xx+xoffset,residual,'r.')
    axis([min(xx)+xoffset max(xx)+xoffset min(residual) max(residual)]);
    xlabel('Residual Plot')
    if NumTrials>1,
       title(['Best of ' num2str(NumTrials) ' fits'])
    else
       title(['Single fit'])
    end
end % if plots

% Put results into a matrix FitResults, one row for each peak, showing peak index number,
% position, amplitude, and width.
FitResults=zeros(NumPeaks,6);
%  FitParameters=FitParameters
switch peakshape(1),
    case {6,7,8}, % equal-width peak models only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
            end
        end
    case {11,12,34,35,36,37}, % Fixed-width shapes only
        for m=1:NumPeaks,
            width(m)=abs(FitParameters(m));
            if peakshape==34,
                gD(m)=width(m);
                gL(m)=extra.*gD(m);
                width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2));
            end
            if m==1,
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS(m) area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS(m) area(m)]];
            end
        end
    case {16,17}, % Fixed-position shapes only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)];
            else
                FitResults=[FitResults ; [round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)]];
            end
        end
    case 28,   % Simple polynomial fit
        FitResults=PEAKHEIGHTS;
    case 29, % Segmented linear fit
        FitResults=PEAKHEIGHTS;
    case {30,31,32,33} % Special case of shapes with 3 iterated variables
        for m=1:NumPeaks,
            width(m)=abs(FitParameters(3*m-1));
            if peakshape==30,
                gD(m)=width(m);
                gL(m)=FitParameters(3*m).*gD(m);
                width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2));
            end
            if m==1,
                FitResults=[round(m) FitParameters(3*m-2) height(m) width(m) area(m) FitParameters(3*m)];
            else
                FitResults=[FitResults ; [round(m) FitParameters(3*m-2) height(m) width(m) area(m)] FitParameters(3*m)];
            end
        end
    otherwise % Normal shapes with 2 iterated variables
        for m=1:NumPeaks,
            width(m)=abs(FitParameters(2*m));
            if peakshape==20,
                gD=width(m);
                gL=extra.*gD;
                width(m) = 2.*(0.5346*gL + sqrt(0.2166*gL.^2 + gD.^2));
            end
            if m==1,
                FitResults=[round(m) FitParameters(2*m-1)+xoffset height(m) width(m) area(m)];
            else
                FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) width(m) area(m)]];
            end % if m==1
            
        end % for m=1:NumPeaks,
end % switch peakshape(1)
  
% Display Fit Results on lower graph
if plots,
    % Display Fit Results on lower  graph
    subplot(2,1,2);
    minxx=min(xx);
    startx=xoffset+min(xx)+(max(xx)-min(xx))./20;
    dxx=(max(xx)-min(xx))./10;
    dyy=((max(residual)-min(residual))./10);
    starty=max(residual)-dyy;
    FigureSize=get(gcf,'Position');
    switch peakshape(1)
        case {9,19,10,23}  % Pulse and sigmoid shapes only
            text(startx,starty+dyy/2,['Peak #          tau1           Height           tau2             Area'] );
        case 28, % Polynomial
            text(startx,starty+dyy/2,['Polynomial coefficients'] );
        case 29 % Segmented linear
             text(startx,starty+dyy/2,['x-axis breakpoints'] );
        case {30,31,32,33} % Special case of shapes with 3 iterated variables
            text(startx,starty+dyy/2,['Peak #          Position        Height         Width             Area       Shape factor'] );            
        otherwise
            text(startx,starty+dyy/2,['Peak #          Position        Height         Width             Area '] );
    end
    % Display FitResults using sprintf
    if peakshape(1)==28||peakshape(1)==29, % Polynomial or segmented linear
        for number=1:length(FitResults),
            column=1;
            itemstring=sprintf('%0.4g',FitResults(number));
            xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)));
            yposition=starty-number.*dyy.*(400./FigureSize(4));
            text(xposition,yposition,['                ' itemstring]);
        end
    else
        for peaknumber=1:NumPeaks,
            for column=1:5,
                itemstring=sprintf('%0.4g',FitResults(peaknumber,column));
                xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)));
                yposition=starty-peaknumber.*dyy.*(400./FigureSize(4));
                text(xposition,yposition,itemstring);
            end
        end
        xposition=startx;
        yposition=starty-(peaknumber+1).*dyy.*(400./FigureSize(4));
        if AUTOZERO==3,
            text(xposition,yposition,[ 'Baseline= ' num2str(baseline) ]);
        end % if AUTOZERO
    end % if peakshape(1)
    if peakshape(1)==30 || peakshape(1)==31 || peakshape(1)==32 || peakshape(1)==33,
        for peaknumber=1:NumPeaks,
            column=6;
            itemstring=sprintf('%0.4g',FitParameters(3*peaknumber));
            xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)));
            yposition=starty-peaknumber.*dyy.*(400./FigureSize(4));
            text(xposition,yposition,itemstring);
        end
    end
end % if plots

if NumArgOut==8,
    if plots,disp('Computing bootstrap sampling statistics.....'),end
    BootstrapResultsMatrix=zeros(6,100,NumPeaks);
    BootstrapErrorMatrix=zeros(1,100,NumPeaks);
    clear bx by
    tic;
    for trial=1:100,
        n=1;
        bx=xx;
        by=yy;
        while n<length(xx)-1,
            if rand>.5,
                bx(n)=xx(n+1);
                by(n)=yy(n+1);
            end
            n=n+1;
        end
        bx=bx+xoffset;
        % fitpeaks called
        [FitResults,BootFitError]=fitpeaks(bx,by,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,FIXEDPARAMETERS,shapesvector);
        for peak=1:NumPeaks,
            switch peakshape(1)
                case {30,31,32,33}
                    BootstrapResultsMatrix(1:6,trial,peak)=FitResults(peak,1:6);
                otherwise
                    BootstrapResultsMatrix(1:5,trial,peak)=FitResults(peak,1:5);
            end
            BootstrapErrorMatrix(:,trial,peak)=BootFitError;
        end
    end
    if plots,toc;end
    for peak=1:NumPeaks,
        if plots,
            disp(' ')
            disp(['Peak #',num2str(peak) '         Position    Height       Width       Area      Shape Factor']);
        end % if plots
        BootstrapMean=mean(real(BootstrapResultsMatrix(:,:,peak)'));
        BootstrapSTD=std(BootstrapResultsMatrix(:,:,peak)');
        BootstrapIQR=iqr(BootstrapResultsMatrix(:,:,peak)');
        PercentRSD=100.*BootstrapSTD./BootstrapMean;
        PercentIQR=100.*BootstrapIQR./BootstrapMean;
        BootstrapMean=BootstrapMean(2:6);
        BootstrapSTD=BootstrapSTD(2:6);
        BootstrapIQR=BootstrapIQR(2:6);
        PercentRSD=PercentRSD(2:6);
        PercentIQR=PercentIQR(2:6);
        if plots,
            disp(['Bootstrap Mean: ', num2str(BootstrapMean)])
            disp(['Bootstrap STD:  ', num2str(BootstrapSTD)])
            disp(['Bootstrap IQR:  ', num2str(BootstrapIQR)])
            disp(['Percent RSD:    ', num2str(PercentRSD)])
            disp(['Percent IQR:    ', num2str(PercentIQR)])
        end % if plots
        BootResults(peak,:)=[BootstrapMean BootstrapSTD PercentRSD BootstrapIQR PercentIQR];
    end % peak=1:NumPeaks,
end % if NumArgOut==8,
if AUTOZERO==3;
else
    baseline=bkgcoef;
end
% ----------------------------------------------------------------------
function [FitResults,LowestError]=fitpeaks(xx,yy,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,fixedparameters,shapesvector)
% Based on peakfit Version 3: June, 2012. 
global PEAKHEIGHTS FIXEDPARAMETERS BIPOLAR MINWIDTH coeff
format short g
format compact
warning off all
FIXEDPARAMETERS=fixedparameters;
if start==0;start=calcstart(xx,NumPeaks,xoffset);end
PEAKHEIGHTS=zeros(1,NumPeaks);
n=length(xx);
newstart=start;
coeff=0;
LOGPLOT=0;

% Perform peak fitting for selected peak shape using fminsearch function
options = optimset('TolX',.00001,'TolFun',.00001,'Display','off','MaxFunEvals',1000 );
LowestError=1000; % or any big number greater than largest error expected
FitParameters=zeros(1,NumPeaks.*2); 
BestStart=zeros(1,NumPeaks.*2); 
height=zeros(1,NumPeaks); 
bestmodel=zeros(size(yy));
for k=1:NumTrials,
    % StartVector=newstart
    switch peakshape(1)
        case 1
            TrialParameters=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 2
            TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 3
            TrialParameters=fminsearch(@(lambda)(fitlogistic(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 4
            TrialParameters=fminsearch(@(lambda)(fitpearson(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 5
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexpgaussian(lambda,zxx,zyy,-extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 6
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitewgaussian(lambda,xx,yy)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 7
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitewlorentzian(lambda,xx,yy)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 8
            cwnewstart(1)=newstart(1);
            for pc=2:NumPeaks,
                cwnewstart(pc)=newstart(2.*pc-1);
            end
            cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
            TrialParameters=fminsearch(@(lambda)(fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(NumPeaks+1)<MINWIDTH,
                    TrialParameters(NumPeaks+1)=MINWIDTH;
                end
            end
        case 9
            TrialParameters=fminsearch(@(lambda)(fitexppulse(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 10
            TrialParameters=fminsearch(@(lambda)(fitupsigmoid(lambda,xx,yy)),newstar,optionst);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 11
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFWGaussian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 12
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFWLorentzian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 13
            TrialParameters=fminsearch(@(lambda)(fitGL(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 14
            TrialParameters=fminsearch(@(lambda)(fitBiGaussian(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 15
            TrialParameters=fminsearch(@(lambda)(fitBWF(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 16
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(Peak)<MINWIDTH,
                    TrialParameters(Peak)=MINWIDTH;
                end
            end
        case 17
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(Peak)<MINWIDTH,
                    TrialParameters(Peak)=MINWIDTH;
                end
            end
        case 18
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexplorentzian(lambda,zxx,zyy,-extra)),newstart,options);
        case 19
            TrialParameters=fminsearch(@(lambda)(alphafunction(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 20
            TrialParameters=fminsearch(@(lambda)(fitvoigt(lambda,xx,yy,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 21
            TrialParameters=fminsearch(@(lambda)(fittriangular(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 22
            TrialParameters=fminsearch(@(lambda)(fitmultiple(lambda,xx,yy,shapesvector,extra)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 23
            TrialParameters=fminsearch(@(lambda)(fitdownsigmoid(lambda,xx,yy)),newstart,optionst);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 24
            TrialParameters=fminsearch(@(lambda)(fitnbinpdf(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 25
            TrialParameters=fminsearch(@(lambda)(fitlognpdf(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 26
            TrialParameters=fminsearch(@(lambda)(fitlinslope(lambda,xx,yy)),polyfit(xx,yy,1),options);
        coeff=TrialParameters;
        case 27
            TrialParameters=fminsearch(@(lambda)(fitd1gauss(lambda,xx,yy)),newstart,options);
            for Peak=1:NumPeaks;
                if TrialParameters(2*Peak)<MINWIDTH,
                    TrialParameters(2*Peak)=MINWIDTH;
                end
            end
        case 28
            TrialParameters=fitpolynomial(xx,yy,extra);
        case 29
            TrialParameters=fminsearch(@(lambda)(fitsegmented(lambda,xx,yy)),newstart,options);
        case 30
            TrialParameters=fminsearch(@(lambda)(fitvoigtv(lambda,xx,yy)),newstart);
        case 31
            zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
            zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
            TrialParameters=fminsearch(@(lambda)(fitexpgaussianv(lambda,zxx,zyy)),newstart);
        case 32
            TrialParameters=fminsearch(@(lambda)(fitpearsonv(lambda,xx,yy)),newstart);
        case 33
            TrialParameters=fminsearch(@(lambda)(fitGLv(lambda,xx,yy)),newstart);
        case 34
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(fitFWVoigt(lambda,xx,yy,extra)),fixedstart,options);
        case 35
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(fitFWGL(lambda,xx,yy,extra)),fixedstart,options);
        case 36
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(fitFWExpGaussian(lambda,xx,yy,extra)),fixedstart,options);
        case 37
            fixedstart=[];
            for pc=1:NumPeaks,
                fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
            end
            TrialParameters=fminsearch(@(lambda)(fitFWPearson(lambda,xx,yy,extra)),fixedstart,options);
        otherwise
    end % switch peakshape
    
for peaks=1:NumPeaks,
     peakindex=2*peaks-1;
     newstart(peakindex)=start(peakindex)-xoffset;
end

    % Construct model from Trial parameters
    A=zeros(NumPeaks,n);
    for m=1:NumPeaks,
        switch peakshape(1)
            case 1
                A(m,:)=gaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 2
                A(m,:)=lorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 3
                A(m,:)=logistic(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 4
                A(m,:)=pearson(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 5
                A(m,:)=expgaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
            case 6
                A(m,:)=gaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
            case 7
                A(m,:)=lorentzian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
            case 8
                A(m,:)=expgaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1),-extra)';
            case 9
                A(m,:)=exppulse(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 10
                A(m,:)=upsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 11
                A(m,:)=gaussian(xx,TrialParameters(m),FIXEDPARAMETERS(m));
            case 12
                A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDPARAMETERS(m));
            case 13
                A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 14
                A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 15
                A(m,:)=BWF(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 16
                A(m,:)=gaussian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
            case 17
                A(m,:)=lorentzian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
            case 18
                A(m,:)=explorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
            case 19
                A(m,:)=alphafunction(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 20
                A(m,:)=voigt(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
            case 21
                A(m,:)=triangular(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 22
                A(m,:)=peakfunction(shapesvector(m),xx,TrialParameters(2*m-1),TrialParameters(2*m),extra(m));
            case 23
                A(m,:)=downsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));      
            case 24
                A(m,:)=nbinpdf(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 25
                A(m,:)=lognormal(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 26
                A(m,:)=linslope(xx,TrialParameters(2*m-1),TrialParameters(2*m));
            case 27
                A(m,:)=d1gauss(xx,TrialParameters(2*m-1),TrialParameters(2*m));       
            case 28
                A(m,:)=polynomial(xx,TrialParameters(2*m-1),TrialParameters(2*m));       
            case 29
                A(m,:)=segmented(xx,yy,PEAKHEIGHTS);
            case 30
                A(m,:)=voigt(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
            case 31
                A(m,:)=expgaussian(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
            case 32
                A(m,:)=pearson(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
            case 33
                A(m,:)=GL(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));
            case 34
                width(m)=abs(FIXEDPARAMETERS(m));

%                 gD(m)=width(m);
%                 gL(m)=extra.*gD(m);
%                 width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2))

            A(m,:)=Voigt(xx,TrialParameters(m), width(m),extra);
            case 35
                A(m,:)=GL(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);
            case 36
                A(m,:)=expgaussian(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);
            case 37
                A(m,:)=pearson(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);
        end % switch
    end % for
    
    % Multiplies each row by the corresponding amplitude and adds them up
    if peakshape(1)==29, % Segmented linear
        model=segmented(xx,yy,PEAKHEIGHTS);
        TrialParameters=coeff;
        Heights=ones(size(coeff));
    else
        if AUTOZERO==3,
            baseline=PEAKHEIGHTS(1);
            Heights=PEAKHEIGHTS(2:1+NumPeaks);
            model=Heights'*A+baseline;
        else
            model=PEAKHEIGHTS'*A;
            Heights=PEAKHEIGHTS;
            baseline=0;
        end
    end
    
    % Compare trial model to data segment and compute the fit error
    MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy));
    % Take only the single fit that has the lowest MeanFitError
    if MeanFitError<LowestError,
        if min(Heights)>=-BIPOLAR*10^100,  % Consider only fits with positive peak heights
            LowestError=MeanFitError;  % Assign LowestError to the lowest MeanFitError
            FitParameters=TrialParameters;  % Assign FitParameters to the fit with the lowest MeanFitError
            height=Heights; % Assign height to the PEAKHEIGHTS with the lowest MeanFitError
        end % if min(PEAKHEIGHTS)>0
    end % if MeanFitError<LowestError
end % for k (NumTrials)
    Rsquared=1-(norm(yy-bestmodel)./norm(yy-mean(yy)));
    SStot=sum((yy-mean(yy)).^2);
    SSres=sum((yy-bestmodel).^2);
    Rsquared=1-(SSres./SStot);
    GOF=[LowestError Rsquared];
for m=1:NumPeaks,
    area(m)=trapz(xx+xoffset,height(m)*A(m,:)); % Compute the area of each component peak using trapezoidal method
end

% Put results into a matrix FitResults, one row for each peak, showing peak index number,
% position, amplitude, and width.
FitResults=zeros(NumPeaks,6);
switch peakshape(1),
    case {6,7,8}, % equal-width peak models only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
            end
        end
    case {11,12,35,36,37}, % Fixed-width shapes only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS(m) area(m)]];
            else
                FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS(m) area(m)]];
            end
        end

    case {16,17}, % Fixed-position shapes only
        for m=1:NumPeaks,
            if m==1,
                FitResults=[round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)];
            else
                FitResults=[FitResults ; [round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)]];
            end
        end
    case 28,   % Simple polynomial fit
        FitResults=PEAKHEIGHTS;
    case 29, % Segmented linear fit
        FitResults=PEAKHEIGHTS;
    case {30,31,32,33} % Special case of shapes with 3 iterated variables
        for m=1:NumPeaks,
            width(m)=abs(FitParameters(3*m-1));
            if peakshape==30,
                gD(m)=width(m);
                gL(m)=FitParameters(3*m).*gD(m);
                width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2));
            end
            if m==1,
                FitResults=[round(m) FitParameters(3*m-2) height(m) width(m) area(m) FitParameters(3*m)];
            else
                FitResults=[FitResults ; [round(m) FitParameters(3*m-2) height(m) width(m) area(m) FitParameters(3*m)]];
            end
        end
    otherwise % Normal shapes with 2 iterated variables
        for m=1:NumPeaks,
            width(m)=abs(FitParameters(2*m));
            if peakshape==20,
                gD=width(m);
                gL=extra.*gD;
                width(m) = 2.*(0.5346*gL + sqrt(0.2166*gL.^2 + gD.^2));
            end
            if m==1,
                FitResults=[round(m) FitParameters(2*m-1)+xoffset height(m) width(m) area(m)];
            else
                FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) width(m) area(m)]];
            end % if m==1
        end % for m=1:NumPeaks,
end % switch peakshape(1)
if peakshape==34,
        DW=2*(0.5346*a*1.2772 + sqrt(0.2166*a*1.2772.^2 + 1.2772.^2))
end
% ----------------------------------------------------------------------
function start=calcstart(xx,NumPeaks,xoffset)
  n=max(xx)-min(xx);
  start=[];
  startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
  for marker=1:NumPeaks,
      markx=startpos(marker);
      start=[start markx n/ (3.*NumPeaks)];
  end % for marker
% ----------------------------------------------------------------------
function err = fitgaussian(lambda,t,y)
% Fitting function for a Gaussian band signal.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
numpeaks=round(length(lambda)/2);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
%    if lambda(2*j)<MINWIDTH,lambda(2*j)=MINWIDTH;end
    A(:,j) = gaussian(t,lambda(2*j-1),lambda(2*j))';
end 
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitewgaussian(lambda,t,y)
% Fitting function for a Gaussian band signal with equal peak widths.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(j),lambda(numpeaks+1))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = FitFWGaussian(lambda,t,y)
%	Fitting function for a fixed width Gaussian
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,lambda(j),FIXEDPARAMETERS(j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = FitFPGaussian(lambda,t,y)
%	Fitting function for fixed-position Gaussians
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = gaussian(t,FIXEDPARAMETERS(j), lambda(j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = FitFPLorentzian(lambda,t,y)
%	Fitting function for fixed-position Lorentzians
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,FIXEDPARAMETERS(j), lambda(j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function err = FitFWLorentzian(lambda,t,y)
%	Fitting function for fixed width Lorentzian
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,lambda(j),FIXEDPARAMETERS(j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitewlorentzian(lambda,t,y)
% Fitting function for a Lorentzian band signal with equal peak widths.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = lorentzian(t,lambda(j),lambda(numpeaks+1))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitlorentzian(lambda,t,y)
%	Fitting function for single lorentzian, lambda(1)=position, lambda(2)=width
%	Fitgauss assumes a lorentzian function 
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = lorentzian(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = lorentzian(x,position,width)
% lorentzian(x,position,width) Lorentzian function.
% where x may be scalar, vector, or matrix
% position and width scalar
% T. C. O'Haver, 1988
% Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);
% ----------------------------------------------------------------------
function err = fitlogistic(lambda,t,y)
%	Fitting function for logistic, lambda(1)=position, lambda(2)=width
%	between the data and the values computed by the current
%	function of lambda.  Fitlogistic assumes a logistic function 
%  T. C. O'Haver, May 2006
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = logistic(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = logistic(x,pos,wid)
% logistic function.  pos=position; wid=half-width (both scalar)
% logistic(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991 
n = exp(-((x-pos)/(.477.*wid)) .^2);
g = (2.*n)./(1+n);
% ----------------------------------------------------------------------
function err = fittriangular(lambda,t,y)
%	Fitting function for triangular, lambda(1)=position, lambda(2)=width
%	between the data and the values computed by the current
%	function of lambda.  Fittriangular assumes a triangular function 
%  T. C. O'Haver, May 2006
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = triangular(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = triangular(x,pos,wid)
%triangle function.  pos=position; wid=half-width (both scalar)
%trianglar(x,pos,wid), where x may be scalar or vector,
%pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991
% Example
% x=[0:.1:10];plot(x,trianglar(x,5.5,2.3),'.')
g=1-(1./wid) .*abs(x-pos);
for i=1:length(x),  
if g(i)<0,g(i)=0;end
end
% ----------------------------------------------------------------------
function err = fitpearson(lambda,t,y,shapeconstant)
%   Fitting functions for a Pearson 7 band signal.
% T. C. O'Haver (toh@umd.edu),   Version 1.3, October 23, 2006.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = pearson(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitpearsonv(lambda,t,y)
% Fitting functions for pearson function with independently variable
% percent Gaussian
% T. C. O'Haver (toh@umd.edu), 2015.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
for j = 1:length(lambda)/3,
    A(:,j) = pearson(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitFWPearson(lambda,t,y,shapeconstant)
%	Fitting function for a fixed width Pearson7
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = pearson(t,lambda(j),FIXEDPARAMETERS(j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = pearson(x,pos,wid,m)
% Pearson VII function. 
% g = pearson(x,pos,wid,m) where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% m=some number
%  T. C. O'Haver, 1990  
g=ones(size(x))./(1+((x-pos)./((0.5.^(2/m)).*wid)).^2).^m;
% ----------------------------------------------------------------------
function err = fitexpgaussian(lambda,t,y,timeconstant)
%   Fitting functions for a exponentially-broadened Gaussian band signal.
%  T. C. O'Haver, October 23, 2006.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = expgaussian(t,lambda(2*j-1),lambda(2*j),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitexplorentzian(lambda,t,y,timeconstant)
%   Fitting functions for a exponentially-broadened lorentzian band signal.
%  T. C. O'Haver, 2013.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = explorentzian(t,lambda(2*j-1),lambda(2*j),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitexpewgaussian(lambda,t,y,timeconstant)
% Fitting function for exponentially-broadened Gaussian bands with equal peak widths.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
numpeaks=round(length(lambda)-1);
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = expgaussian(t,lambda(j),lambda(numpeaks+1),timeconstant);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitexpgaussianv(lambda,t,y)
% Fitting functions for  exponentially-broadened Gaussians with
% independently variable time constants
% T. C. O'Haver (toh@umd.edu), 2015.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
for j = 1:length(lambda)/3,
    A(:,j) = expgaussian(t,lambda(3*j-2),lambda(3*j-1),-lambda(3*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitFWExpGaussian(lambda,t,y,shapeconstant)
%	Fitting function for a fixed width Exponentially-broadened gaussian
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = expgaussian(t,lambda(j),FIXEDPARAMETERS(j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = expgaussian(x,pos,wid,timeconstant)
%  Exponentially-broadened gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 2006
g = exp(-((x-pos)./(0.6005615.*wid)) .^2);
g = ExpBroaden(g',timeconstant);
% ----------------------------------------------------------------------
function g = explorentzian(x,pos,wid,timeconstant)
%  Exponentially-broadened lorentzian(x,pos,wid) = lorentzian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 2013
g = ones(size(x))./(1+((x-pos)./(0.5.*wid)).^2);
g = ExpBroaden(g',timeconstant);
% ----------------------------------------------------------------------
function yb = ExpBroaden(y,t)
% ExpBroaden(y,t) zero pads y and convolutes result by an exponential decay
% of time constant t by multiplying Fourier transforms and inverse
% transforming the result.
hly=round(length(y)./2);
ey=[y(1).*ones(1,hly)';y;y(length(y)).*ones(1,hly)'];
% figure(2);plot(ey);figure(1);
fy=fft(ey);
a=exp(-(1:length(fy))./t);
fa=fft(a);
fy1=fy.*fa';
ybz=real(ifft(fy1))./sum(a);
yb=ybz(hly+2:length(ybz)-hly+1);
% ----------------------------------------------------------------------
function err = fitexppulse(tau,x,y)
% Iterative fit of the sum of exponential pulses
% of the form Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = exppulse(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = exppulse(x,t1,t2)
% Exponential pulse of the form 
% g = (x-spoint)./pos.*exp(1-(x-spoint)./pos);
e=(x-t1)./t2;
p = 4*exp(-e).*(1-exp(-e));
p=p .* (p>0);
g = p';
% ----------------------------------------------------------------------
function err = fitalphafunction(tau,x,y)
% Iterative fit of the sum of alpha funciton
% of the form Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = alphafunction(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = alphafunction(x,pos,spoint)
% alpha function.  pos=position; wid=half-width (both scalar)
% alphafunction(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% Taekyung Kwon, July 2013  
g = (x-spoint)./pos.*exp(1-(x-spoint)./pos);
for m=1:length(x);if g(m)<0;g(m)=0;end;end
% ----------------------------------------------------------------------
function err = fitdownsigmoid(tau,x,y)
% Fitting function for iterative fit to the sum of
% downward moving sigmiods 
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = downsigmoid(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitupsigmoid(tau,x,y)
% Fitting function for iterative fit to the sum of
% upwards moving sigmiods
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = upsigmoid(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g=downsigmoid(x,t1,t2)
 % down step sigmoid
g=.5-.5*erf(real((x-t1)/sqrt(2*t2)));
% ----------------------------------------------------------------------
function g=upsigmoid(x,t1,t2)
% up step sigmoid
g=1/2 + 1/2* erf(real((x-t1)/sqrt(2*t2))); 
% ----------------------------------------------------------------------
function err = fitGL(lambda,t,y,shapeconstant)
%   Fitting functions for Gaussian/Lorentzian blend.
% T. C. O'Haver (toh@umd.edu), 2012.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = GL(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitFWGL(lambda,t,y,shapeconstant)
%	Fitting function for a fixed width Gaussian/Lorentzian blend
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = GL(t,lambda(j),FIXEDPARAMETERS(j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitGLv(lambda,t,y)
% Fitting functions for Gaussian/Lorentzian blend function with
% independently variable percent Gaussian
% T. C. O'Haver (toh@umd.edu), 2015.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
for j = 1:length(lambda)/3,
    A(:,j) = GL(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = GL(x,pos,wid,m)
% Gaussian/Lorentzian blend. m = percent Gaussian character
% pos=position; wid=half-width
% m = percent Gaussian character.
%  T. C. O'Haver, 2012
% sizex=size(x)
% sizepos=size(pos)
% sizewid=size(wid)
% sizem=size(m)
g=2.*((m/100).*gaussian(x,pos,wid)+(1-(m(1)/100)).*lorentzian(x,pos,wid))/2;
% ----------------------------------------------------------------------
function err = fitvoigt(lambda,t,y,shapeconstant)
% Fitting functions for Voigt profile function
% T. C. O'Haver (toh@umd.edu), 2013.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = voigt(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitFWVoigt(lambda,t,y,shapeconstant)
%	Fitting function for a fixed width Voigt
global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
numpeaks=round(length(lambda));
A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    A(:,j) = voigt(t,lambda(j),FIXEDPARAMETERS(j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitvoigtv(lambda,t,y)
% Fitting functions for Voigt profile function with independently variable
% alphas
% T. C. O'Haver (toh@umd.edu), 2015.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/3));
for j = 1:length(lambda)/3,
    A(:,j) = voigt(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g=voigt(xx,pos,gD,alpha)
% Voigt profile function. xx is the independent variable (energy,
% wavelength, etc), gD is the Doppler (Gaussian) width, and alpha is the
% shape constant (ratio of the Lorentzian width gL to the Doppler width gD.
% Based on Chong Tao's "Voigt lineshape spectrum simulation", 
% File ID: #26707
% alpha=alpha
gL=alpha.*gD;
gV = 0.5346*gL + sqrt(0.2166*gL.^2 + gD.^2);
x = gL/gV;
y = abs(xx-pos)/gV;
g = 1/(2*gV*(1.065 + 0.447*x + 0.058*x^2))*((1-x)*exp(-0.693.*y.^2) + (x./(1+y.^2)) + 0.016*(1-x)*x*(exp(-0.0841.*y.^2.25)-1./(1 + 0.021.*y.^2.25)));
g=g./max(g);
% ----------------------------------------------------------------------
function err = fitBiGaussian(lambda,t,y,shapeconstant)
%   Fitting functions for BiGaussian.
% T. C. O'Haver (toh@umd.edu),  2012.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = BiGaussian(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = BiGaussian(x,pos,wid,m)
% BiGaussian (different widths on leading edge and trailing edge).
% pos=position; wid=width 
% m determines shape; symmetrical if m=50.
%  T. C. O'Haver, 2012
lx=length(x);
hx=val2ind(x,pos);
g(1:hx)=gaussian(x(1:hx),pos,wid*(m/100));
g(hx+1:lx)=gaussian(x(hx+1:lx),pos,wid*(1-m/100));
% ----------------------------------------------------------------------
function err = fitBWF(lambda,t,y,shapeconstant)
%   Fitting function for Breit-Wigner-Fano.
% T. C. O'Haver (toh@umd.edu),  2014.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = BWF(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = BWF(x,pos,wid,m)
% BWF (Breit-Wigner-Fano) http://en.wikipedia.org/wiki/Fano_resonance
% pos=position; wid=width; m=Fano factor
%  T. C. O'Haver, 2014
y=((m*wid/2+x-pos).^2)./(((wid/2).^2)+(x-pos).^2);
% y=((1+(x-pos./(m.*wid))).^2)./(1+((x-pos)./wid).^2);
g=y./max(y);
% ----------------------------------------------------------------------
function err = fitnbinpdf(tau,x,y)
% Fitting function for iterative fit to the sum of
% Negative Binomial Distributions
% (http://www.mathworks.com/help/stats/negative-binomial-distribution.html)
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = nbinpdf(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function err = fitlognpdf(tau,x,y)
% Fitting function for iterative fit to the sum of
% Lognormal Distributions
% (http://www.mathworks.com/help/stats/lognormal-distribution.html)
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = lognormal(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = lognormal(x,pos,wid)
% lognormal function.  pos=position; wid=half-width (both scalar)
% lognormal(x,pos,wid), where x may be scalar, vector, or matrix
% pos=position; wid=half-width (both scalar)
% T. C. O'Haver, 1991  
g = exp(-(log(x/pos)/(0.01.*wid)) .^2);
% ----------------------------------------------------------------------
function err = fitsine(tau,x,y)
% Fitting function for iterative fit to the sum of
% sine waves (alpha test, NRFPT)
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    A(:,j) = sine(x,tau(2*j-1),tau(2*j));
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g=sine(x,f,phase) 
% Sine wave (alpha test)
g=sin(2*pi*f*(x+phase));
% ----------------------------------------------------------------------
function err = fitd1gauss(lambda,t,y)
%   Fitting functions for the first derivative of a Gaussian
%  T. C. O'Haver, 2014
global PEAKHEIGHTS AUTOZERO BIPOLAR
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = d1gauss(t,lambda(2*j-1),lambda(2*j))';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
err = norm(z-y');
% ----------------------------------------------------------------------
function y=d1gauss(x,p,w)
% First derivative of Gaussian (alpha test)
y=-(5.54518.*(x-p).*exp(-(2.77259.*(p-x).^2)./w^2))./w.^2;
y=y./max(y);
% ----------------------------------------------------------------------
function coeff = fitpolynomial(t,y,order)
coeff=polyfit(t,y,order);
% order=order
% coeff=coeff
% ----------------------------------------------------------------------
function y=polynomial(t,coeff)
y=polyval(coeff,t);
% ----------------------------------------------------------------------
function err = fitsegmented(lambda,t,y)
%   Fitting functions for articulated segmented linear
%  T. C. O'Haver, 2014
global LOGPLOT
breakpoints=[t(1) lambda max(t)];
z = segmented(t,y,breakpoints);
% lengthz=length(z);
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y);
end
% ----------------------------------------------------------------------
function yi=segmented(x,y,segs)
global PEAKHEIGHTS
clear yy
for n=1:length(segs)
  yind=val2ind(x,segs(n));
  yy(n)=y(yind(1));
end
yi=INTERP1(segs,yy,x);
PEAKHEIGHTS=segs;
% ----------------------------------------------------------------------
function err = fitlinslope(tau,x,y)
% Fitting function for iterative fit to linear function
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(x),round(length(tau)/2));
for j = 1:length(tau)/2,
    z = (x.*tau(2*j-1)+tau(2*j))';
    A(:,j) = z./max(z);
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function y=linslope(x,slope,intercept)
y=x.*slope+intercept;
% y=y./max(y);
% ----------------------------------------------------------------------
function b=iqr(a)
% b = IQR(a)  returns the interquartile range of the values in a.  For
%  vector input, b is the difference between the 75th and 25th percentiles
%  of a.  For matrix input, b is a row vector containing the interquartile
%  range of each column of a.
%  T. C. O'Haver, 2012
mina=min(a);
sizea=size(a);
NumCols=sizea(2);
for n=1:NumCols,b(:,n)=a(:,n)-mina(n);end
Sorteda=sort(b);
lx=length(Sorteda);
SecondQuartile=round(lx/4);
FourthQuartile=3*round(lx/4);
b=abs(Sorteda(FourthQuartile,:)-Sorteda(SecondQuartile,:));
% ----------------------------------------------------------------------
function err = fitmultiple(lambda,t,y,shapesvector,m)
% Fitting function for a multiple-shape band signal.
% The sequence of peak shapes are defined by the vector "shape".
% The vector "m" determines the shape of variable-shape peaks.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT coeff
numpeaks=round(length(lambda)/2);

A = zeros(length(t),numpeaks);
for j = 1:numpeaks,
    if shapesvector(j)==28,
        coeff=polyfit(t,y,m(j));
        A(:,j) = polyval(coeff,t);
    else
        A(:,j) = peakfunction(shapesvector(j),t,lambda(2*j-1),lambda(2*j),m(j))';
    end
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function p=peakfunction(shape,x,pos,wid,m,coeff)
% function that generates any of 20 peak types specified by number. 'shape'
% specifies the shape type of each peak in the signal: "peakshape" = 1-20.
% 1=Gaussian 2=Lorentzian, 3=logistic, 4=Pearson, 5=exponentionally
% broadened Gaussian; 9=exponential pulse, 10=up sigmoid,
% 13=Gaussian/Lorentzian blend; 14=BiGaussian, 15=Breit-Wigner-Fano (BWF) ,
% 18=exponentionally broadened Lorentzian; 19=alpha function; 20=Voigt
% profile; 21=triangular; 23=down sigmoid; 25=lognormal. "m" is required
% for variable-shape peaks only.
switch shape,
    case 1
        p=gaussian(x,pos,wid);
    case 2
        p=lorentzian(x,pos,wid);
    case 3
        p=logistic(x,pos,wid);
    case 4
        p=pearson(x,pos,wid,m);
    case 5
        p=expgaussian(x,pos,wid,m);
    case 6
        p=gaussian(x,pos,wid);
    case 7
        p=lorentzian(x,pos,wid);
    case 8
        p=expgaussian(x,pos,wid,m)';
    case 9
        p=exppulse(x,pos,wid);
    case 10
        p=upsigmoid(x,pos,wid);
    case 11
        p=gaussian(x,pos,wid);
    case 12
        p=lorentzian(x,pos,wid);
    case 13
        p=GL(x,pos,wid,m);
    case 14
        p=BiGaussian(x,pos,wid,m);
    case 15
        p=BWF(x,pos,wid,m);
    case 16
        p=gaussian(x,pos,wid);
    case 17
        p=lorentzian(x,pos,wid);
    case 18
        p=explorentzian(x,pos,wid,m)';
    case 19
        p=alphafunction(x,pos,wid);
    case 20
        p=voigt(x,pos,wid,m);
    case 21
        p=triangular(x,pos,wid);    
    case 23
        p=downsigmoid(x,pos,wid);
    case 25
        p=lognormal(x,pos,wid);
    case 26
        p=linslope(x,pos,wid);
    case 27
        p=d1gauss(x,pos,wid);
    case 28
        p=polynomial(x,coeff);
    otherwise
end % switch

function a=rmnan(a)
% Removes NaNs and Infs from vectors, replacing with nearest real numbers.
% Example:
%  >> v=[1 2 3 4 Inf 6 7 Inf  9];
%  >> rmnan(v)
%  ans =
%     1     2     3     4     4     6     7     7     9
la=length(a);
if isnan(a(1)) || isinf(a(1)),a(1)=0;end
for point=1:la,
    if isnan(a(point)) || isinf(a(point)),
        a(point)=a(point-1);
    end
end

function P=findpeaksG(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% function P=findpeaksG(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% Function to locate the positive peaks in a noisy x-y time series data
% set.  Detects peaks by looking for downward zero-crossings
% in the first derivative that exceed SlopeThreshold.
% Returns list (P) containing peak number and position,
% height, width, and area of each peak. Arguments "slopeThreshold",
% "ampThreshold" and "smoothwidth" control peak sensitivity.
% Higher values will neglect smaller features. "Smoothwidth" is
% the width of the smooth applied before peak detection; larger
% values ignore narrow peaks. If smoothwidth=0, no smoothing
% is performed. "Peakgroup" is the number points around the top
% part of the peak that are taken for measurement. If Peakgroup=0
% the local maximum is takes as the peak height and position.
% The argument "smoothtype" determines the smooth algorithm:
%   If smoothtype=1, rectangular (sliding-average or boxcar)
%   If smoothtype=2, triangular (2 passes of sliding-average)
%   If smoothtype=3, pseudo-Gaussian (3 passes of sliding-average)
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html and
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% (c) T.C. O'Haver, 1995, 2014.  Version 5, Last revised January, 2016  
%
%
% Examples:
% findpeaksG(0:.01:2,humps(0:.01:2),0,-1,5,5)
% x=[0:.01:50];findpeaks(x,cos(x),0,-1,5,5)
% x=[0:.01:5]';findpeaks(x,x.*sin(x.^2).^2,0,-1,5,5)
% x=[-10:.1:10];y=exp(-(x).^2);findpeaks(x,y,0.005,0.3,3,5,3);
% Find, measure, and plot noisy peak with unknown position
% x=[-10:.2:10];
% y=exp(-(x+5*randn()).^2)+.1.*randn(size(x));
% P=findpeaksG(x,y,0.003,0.5,7,9,3);
% xx=linspace(min(x),max(x));
% yy=P(3).*gaussian(xx,P(2),P(4));
% plot(x,y,'.',xx,yy)
%
% Related functions:
% findvalleys.m, findpeaksL.m, findpeaksb.m, findpeaksb3.m,
% findpeaksplot.m, peakstats.m, findpeaksnr.m, findpeaksGSS.m,
% findpeaksLSS.m, findpeaksfit.m, findsteps.m, findsquarepulse.m, idpeaks.m

% Copyright (c) 2013, 2014 Thomas C. O'Haver
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
%
if nargin~=7;smoothtype=1;end  % smoothtype=1 if not specified in argument
if smoothtype>3;smoothtype=3;end
if smoothtype<1;smoothtype=1;end
if smoothwidth<1;smoothwidth=1;end
smoothwidth=round(smoothwidth);
peakgroup=round(peakgroup);
if smoothwidth>1,
    d=fastsmooth(deriv(x,y),smoothwidth,smoothtype);
else
    d=deriv(x,y);
end
n=round(peakgroup/2+1);
P=[0 0 0 0 0];
vectorlength=length(y);
peak=1;
AmpTest=AmpThreshold;
for j=2*round(smoothwidth/2)-1:length(y)-smoothwidth-1,
    if sign(d(j)) > sign (d(j+1)), % Detects zero-crossing
        if d(j)-d(j+1) > SlopeThreshold*y(j), % if slope of derivative is larger than SlopeThreshold
            if or(y(j) > AmpTest, y(j+1) > AmpTest),  % if height of peak is larger than AmpThreshold (new version by Anthony Willey)
          % if y(j) > AmpTest,  % if height of peak is larger than AmpThreshold (old version)
                xx=zeros(size(peakgroup));yy=zeros(size(peakgroup));
                for k=1:peakgroup, % Create sub-group of points near peak
                    groupindex=j+k-n+2;
                    if groupindex<1, groupindex=1;end
                    if groupindex>vectorlength, groupindex=vectorlength;end
                    xx(k)=x(groupindex);yy(k)=y(groupindex);
                end
                if peakgroup>2,
                    [Height, Position, Width]=gaussfit(xx,yy);
                    PeakX=real(Position);   % Compute peak position and height of fitted parabola
                    PeakY=real(Height);
                    MeasuredWidth=real(Width);
                    % if the peak is too narrow for least-squares technique to work
                    % well, just use the max value of y in the sub-group of points near peak.
                else
                    PeakY=max(yy);
                    pindex=val2ind(yy,PeakY);
                    PeakX=xx(pindex(1));
                    MeasuredWidth=0;
                end
                % Construct matrix P. One row for each peak detected,
                % containing the peak number, peak position (x-value) and
                % peak height (y-value). If peak measurement fails and
                % results in NaN, or if the measured peak height is less
                % than AmpThreshold, skip this peak
                if isnan(PeakX) || isnan(PeakY) || PeakY<AmpThreshold,
                    % Skip this peak
                else % Otherwise count this as a valid peak
                    P(peak,:) = [round(peak) PeakX PeakY MeasuredWidth  1.0646.*PeakY*MeasuredWidth];
                    peak=peak+1; % Move on to next peak
                end
            end
        end
    end
end
% ------------------------------------------------------------------------
function [Height, Position, Width]=gaussfit(x,y)
% Converts y-axis to a log scale, fits a parabola
% (quadratic) to the (x,ln(y)) data, then calculates
% the position, width, and height of the
% Gaussian from the three coefficients of the
% quadratic fit.  This is accurate only if the data have
% no baseline offset (that is, trends to zero far off the
% peak) and if there are no zeros or negative values in y.
%
% Example 1: Simplest Gaussian data set
% [Height, Position, Width]=gaussfit([1 2 3],[1 2 1]) 
%    returns Height = 2, Position = 2, Width = 2
%
% Example 2: best fit to synthetic noisy Gaussian
% x=50:150;y=100.*gaussian(x,100,100)+10.*randn(size(x));
% [Height,Position,Width]=gaussfit(x,y) 
%   returns [Height,Position,Width] clustered around 100,100,100.
%
% Example 3: plots data set as points and best-fit Gaussian as line
% x=[1 2 3 4 5];y=[1 2 2.5 2 1];
% [Height,Position,Width]=gaussfit(x,y);
% plot(x,y,'o',linspace(0,8),Height.*gaussian(linspace(0,8),Position,Width))

% Copyright (c) 2012, Thomas C. O'Haver
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
maxy=max(y);
for p=1:length(y),
    if y(p)<(maxy/100),y(p)=maxy/100;end
end % for p=1:length(y),
z=log(y);
coef=polyfit(x,z,2);
a=coef(3);
b=coef(2);
c=coef(1);
Height=exp(a-c*(b/(2*c))^2);
Position=-b/(2*c);
Width=2.35482/(sqrt(2)*sqrt(-c));
% ------------------------------------------------------------------------
function FWHM=halfwidth(x,y)
% function FWHM=halfwidth(x,y) computes the approximate full width at half
% maximum of the central maximum of any shape peak that has a zero
% baseline. Not highly accurate if the function is too sparsely sampled.
% Tom O'Haver (toh@umd.edu) April 2016
%
% Example 1:
% x=-5:.1:5;y=sinc(x);
% plot(x,y);
% FWHM=halfwidth(x,y)
%
% Example 2:
% x=[0:.1:10];
% W=3; % W is the true half-width
% y=gaussian(x,5,W);
% FWHM=halfwidth(x,y);
% plot(x,y);
% hold on;plot([5-FWHM/2 5+FWHM/2],[max(y)/2 max(y)/2],'r');hold off
%
try
    maxy=max(y);
    halfy=maxy/2;
    indmax=round(val2ind(y,maxy));
    % oy=y-maxy/2;
    n=indmax(1);
    while y(n)>halfy,
        n=n-1;
    end
    y1=interp1([y(n) y(n+1)],[x(n) x(n+1)],max(y)/2);
    n=indmax(1);
    while y(n)>halfy,
        n=n+1;
    end
    y2= interp1([y(n-1) y(n)],[x(n-1) x(n)],max(y)/2);
    FWHM=y2-y1;
catch
    FWHM=NaN;
end
% ------------------------------------------------------------------------
function ydc=deconvgauss(x,y,w)
% function ydc=deconvgauss(y,tc) deconvolutes a Gaussian function of
% width 'w' from vector y, returning the deconvoluted result
c=gaussian(x,0,w)+gaussian(x,max(x),w);   % Gaussian convolution function, c
ydc=ifft(fft(y)./fft(c)).*sum(c);
% ------------------------------------------------------------------------
function ydc=deconvexp(x,y,tc)
% function ydc=deconvgauss(x,y,tc) deconvolutes an exponential function of
% time constant 'tc' from vector y, returning the deconvoluted result
c=exp(-x./tc);    % exponential trailing convolution function, c
ydc=ifft(fft(y)./fft(c)).*sum(c);
% ------------------------------------------------------------------------
function s=addup(x)
s=sum(x);
