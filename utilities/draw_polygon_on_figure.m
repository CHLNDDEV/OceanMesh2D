function [polygon] = draw_polygon_on_figure(ptol)
% Draw a polygon on an open figure by clicking the screen. 
% Type the character `h` on the keyboard to finish drawing the polygon. 
% By default, the function merges points less than `ptol` distance apart
% from one another. By default `ptol` is set to 0.1 units. 
%
% Usage:
%   polygon = draw_polygon_figure()
%
%   polygon = draw_polygon_figure(0.2)
%
% Inputs:
%   ptol: distance to merge nearby points (default=0.1). 
% 
% Outputs:
%   polygon: an array of coordinates that the user has drawn.  
%
% Author: Keith J. Roberts
% Created: July 2021

 if nargin < 1
     ptol=0.01;
 end
 h = gcf();
 first = 1;
 while(1)
     [x,y] = ginput(2);
     line(x,y)
     if first
         line_coordinates = [x,y];
         first = 0;
     else
         line_coordinates = [line_coordinates; x,y];
     end
     isKeyPressed = ~isempty(get(h,'CurrentCharacter'));
     if isKeyPressed
         break
     end
 end
 snap = max(max(line_coordinates,[],1)-min(line_coordinates,[],1),[],2)*ptol;
 [~,ix,~] = unique(round(line_coordinates/snap)*snap,'rows','stable');
 polygon = line_coordinates(ix,:);
 hold on; plot(polygon(:,1),polygon(:,2),'r--');
 end
