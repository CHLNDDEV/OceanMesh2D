function [p] = draw_polygon_on_figure(ptol)
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
p = line_coordinates(ix,:);
hold on; plot(p(:,1),p(:,2),'r--');
end
