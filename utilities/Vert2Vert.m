function [adj,xadj] = Vert2Vert( t )
%   VERT2VERT Complicated output for Vert2Vert
    bars = [t(:,[1,2]); t(:,[1,3]); t(:,[2,3])];        
    bars = [bars; fliplr(bars)]; 
    bars = sortrows(bars,1); 
    xadj = find(diff(bars(:,1))==1); 
    adj  = bars(:,2); 
    xadj = [1;xadj];
    xadj = [xadj; numel(adj)];
end

