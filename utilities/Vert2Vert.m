function [adj,xadj,nnv,vtov] = Vert2Vert( t )
%   VERT2VERT Complicated output for Vert2Vert
    bars = [t(:,[1,2]); t(:,[1,3]); t(:,[2,3])];        
    bars = [bars; fliplr(bars)]; 
    bars = sortrows(bars,1); 
    xadj = find(diff(bars(:,1))==1); 
    adj  = bars(:,2); 
    xadj = [1;xadj];
    xadj = [xadj; numel(adj)];
    if nargout > 2
        if nargout > 3
            vtov = zeros(max(diff(xadj))-1,length(xadj)-1);
        end
        nnv = zeros(length(xadj)-1,1);
        for i = 1:length(xadj)-1
            nei = unique(adj(xadj(i):xadj(i+1)-1));
            nnv(i) = length(nei);
            if nargout > 3
                vtov(1:nnv(i),i) = nei;
            end
        end
        if nargout > 3 
            vtov(max(nnv)+1:end,:) = [];
        end
    end
end

