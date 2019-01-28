function [obj,ind] = ExtractSubDomain(obj,bou,keep_inverse)
p = obj.p; t = obj.t; 
if nargin == 1 || (nargin == 3 && isempty(bou))
    plot(p(:,1),p(:,2),'k.');
    h = impoly;
    bou  = h.getPosition;
end
if nargin == 2
    keep_inverse = 0;
end
if size(bou,1) == 2
     bou = [bou(1,1) bou(2,1);
            bou(1,1) bou(2,2); ...
            bou(1,2) bou(2,2);
            bou(1,2) bou(2,1); ...
            bou(1,1) bou(2,1)];
end
bxy1 = p(t(:,1),:); bxy2 = p(t(:,2),:); bxy3 = p(t(:,3),:); 
in1 = inpoly(bxy1,bou); in2 = inpoly(bxy2,bou); in3 = inpoly(bxy3,bou);
in = in1 & in2 & in3;

if keep_inverse == 0
    t(~in,:) = [];
else
    t(in,:) = [];
end
[p1,t] = fixmesh(p,t);
[~,~,ind] =  intersect(p1,p,'rows','stable');
obj.p = p1; obj.t = t;  
if ~isempty(obj.b)
    obj.b = obj.b(ind); 
end
if ~isempty(obj.bx)
    obj.bx = obj.bx(ind); 
end
if ~isempty(obj.by)
    obj.by = obj.by(ind);
end
end