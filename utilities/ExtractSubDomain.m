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
bcx  = (p(t(:,1),1) + p(t(:,2),1) + p(t(:,3),1))/3;
bcy  = (p(t(:,1),2) + p(t(:,2),2) + p(t(:,3),2))/3;
bcxy = [bcx,bcy];
[in] = inpoly(bcxy,bou);

if keep_inverse == 0
    t(~in,:) = [];
else
    t(in,:) = [];
end
[p1,t] = fixmesh(p,t);
[~,~,ind] =  intersect(p1,p,'rows');
obj.p = p1; obj.t = t; 
end