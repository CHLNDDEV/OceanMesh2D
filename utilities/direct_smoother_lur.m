function [p,t] = direct_smoother_lur(p,t,pfix,nscreen)
% a direct solve for smoothing triangulations...edited by kjr,und 2017
if nscreen
    disp('ALERT: beginning implicit smoothing of mesh..');
end

% Set numerical smoothing parameters
mu = 1;
kinf = 10^12;

nnodes = max(t(:));

% Construct element stiffness matrix
D = 2*eye(2);
Q = [ -1 -mu; mu -1 ];
T = [ -1 -sqrt(3); sqrt(3) -1 ];
O = sparse(2,2);
ke{2} = { D, Q,  O,  Q'; Q', D,  Q,  O; O,  Q', D,  Q;  Q,  O,  Q', D }; %--not used
ke{1} = { D, T, T'; T' D T; T T' D};


I = 1:2:2*nnodes;
% Global assembly
nT = length(t);
xgs = zeros(4*nT*9,1); ygs = zeros(4*nT*9,1);
jjj=0;
% unrolled inner loops for speed
for n = 1:nT
    nodelist = t(n,:);
    i = I(nodelist(:))';
    j = I(nodelist(:))';
    
    jjj=jjj+1;
    ttt1=repmat([i(1),i(1),i(1)+1,i(1)+1],3,1);
    ttt2=repmat([i(2),i(2),i(2)+1,i(2)+1],3,1);
    ttt3=repmat([i(3),i(3),i(3)+1,i(3)+1],3,1);
    ttt=[ttt1,ttt2,ttt3];
    xgs(jjj:jjj+35,1)=reshape(ttt',[],1);
        
    zzz1=[j,j+1,j,j+1];
    zzz=[zzz1,zzz1,zzz1];
    ygs(jjj:jjj+35,1)=reshape(zzz',[],1); 
    jjj = jjj+35;
end
% exploit that lds is repeated (length of chunk is 36).
jjj=0;
for r = 1 : 3
    for s = 1 : 3
        jjj=jjj+1;
        ld = ke{1}{r,s};
        ld = ld(:);
        lds(jjj:jjj+3,1)=ld;
        jjj=jjj+3;
    end
end
rp  = length(xgs)/36;
lds = repmat(lds,[rp,1]);
K   = sparse(xgs,ygs,lds);

% Enforce BCs
[bnde,~]=extdom_edges2(t,p);
iFixed = unique(bnde(:));
if ~isempty(pfix)
   ifix   = ourKNNsearch(p',pfix',1); 
   iFixed = unique([ifix; iFixed]); 
end
%F = sparse(2*nnodes,1);
F = zeros(2*nnodes,1); 
% iFixed is small so assembly approach is trivial.
for r = 1:length(iFixed)
    i = I(iFixed(r));
    F(i:i+1,1) = F(i:i+1,1) + kinf*p(iFixed(r),1:2)';
    K(i,i) = kinf;
    K(i+1,i+1) = kinf;
end
F = sparse(F); 

% Solve for new nodal positions and save to mesh points
c = K\F;
p(:,1:2) = reshape(c,[2,nnodes])';
if nscreen
    disp('ALERT: finished implicit smoothing of mesh..');
end
end






