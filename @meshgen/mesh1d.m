function [pout,t]=mesh1d(poly,fh,h,fix,varargin)
%OCEANMESH1D Mesh Generator using Distance Functions.
%   [P,T]=DISTMESHND(FDIST,FH,H,BOX,FIX,FDISTPARAMS)
%
%   POUT:           Resampled Node positions (N X 2)
%      T:           Edge indices (NTx2))
%   POLY:           The polygon (NX X 2
%      FH:          Edge length function
%      H:           Smallest edge length
%      FIX:         Indices of Poly Fixed node positions (NFIX X 1)

% Preprocessing steps - move into parametric space
X = poly(:,1);
Y = poly(:,2);
xd = diff(X);
yd = diff(Y);
dist = sqrt(xd.^2 + yd.^2);
u = cumsum(dist);
u = [0; u];
totaldist = u(end);

np = totaldist / h;
t = linspace(0,max(u),ceil(np));

xn = interp1(u,X,t);
yn = interp1(u,Y,t);

A = 0;
B = totaldist;

box  = [A,B]';

fdist = @(x) my_1d_sdf(x, A, B);

tmph = fh([xn',yn']);

F = griddedInterpolant(t,tmph,'linear','linear');

fh = @(x) F(x);

dim=size(box,2);
ptol=.001;
ttol=.1;
L0mult=1+.4/2^(dim-1);
deltat=.1; geps=1e-1*h;
deps=sqrt(eps)*h;

% 1. Create initial distribution in bounding box
if dim==1
    p=(box(1):h:box(2))';
else
    disp('Only 1D is supported')
    %error
end

% 2. Remove points outside the region, apply the rejection method
p=p(feval(fdist,p,varargin{:})<geps,:);
r0=feval(fh,p);
p=[fix; p(rand(size(p,1),1)<min(r0)^dim./r0.^dim,:)];
N=size(p,1);

count=0;
p0=inf;
while 1
    % 3. Retriangulation by Delaunay
    if max(sqrt(sum((p-p0).^2,2)))>ttol*h
        p0=p;
        %if dim==1
        %    t = [(1:length(p)-1)',(2:length(p))'];
        %else
        t=delaunayn(p);
        pmid=zeros(size(t,1),dim);
        for ii=1:dim+1
            pmid=pmid+p(t(:,ii),:)/(dim+1);
        end
        t=t(feval(fdist,pmid,varargin{:})<-geps,:);
        % 4. Describe each edge by a unique pair of nodes
        pair=zeros(0,2);
        localpairs=nchoosek(1:dim+1,2);
        for ii=1:size(localpairs,1)
            pair=[pair;t(:,localpairs(ii,:))];
        end
        pair=unique(sort(pair,2),'rows');
        % 5. Graphical output of the current mesh
        fprintf('Retriangulation #%d\n',count)
        count=count+1;
    end
    
    % 6. Move mesh points based on edge lengths L and forces F
    bars=p(pair(:,1),:)-p(pair(:,2),:);
    L=sqrt(sum(bars.^2,2));
    L0=feval(fh,(p(pair(:,1),:)+p(pair(:,2),:))/2);
    L0=L0*L0mult*(sum(L.^dim)/sum(L0.^dim))^(1/dim);
    F=max(L0-L,0);
    Fbar=[bars,-bars].*repmat(F./L,1,2*dim);
    dp=full(sparse(pair(:,[ones(1,dim),2*ones(1,dim)]), ...
        ones(size(pair,1),1)*[1:dim,1:dim], ...
        Fbar,N,dim));
    dp(1:size(fix,1),:)=0;
    p=p+deltat*dp;
    
    % 7. Bring outside points back to the boundary
    d=feval(fdist,p,varargin{:}); ix=d>0;
    gradd=zeros(sum(ix),dim);
    for ii=1:dim
        a=zeros(1,dim);
        a(ii)=deps;
        d1x=feval(fdist,p(ix,:)+ones(sum(ix),1)*a,varargin{:});
        gradd(:,ii)=(d1x-d(ix))/deps;
    end
    p(ix,:)=p(ix,:)-d(ix)*ones(1,dim).*gradd;
    
    % 8. Termination criterion
    maxdp=max(deltat*sqrt(sum(dp(d<-geps,:).^2,2)));
    if maxdp<ptol*h
        % put back in the real space
        pout(:,1) = interp1(u,X,p);
        pout(:,2) = interp1(u,Y,p);
        break
    end
    if count > 250
        % put back in the real space
        pout(:,1) = interp1(u,X,p);
        pout(:,2) = interp1(u,Y,p);
        break
    end
end
