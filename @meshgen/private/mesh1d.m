
function [pout,t,converged]=mesh1d(poly,fh0,h,fix,boubox,box_num0,varargin)
% Mesh Generator using Distance Functions.
%   [P,T]=mesh1d(FDIST,FH,H,BOX,FIX,FDISTPARAMS)
%
%   POUT:           Resampled Node positions (N X 2)
%      T:           Edge indices (NTx2))
%   POLY:           The polygon (NX X 2
%      FH:          Edge length function
%      H:           Smallest edge length
%      FIX:         Indices of Poly Fixed node positions (NFIX X 1)

% Preprocessing steps - move into parametric space
%poly = add_flairs(poly, 50); 

X = poly(:,1);
Y = poly(:,2);

% find sharp corners & add flairs 
% [sharp_corners] = find_sharp_corners(poly, 40);
% if ~isempty(sharp_corners)
%     fix = [fix; sharp_corners]; 
% end


% add ending point to emulate periodic mesh generation
%X(end+1) = X(1,:) + eps;
%Y(end+1) = Y(1,:) + eps;

xd = diff(X);
yd = diff(Y);
dist = sqrt(xd.^2 + yd.^2);
u = cumsum(dist);
u = [0; u];
totaldist = u(end);

np = totaldist / min(h);
% kjr check if total np will be too little 
% i.e., less than 5 
if np < 5
    pout = [];
    t = [];
    converged=1; 
    return 
end 
t = linspace(0,max(u),ceil(np));

xn = interp1(u,X,t);
yn = interp1(u,Y,t);

u_fix = [];
if ~isempty(fix) 
    % Convert fixed points to parametric space
    % determine the distance these are points are from the zero point. 
    u_fix = [];
    for i = 1:length(fix)
        xd_tmp = diff(X(1:fix(i)));
        yd_tmp = diff(Y(1:fix(i)));
        dist_tmp = sqrt(xd_tmp.^2 + yd_tmp.^2);
        if isempty(dist_tmp)
            continue
        elseif cumsum(dist_tmp) == totaldist
            continue
        else
            csum = cumsum(dist_tmp);
        end
        u_fix = [u_fix; csum(end)];
    end
end

A = 0;
B = totaldist;

box  = [A,B]';

fdist = @(x) my_1d_sdf(x, A, B);

tmph = nestedHFx([xn',yn'],boubox,fh0,h);

    function [h_inner] = nestedHFx(x,boubox,fh0,h)
        h_inner= x(:,1)*0;
        for box_num = 1:length(h)  % For each bbox, find the points that are in and calculate
            fh_l = fh0{box_num};
            if box_num > 1
                iboubox = boubox{box_num}(1:end-1,:) ;
                inside = inpoly(x,iboubox) ;
            else
                inside = true(size(h_inner));
            end
            h_inner(inside) = feval(fh_l,x(inside,:))./111e3; 
        end
    end

F = griddedInterpolant(t,tmph,'linear','linear');

fh = @(x) F(x);

dim=size(box,2);
ptol=.01;
L0mult=1+.4/2^(dim-1);
deltat=.10; geps=1e-3*h(box_num0);
deps=sqrt(eps)*h(box_num0);

% 1. Create initial distribution in bounding box
if dim==1
    p=(box(1):h(box_num0):box(2))';
else
    disp('Only 1D is supported')
    %error
end

% 2. Remove points outside the region, apply the rejection method
p=p(feval(fdist,p)<geps,:);
r0=feval(fh,p);
% Add the fixed points.
if ~isempty(u_fix)
    p=[u_fix; p(rand(size(p,1),1)<min(r0)^dim./r0.^dim,:)];
else 
    p = p(rand(size(p,1),1)<min(r0)^dim./r0.^dim,:);
end
% find their indicies
p = unique(p,'stable'); 
[p,I] = sort(p);
[~,fixed_point_ind] = sort(I);
fixed_point_ind = fixed_point_ind(2:length(u_fix)-1);
% find the position of the fixed points  
N=size(p,1);

count=0;
while 1
    % 3. Retriangulation
    t=[(1:length(p)-1)',(2:length(p))'];
    pmid=zeros(size(t,1),dim);
    for ii=1:dim+1
        pmid=pmid+p(t(:,ii),:)/(dim+1);
    end
    t=t(feval(fdist,pmid)<-geps,:);
    % 4. Describe each edge by a unique pair of nodes
    pair=zeros(0,2);
    localpairs=nchoosek(1:dim+1,2);
    for ii=1:size(localpairs,1)
        pair=[pair;t(:,localpairs(ii,:))];
    end
    pair=unique(sort(pair,2),'rows');
    % 5. Graphical output of the current mesh
    %fprintf('Retriangulation #%d\n',count)
    count = count + 1;
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
    %dp(1:size(fix,1),:)=0;
    dp(fixed_point_ind,:)=0;
    % lock the first and last point 
    dp(1) = 0; dp(end) = 0;
    
    p=p+deltat*dp;
    
    % 7. Bring outside points back to the boundary
    d=feval(fdist,p); 
    ix=d>0;
    %ix(1:length(fix),:) = 0;
    ix(fixed_point_ind,:) = 0;
    gradd=zeros(sum(ix),dim);
    for ii=1:dim
        a=zeros(1,dim);
        a(ii)=deps;
        d1x=feval(fdist,p(ix,:)+ones(sum(ix),1)*a);
        gradd(:,ii)=(d1x-d(ix))/deps;
    end
    p(ix,:)=p(ix,:)-d(ix)*ones(1,dim).*gradd;

    % 8. Termination criterion
    maxdp=max(deltat*sqrt(sum(dp(d<-geps,:).^2,2)));

    %disp(maxdp)
    if maxdp<ptol*h
        %disp('Converged...')
        converged=1;
        % put back in the real space
        pout(:,1) = interp1(u,X,p);
        pout(:,2) = interp1(u,Y,p);
        %figure; plot(pout(:,1),pout(:,2),'rs')
        if any(isnan(pout(:,1)))
            warning('Nans detected in 1d mesh')
        end
        break
    end
    if count > 5000
        warning('Mesh1d: some line segments did not converge...')
        converged=0;
        pout(:,1) = interp1(u,X,p);
        pout(:,2) = interp1(u,Y,p);
        if any(isnan(pout(:,1)))
            warning('Nans detected in 1d mesh')
        end
        break
    end

end
end
