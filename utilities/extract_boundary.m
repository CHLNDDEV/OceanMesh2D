function [poly,poly_idx,opendat,boudat] = extract_boundary(v_start,v_end,bnde,pts,order,opendat,boudat)
% DESCRIPTION: Given a set of boundary edges and a starting and ending index
%              of a singly- or multi-polygonal region, organize them in a
%              winding order and/or add them to an existing opendat/boudat
%              structure. The program breaks up the boundary so that it
%              will never exceed 100,000 nodes to avoid memory allocation
%              issues adcprep.
%
% INPUTS:
%      v_start: the starting index of the boundary you want to trace
%        v_end: the ending index of the boundary you want to trace.
%         bnde: the indices of each boundary edge as a nbnde x 2 matrix
%         pts:  the x,y locations of all the points in the region
%               stored as an np x 2 matrix.
%         order:the order in which the traversal takes place
%               counter-clockwise (0) or clockwise (1).
%         opendat: open boundary information from a pre-existing grid
%          boudat: land boundary information from a pre-exist
% OUTPUTS:
%          poly: the boundary of each enclosing polygon sorted in winding-order
%                poly is returned as a cell-array of length number of polys.
%      poly_idx: indices of the polygon coordinates in the same format as
%                poly
%       opendat: open boundary information with appened open bou
%       boudat: land boundary information with appended land bou
%
% kjr,UND,CHL,2017
%
%                                           TRAVERSAL METHOD
% Pick any unvisited edge segment [v_start,v_next] and add these vertices to the polygon loop.
% Find the unvisited edge segment [v_i,v_j] that has either v_i = v_next or v_j = v_next and add the other vertex (the one not equal to v_next) to the polygon loop.
% Reset v_next as this newly added vertex, mark the edge as visited and continue from 2.
% Traversal is done when we get back to v_start.
% NOTE: that the signed area will be positive if the vertices are
% oriented counterclockwise, and will be negative if it is oriented clockwise
bnde= unique(bnde,'rows');
active = true(size(bnde,1),1);
p = 0;

[rt,dmy] = find(v_start==bnde);
if isempty(rt), disp('v_start does not exist on boundary, check numbering'); return; end
r  = rt(order+1); % change this only here to from 1 to 2 or 2 to 1 to go left or right
tsel = bnde(r,:);
sel  = tsel(tsel~=v_start);
v_next = sel;
active(r) = 0;
% cut up boundary so land boundary never exceeds 100k nodes.
cut = true;
while cut
    p = p + 1;
    if(p > 1 )
        [rt,dmy] = find(v_next==bnde & active);
        r  = rt(1);
        tsel = bnde(r,:);
        sel  = tsel(tsel~=v_next);
        v_next = sel;
        active(r) = 0;
    end
    
    temp = [];
    temp2= [];
    
    temp  = pts(bnde(r,:)',:);
    temp2 = bnde(r,:)';
    k = 2;
    while v_next~=v_end % terminates when we reach v_end
        rt= (v_next==bnde(:,1) | v_next==bnde(:,2)) &  active;
        r = find(rt,1);
        tsel = bnde(r,:);
        sel=tsel(tsel~=v_next);
        k = k + 1;
        temp(k,:)= pts(sel,:);
        temp2(k,:)= sel;
        active(r) = 0;
        v_next = sel;
        % exceeded max xize, break
        if(k > 100e3),disp('exceed'); break, end
        % reached ending vertex, break
        if(v_next==v_end), cut=false; disp('reached ending'); end
        % exhausted all edges and couldn't connect
        if(~any(active)), cut=false; disp('coudln''t conntect'); break, end
    end
    poly{p}     = temp;
    poly_idx{p} = temp2;
    [area]=parea(poly{p}(:,1),poly{p}(:,2));
    if(order==0) % ccw
        if sign(area)<0
            poly{p} = flipud(poly{p});
            poly_idx{p} = flipud(poly_idx{p});
        end
    else % cw
        if sign(area)>0
            poly{p} = flipud(poly{p});
            poly_idx{p} = flipud(poly_idx{p});
        end
    end
end
for ii = 1 : p
    hold on; plot(poly{ii}(:,1),poly{ii}(:,2),'r-','linewi',2);
end

type = input('What kind of boundary is this, 1 (flux) or 2 (elevation)?');

% if populated
if ~isempty(boudat)
    if(type==1)
        nbou = boudat.nbou;
        nvel = boudat.nvel;
        nvell= boudat.nvell;
        nbvv = boudat.nbvv;
        ibtype = boudat.ibtype;
        type2 = input('What kind of flux boundary is it, 20(island),22(River)?');
        for ii = 1 : length(poly)
            nbou = nbou + 1;
            nvell(nbou) = length(poly{ii}(:,1));
            nvel = nvel + nvell(nbou);
            nbvv(1:nvell(nbou),nbou) = poly_idx{ii}(:);
            ibtype(nbou) = type2;
        end
        boudat.nbou = nbou ;
        boudat.nvel = nvel ;
        boudat.nvell = nvell ;
        boudat.ibtype = ibtype ;
        boudat.nbvv = nbvv ;
    end
end

% if populated 
if ~isempty(opendat)
    if type==2
        nope = opendat.nope;
        nvdll= opendat.nvdll;
        neta = opendat.neta;
        ibtype=opendat.ibtype;
        nbdv  = opendat.nbdv;
        
        for ii = 1 : length(poly)
            nope = nope + 1;
            nvdll(nope) = length(poly_idx{ii}(:,1));
            neta = neta + nvdll(nope);
            ibtype(nope) = 0;
            nbdv(1:nvdll(nope),nope) = poly_idx{ii}(:,1);
        end
        % ocean boundary
        opendat.nope = nope ;
        opendat.neta = neta ;
        opendat.nvdll = nvdll ;
        opendat.ibtype = ibtype ;
        opendat.nbdv = nbdv ;
        
    end
end

% if empty 
if isempty(opendat)
    if(type==2)
        % nothing populated
        nope = 0 ; 
        nvdll = []  ; 
        neta = 0 ; 
        nbdv = [] ; 
        for ii = 1 : length(poly)
            nope = nope + 1;
            nvdll(nope) = length(poly_idx{ii}(:,1));
            neta = neta + nvdll(nope);
            ibtype(nope) = 0;
            nbdv(1:nvdll(nope),nope) = poly_idx{ii}(:,1);
        end
        % ocean boundary
        opendat.nope = nope ;
        opendat.neta = neta ;
        opendat.nvdll = nvdll ;
        opendat.ibtype = ibtype ;
        opendat.nbdv = nbdv ;
    end
end

% if empty 
if isempty(boudat)
    if(type==1)
        type2 = input('What kind of flux boundary is it, 20(island),2(River)?');
        for ii = 1 : length(poly)
            nbou = nbou + 1;
            nvell(nbou) = length(poly{ii}(:,1));
            nvel = nvel + nvell(nbou);
            nbvv(1:nvell(nbou),nbou) = poly_idx{ii}(:);
            ibtype(nbou) = type2;
        end
        boudat.nbou = nbou ;
        boudat.nvel = nvel ;
        boudat.nvell = nvell ;
        boudat.ibtype = ibtype ;
        boudat.nbvv = nbvv ;
    end
end

end
% helper function, computes area of polygon
function [area]=parea(x,y)
n    = length(x);
xp   = [x; x(1)];
yp   = [y; y(1)];
area = 0;
for i = 1:n
    area = area + det([xp(i), xp(i+1); yp(i), yp(i+1)]);
end
area = 1/2*area;
end