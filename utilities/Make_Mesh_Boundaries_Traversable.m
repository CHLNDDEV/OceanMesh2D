function obj = Make_Mesh_Boundaries_Traversable( obj, dj_cutoff, nscreen )
%  obj =  Make_Mesh_Boundaries_Traversable(obj,dj_cutoff,nscreen)
%  A msh object (containing p and t) is  "cleaned" and returned. 
%  ncscreen ~= 0 will display info to screen
%  See "Exterior Check" description below for definition of dj_cutoff
%
%  Alternates between checking interior and exterior portions 
%  of the graph exhaustively until convergence is obtained, defined as:
%  Having no nodes connected to more than 2 boundary edges. 
%
%  Interior Check: Deletes elements that are within the interior of the 
%  mesh so that no nodes are connected to more than 2 boundary edges. For
%  example, a spit could become very thin in a middle portion so that you
%  a node is connected to two elements but four boundary edges, in a
%  bow-tie type formation. This code will delete one of those connecting
%  elements to ensure the spit is continous and only two boundary edges
%  are connected to that node. In the case of a choice between elements to 
%  delete, the one with the lowest quality is chosen. 
%
%  Exterior Check: Finds small disjoint portions of the graph and removes 
%  them using a breadth-first search. The individual disjoint portions are 
%  removed based on dj_cutoff. 
%  dj_cutoff <= 1 indicates the proportional of the total area of the mesh
%  that, below which, a disjoint portion is discarded.
%  dj_cutoff > 1 indicates the absolute area in km2 that, below which, a 
%  disjoint portion is discarded.
%  dj_cutoff can be used to keep only the largest area of the mesh 
%  (e.g. dj_cutoff = 0.5) or a few largest areas (e.g. dj_cutoff = 0.1), 
%  or it can be set to ensure that lakes or polders of certain sizes are 
%  kept (e.g. dj_cutoff = 1 [km2])
%  
%  Copyright (C) 2018  Keith Roberts & William Pringle
%  Algorithm written by William Pringle, CHL, UND 2017
%  Improvements by Keith Roberts, CHL, UND 2017
%  WJP: Organisation of code Jan 2018
%  KJR: Speed-up, CHL, UND 2018. 
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.  

% Entering the code
disp('Making mesh boundaries traversable...');  
% This is to avoid heaps of warnings in the triangulation call which are
% unneccesary
warning('off','all')

% Get p and t out of obj
p = obj.p; t = obj.t;

% Delete disjoint nodes
[p,t] = fixmesh(p,t);

% Get boundary edges and nodes
[etbv,vxe] = extdom_edges2( t, p ) ;

% WJP sometimes we may wanna delete some exterior portions even with a valid
% mesh so allow entry even in this case.
if numel(etbv) == numel(vxe)
    etbv(end+1,:) = 1; 
end
% Loop until all the nodes only have two boundary edges
%(the number of boundary edges will equal number of boundary nodes)

while numel(etbv) > numel(vxe)
    
    % Delete elements in the exterior of the mesh
    t = delete_exterior_elements(p,t,dj_cutoff,nscreen);
    
    % Delete disjoint nodes
    [p,t] = fixmesh(p,t);

    % Delete elements in the interior of the mesh
    t = delete_interior_elements(p,t,nscreen);
     
    % Delete disjoint nodes
    [p,t] = fixmesh(p,t);
    
    % Get boundary edges and nodes
    [etbv,vxe] = extdom_edges2( t, p ) ;
    
    %if numel(vxe) > numel(etbv)
    %   error(['number of boundary vertices larger than boundary edges', ...
    %          ', try a larger dj_cutoff']) 
    %end
end
% Finished cleaning
disp('ALERT: finished cleaning up mesh..'); 
% Turn warnings back on
warning('on','all')
% Put back into the msh object
obj.p = p; obj.t = t;
end
% The sub-functions...
%% Delete elements outside the main mesh depending on dj_cutoff input
% dj_cutoff >= 1
%    area in km2
% dj_cutoff < 1
%    proportion of the total mesh area
function t = delete_exterior_elements(p,t,dj_cutoff,nscreen)
L = size(t,1); 
t1 = t; t = [];
global MAP_PROJECTION
if isempty(MAP_PROJECTION)
    % need to project
    m_proj('Azimuthal Equal-area','lon',[min(p(:,1)),max(p(:,1))],...
           'lat',[min(p(:,2)),max(p(:,2))]) ;
    % Do the transformation
    [X,Y] = m_ll2xy(p(:,1),p(:,2)); 
else
    % already projected
    X = p(:,1); Y = p(:,2);  
end
A = sum(polyarea(X(t1(:,1:3))',Y(t1(:,1:3))'));  An = A;
if dj_cutoff >= 1
    Re2 = (6378.137)^2; An = Re2*An;
    % Absolute area
    while An > dj_cutoff

        % Peform the Breadth-First-Search to get nflag
        nflag = BFS(p,t1);

        % Get new triangulation and its area
        t2 = t1(nflag == 1,:);
        An = Re2*sum(polyarea(X(t2(:,1:3))',Y(t2(:,1:3))')); 
        
        % If large enough at t2 to the triangulation
        if An > dj_cutoff
            t = [t; t2]; 
        end
        % Delete where nflag == 1 since this patch didn't meet the fraction
        % limit criterion.
        t1(nflag == 1,:) = []; 
        % Calculate the remaining area       
        An = Re2*sum(polyarea(X(t1(:,1:3))',Y(t1(:,1:3))')); 
        
    end
elseif dj_cutoff > 0
    % Area proportion
    while An/A > dj_cutoff

        % Peform the Breadth-First-Search to get nflag
        nflag = BFS(p,t1);

        % Get new triangulation and its area
        t2 = t1(nflag == 1,:);
        An = sum(polyarea(X(t2(:,1:3))',Y(t2(:,1:3))')); 
        
        % If large enough at t2 to the triangulation
        if An/A > dj_cutoff
            t = [t; t2]; 
        end
        % Delete where nflag == 1 since this patch didn't meet the fraction
        % limit criterion.
        t1(nflag == 1,:) = []; 
        % Calculate the remaining area       
        An = sum(polyarea(X(t1(:,1:3))',Y(t1(:,1:3))')); 
        
    end
elseif dj_cutoff < 0
    error('Keep cannot be negative 0')
else
    % keep is zero, do nothing
end
if nscreen
    disp(['  ACCEPTED: deleting ' num2str(L-size(t,1)) ...
          ' elements outside main mesh']) ;
end
if size(t,1) < 1
    error(['All elements have been deleted... something wrong? ' ...
           'dj_cutoff is set to' num2str(dj_cutoff)])
end
    
end

%% Delete interior elements
function t = delete_interior_elements(p,t,nscreen)
% Get updated boundary edges.
etbv = extdom_edges2( t, p ) ;
% Get all nodes that are on edges.
[nodes_on_edge,~,n] = unique(etbv(:));
% Count how many edges a node appears in.
I     = accumarray(n,1:numel(n),[],@(x){x});
count = cellfun('length',I);
%
[vtoe,nne] = VertToEle(t);
% Get the nodes which appear more than twice and delete element connected
% to these nodes where all nodes of element are on boundary edges
del_elem_idx = [];
for i = nodes_on_edge(count > 2)'
    con_elem = vtoe(1:nne(i),i);
    n = 0; del_elem = [];
    for elem = con_elem'
        I = etbv(:) == t(elem,1); 
        J = etbv(:) == t(elem,2); 
        K = etbv(:) == t(elem,3); 
        % All nodes on element are boundary edges
        if any(I) && any(J) && any(K) 
            n = n + 1;
            del_elem(n) = elem;
        end
    end
    if n == 1
        % Only one element to delete.
        del_elem_idx(end+1) = del_elem;
    elseif n > 1
        % Delete worst quality qualifying element.
        tq = gettrimeshquan( p, t(del_elem,:));
        [~,idx] = min(tq.qm);
        del_elem_idx(end+1) = del_elem(idx);
    else
        % No connected elements have all nodes on boundary edge so we
        % select the worst quality connecting element.
        tq = gettrimeshquan( p, t(con_elem,:));
        [~,idx] = min(tq.qm);
        del_elem_idx(end+1) = con_elem(idx);
    end
end

if nscreen
    disp(['  ACCEPTED: deleting ' num2str(length(del_elem_idx)) ...
          ' elements inside main mesh'])
end
t(del_elem_idx,:) = [];

end

function nflag = BFS(p,t1)

    % Select a random element.
    EToS =  randi(size(t1,1),1);

    % Get element-to-element connectivity.
    tri = triangulation(t1,p);
    nei = tri.neighbors;

    % Traverse grid deleting elements outside.
    ic = zeros(ceil(sqrt(size(t1,1))*2),1);
    ic0 = zeros(ceil(sqrt(size(t1,1))*2),1);
    nflag = zeros(size(t1,1),1);
    ic(1) = EToS;
    icc  = 1;

    % Using BFS loop over until convergence is reached (i.e., we
    % obtained a connected region).
    while icc
        ic0(1:icc) = ic(1:icc);
        icc0 = icc;
        icc = 0;
        for nn = 1:icc0
            i = ic0(nn);
            % Flag the current element as OK
            nflag(i) = 1;
            % Search neighbouring elements
            nb = nei(i,:); 
            nb(isnan(nb)) = []; 
            % Flag connected neighbors as OK
            for nnb = 1:length(nb)
                if ~nflag(nb(nnb))
                    icc = icc + 1;
                    ic(icc) = nb(nnb);
                    nflag(nb(nnb)) = 1;
                end
            end
        end
    end
end
