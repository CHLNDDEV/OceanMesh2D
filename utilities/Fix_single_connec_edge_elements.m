function obj = Fix_single_connec_edge_elements(obj,maxit,nscreen)
% obj = Fix_single_connec_edge_elements(obj,maxit,nscreen)
% Elements in a msh object that have only one connecting edge are 
% deleted (if maxit = [] or inf, then deleted exhaustively, 
% otherwise will delete a max number of maxit times)
% The new msh obj is returned. ncscreen ~= 0 will display info to screen
%
%  Copyright (C) 2018  Keith Roberts & William Pringle
%  Algorithm invented and written by William Pringle, CHL, UND 2017
%  Improvements by Keith Roberts, CHL, UND 2017
%  WJP: Help and organisation of code made clearer Jan 2018
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
if nargin < 2 || isempty(maxit)
    maxit = inf;
end
if nargin < 3
    nscreen = 1;
end

if maxit == 0; return; end 

% Get t and p from mesh object
t = obj.t; p = obj.p;

% Now start the loop
del = 1;
old = size(t,1);
it = 0;
while it < maxit && ~isempty(del)
    tri  = triangulation(t,p) ; 
    nei  = tri.neighbors; 
    conn = ~isnan(nei);
    nnei = sum(conn,2);
    del  = find(nnei==1);
    t(del,:) = [];
    % delete disjoint nodes
    [p,t] = fixmesh(p,t);
    it = it + 1;
end
new = size(t,1);

if nscreen
    disp(['  ACCEPTED: deleted ' num2str(old-new) ' bad' ...
          ' elements that are connected to a single neighboring element '])
end
% Put back into msh obj
obj.p = p; obj.t = t;
%EOF
end

function nn = get_small_connectivity(p,t)
    % Get node connectivity (look for 4)
    [vtoe, enum] = VertToEle(t);
    % Make sure they are not boundary nodes
    bdbars = extdom_edges2(t, p);
    bdnodes = unique(bdbars(:));
    I = find(enum <= 4);
    nn = setdiff(I',bdnodes); % and don't destroy fix_p or bnde!
    % make sure they are not connected to more than one boundary node
    nflag = [];
    for node_check = 1:length(nn)
        ineigh = vtoe(1:enum(nn(node_check)),nn(node_check)) ;
        val = unique(t(ineigh,:)) ;
        bdconect = intersect(val,bdnodes);
        if length(bdconect) > 1
            % connected to more than one bdnode
            nflag(end+1) = node_check;
        end
    end
    nn(nflag) = [];
    return;
end