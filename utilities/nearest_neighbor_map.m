function ind = nearest_neighbor_map(m_old, m_new, type)
% ind = nearest_neighbor_map(m_old, m_new, type)
%
% Determine the indices of the nearest neighbors from m_old.p to m_new.p 
% to do for example a trivial transfer of nodal attributes or bathymetry. 
%
% Inputs
% m_old: the old mesh object 
% m_new: the new mesh object
% type (optional): 'approx' to use ourKNNsearch ANN wrapper
%                  'precise' to use built in MATLAB knnsearch
%
% Output
% The indices of points from m_old to m_new 
%
% Example
% Transfer bathymetry data from one grid to a new one.
% m_new.b = m_old.b(ind)

% checking type input
if nargin < 3 || isempty(type)
    type = 'approx';
end
% check that knnsearch built-in exists if precise selection, 
% otherwise use approx
if strcmp(type,'precise') && ~exist('knnsearch') 
    type = 'approx';
end

if strcmp(type,'approx')
   ind = ourKNNsearch(m_old.p',m_new.p',1);
elseif strcmp(type,'precise')
   ind = knnsearch(m_old.p,m_new.p); 
else
   error(['Unknown selection for type: ' type])
end
    
end