function ind = nearest_neighbor_map(m_old, m_new)
% Determine the indices of the nearest neighbors from m_old.p to m_new.p 
% to do for example a trivial transfer of nodal attributes or bathymetry. 
%
% Inputs
% m_old: the old mesh object 
% m_new: the new mesh object
%
% Output
% The indices of points from m_old to m_new 
%
% Example
% Transfer bathymetry data from one grid to a new one.
% m_new.b = m_old.b(ind)
ind = ourKNNsearch(m_old.p',m_new.p',1);
end