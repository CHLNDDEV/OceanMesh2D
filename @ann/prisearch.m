function [idx dst] = prisearch(ann, query, k, epsl, asm)
% 
% kNN search
%   Usage:
%     [idx dst] = prisearch(ann, query, k, asm, eps)
%
% Inputs:
%   ann - ann class object
%   query - (d)x(N) query points
%   k - number of nearest nieghbors (if points in ann < k than less than k
%                                    points are returned)
%   epsl - epsilon search precision 
%   asm - allow self match flag, if false points with dst = 0 are ignored 
%         (defualt is true)
%
if nargin == 4
    asm = true;
end

if ~asm
    k = k+1;
end

if ~isa(query, ann.ccls)
    query = ann.cfun(query);
end

[idx dst] = annmex(ann.modes.PRISEARCH, ann.kd_ptr, query, k, epsl);

if ~asm
    gsm = dst(1,:)==0;
    dst(1:end-1,gsm) = dst(2:end,gsm);
    idx(1:end-1,gsm) = idx(2:end,gsm);
    dst(end,:) = [];
    idx(end,:) = [];
end
idx = idx + 1; % fix zero indexing of ann