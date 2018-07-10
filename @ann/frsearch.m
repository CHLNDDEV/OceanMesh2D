function [idx dst inr] = frsearch(ann, query, r, k, epsl, asm)
% 
% kNN FR search
%   Usage:
%     [idx dst inr] = frsearch(ann, query, rad, k, asm, eps)
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
if nargin == 5
    asm = true;
end

if ~asm
    k = k+1;
end

if ~isa(query, ann.ccls)
    query = ann.cfun(query);
end

[idx dst inr] = annmex(ann.modes.FRSEARCH, ann.kd_ptr, query, k, epsl, r);

if ~asm 
    if dst(1) == 0
        dst(1) = [];
        idx(1) = [];
        inr = inr - 1;
    else
        dst(end) = [];
        idx(end) = [];
    end
end
idx = idx + 1; % fix zero indexing of ann