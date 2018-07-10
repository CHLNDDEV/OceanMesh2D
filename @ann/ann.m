function anno = ann(data)
%
% ANN class 
% 
% Usage:
%   anno = ann(pts)
%
% Input:
%   pts - (d)x(N) matrix of d-dimensional vectors representing N points
%
% Output:
%   anno - object handle to be used in class' methods
%
% Methods:
%   ksearch:
%       k-nn search
%       [idx dst] = ksearch(anno, pt, k, eps, [asm])
%       
%       idx - indices of matching points:  (k)x(#points)
%       dst - square distance of the points
%
%       anno - active ann object handle
%       pt   - query point(s): a (d)x(#points) matrix
%       k    - required k aprox. nearest neighbors
%       eps  - accuracy multiplicative factor
%       asm  - allow self match flag (optional, default: true)
%
%   prisearch:
%       priority search
%       [idx dst] = ksearch(anno, pt, k, eps, [asm])
%
%       see ksearch for inputs/outputs explaination
%
%   frsearch:
%       fixed radius search
%       [idx dst inr] = frsearch(anno, pt, r, k, eps, [asm])
%
%       idx - indices of matching points
%       dst - square distance of the points
%       inr - how many points are in r (range query)
%
%       anno - active ann object handle
%       pt   - query point (only one this time, (d)x1 vector)
%       r    - search radius (_not_ squared)
%       k    - required k aprox. nearest neighbors. Will return at most k nieghbors 
%              of distance at most r.
%       eps  - accuracy multiplicative factor
%       asm  - allow self match flag (optional, default: true)
%       
%       note that numel(idx) and numel(dst) = min(k, inr)
%       inr might be larger than k reflecting how many points are actually
%       in the search radius.
%
%   close:
%      you MUST explicitly close the ann class handle when you are done with it.
%      This is crucial for releasing memory and other resources.
%
%      anno = close(anno);
%
%
%   This wrapper for Matlab was written by Shai Bagon (shai.bagon@weizmann.ac.il).
%   Department of Computer Science and Applied Mathmatics
%   Wiezmann Institute of Science
%   http://www.wisdom.weizmann.ac.il/~bagon
%
%	The core cpp application was written by David M. Mount and Sunil Arya
%	(available from http://www.cs.umd.edu/~mount/ANN/):
%
%	ANN is a library for approximate nearest neighbor searching,
%	based on the use of standard and priority search in kd-trees
%	and balanced box-decomposition (bbd) trees. Here are some
%	references to the main algorithmic techniques used here:
%
%		kd-trees:
%			Friedman, Bentley, and Finkel, ``An algorithm for finding
%				best matches in logarithmic expected time,'' ACM
%				Transactions on Mathematical Software, 3(3):209-226, 1977.
%
%		Priority search in kd-trees:
%			Arya and Mount, ``Algorithms for fast vector quantization,''
%				Proc. of DCC '93: Data Compression Conference, eds. J. A.
%				Storer and M. Cohn, IEEE Press, 1993, 381-390.
%
%		Approximate nearest neighbor search and bbd-trees:
%			Arya, Mount, Netanyahu, Silverman, and Wu, ``An optimal
%				algorithm for approximate nearest neighbor searching,''
%				5th Ann. ACM-SIAM Symposium on Discrete Algorithms,
%				1994, 573-582.
% 
%   This software is provided under the provisions of the 
%   Lesser GNU Public License (LGPL).
%   This software can be used only for research purposes, you should cite 
%   the aforementioned papers in any resulting publication.
%
%
%   The Software is provided "as is", without warranty of any kind.
%
%

if ~isnumeric(data)
    error('ann:open','data must be numeric 2D matrix');
end

if ndims(data) > 2 || numel(data) == 0
    error('ann:open','Data must be 2D');
end

[anno.ccls anno.cfun] = get_data_class();
anno.modes = modes();

[anno.dim anno.npts] = size(data);

% convert data to be of supported class
if ~isa(data, anno.ccls)
    data = anno.cfun(data);
end

anno.kd_ptr = annmex(anno.modes.OPEN, data);

anno.working_flag = 1;
anno = class(anno, 'ann');
