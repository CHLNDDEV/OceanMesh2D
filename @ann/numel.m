function varargout = numel(anno)
if ~strcmp(class(anno),'ann')
    error('ann:numel','invalid handle');
end
if nargout >= 1
    varargout{1} = anno.npts;
end
if nargout == 2
    varargout{2} = anno.dim;
end