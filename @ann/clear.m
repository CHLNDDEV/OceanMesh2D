function clear(anno)
if ~strcmp(class(anno),'ann')
    error('ann:close','invalid handle');
end
if isfield(anno,'working_flag') && anno.working_flag
    anno.kd_ptr = annmex(anno.modes.CLOSE, anno.kd_ptr);
end
close(anno);