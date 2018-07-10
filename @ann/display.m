function display(anno)
fprintf(1,'\n %s = \n', inputname(1));
fprintf(1,'\tdim = %d\n', anno.dim);
fprintf(1,'\t#points = %d\n', anno.npts);
fprintf(1,'\tclass = %s\n', anno.ccls);
fprintf(1,'\tkd_ptr = %lu\n', anno.kd_ptr);
