% what type of architecture we are running on? 64 or 32 bits?
cd @ann/private
try
    if ~isempty(findstr(computer,'64'))
        
        mex -O -DUSE64BITS  annmex.cpp ANN.cpp   ...
            kd_pr_search.cpp kd_split.cpp kd_util.cpp  ...
            kd_fix_rad_search.cpp kd_search.cpp kd_tree.cpp my_mex_mem.cpp ...
            kd_dump.cpp bd_pr_search.cpp bd_tree.cpp perf.cpp bd_fix_rad_search.cpp bd_search.cpp brute.cpp
        
    else
        mex -g  annmex.cpp ANN.cpp bd_pr_search.cpp bd_tree.cpp kd_dump.cpp ...
            kd_pr_search.cpp kd_split.cpp kd_util.cpp perf.cpp ...
            bd_fix_rad_search.cpp bd_search.cpp brute.cpp ...
            kd_fix_rad_search.cpp kd_search.cpp kd_tree.cpp my_mex_mem.cpp
    end
catch
    lasterr
end

cd ../../