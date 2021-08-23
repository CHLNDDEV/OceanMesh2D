addpath(genpath('Tests'));
addpath(genpath('../datasets/'));

if ~exist('../m_map/') || ~exist('GSHHS_f_L1.shp')

   if ~exist('../m_map/')
      error('OceanMesh2D requires m_map package! Run setup.sh')
   end
   if ~exist('GSHHS_f_L1.shp')
      error('Need GSHHS global shoreline to run tests. Run setup.sh')
   end

else

   TestSanity 

   TestEleSizes
   
   TestInterp

   if exist('SRTM15+.nc')
   
      TestECGC 

   else
   
      warning('Need to download SRTM15+.nc to run TestECGC. Run setup.sh') 

   end

   if exist('PostSandyNCEI.nc') && exist('PostSandyNCEI.shp')

      TestJBAY

   else

      warning('Need to download PostSandyNCEI data to run TestJBAY. Available from the google drive referenced on the README page') 

   end

end

exit
