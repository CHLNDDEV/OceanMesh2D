run('../setup_oceanmesh2d.m')

if ~(exist('m_proj','file')==2)
    
    error('OceanMesh2D requires m_map package! Run setup.sh/setup.bat')
    
elseif ~(exist('GSHHS_f_L1.shp','file')==2)

    error('Need GSHHS global shoreline to run tests. Run setup.sh/setup.bat')

else

   TestSanity 

   TestEleSizes
   
   TestInterp

   if exist('SRTM15+.nc','file')==2
   
      TestECGC 

   else
   
      warning('Need to download SRTM15+.nc to run TestECGC. Run setup.sh/setup.bat') 

   end

   if exist('PostSandyNCEI.nc','file')==2 && exist('PostSandyNCEI.shp','file')==2

      TestJBAY

   else

      warning('Need to download PostSandyNCEI data to run TestJBAY. Available from the google drive referenced on the README page')

   end

end

exit
