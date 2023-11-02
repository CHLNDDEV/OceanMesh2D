@echo off

:: ####################################################
:: # This script downloads m_map and some required    #
:: # datasets for the examples                        #
:: #                                                  #
:: # Requires wget and unzip                          #
:: ####################################################
:: # Can toggle parameters if desired but the script
:: # checks for their existance too

:: set to 'true' to download
:: m_map mapping toolbox (reqd)
set m_map=true
:: global shoreline
set gshhs=true
:: SRTM15+V2.1 global bathymetry -> {user can select}
set srtm=true
:: GEBCO_2020 global bathymetry  -> {desired source}
set gebco=false

:: download m_map archive and unzip
if "%m_map%" == "true" (
  if EXIST "m_map" (
    echo m_map directory already exists
  ) else (
    wget --no-check-certificate "https://www.eoas.ubc.ca/~rich/m_map1.4.zip" -O m_map.zip
    powershell -command "Expand-Archive -Force 'm_map.zip' -DestinationPath ."
    del m_map.zip > nul
  )
)

:: move to datasets directory
if NOT EXIST "datasets" mkdir datasets
pushd datasets

:: download GSHHS_shp shoreline
if "%gshhs%" == "true" (
  if EXIST "GSHHS_shp" (
    echo GSHHS_shp directory ^(global shoreline^) already exists
  ) else (
    wget --no-check-certificate "http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip" -O gshhs.zip
    powershell -command "Expand-Archive -Force 'gshhs.zip' -DestinationPath ."
    del gshhs.zip > nul
  )
)

:: download SRTM15+ bathymetry
if "%srtm%" == "true" (
  if EXIST "SRTM15+.nc" (
    echo SRTM15+.nc global bathymetry file already exists
  ) else (
    wget --no-check-certificate "https://topex.ucsd.edu/pub/srtm15_plus/SRTM15_V2.4.nc" -O SRTM15+.nc
  )
)

:: download GEBCO_2020.nc bathymetry
if "%gebco%" == "true" (
  if EXIST "GEBCO_2020.nc" (
    echo "GEBCO_2020.nc global bathymetry file already exists"
  ) else (
    wget --no-check-certificate "https://www.bodc.ac.uk/data/open_download/gebco/gebco_2021/zip/" -O GEBCO_2020.nc
  )
)

popd 

pause
