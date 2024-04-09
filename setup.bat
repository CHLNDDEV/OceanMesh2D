@echo off

:: ####################################################
:: # This script downloads m_map and some required    #
:: # datasets for the examples                        #
:: #                                                  #
:: # Uses curl and tar from Windows                   #
:: ####################################################

:: Can toggle parameters if desired but the script
:: checks for their existence too
:: set to 'true' to download

:: m_map mapping toolbox (reqd)
set m_map=true
:: global shoreline
set gshhs=true
:: SRTM15+V2.1 global bathymetry -> {user can select}
set srtm=true
:: GEBCO_2020 global bathymetry  -> {desired source}
set gebco=false

:: keep downloaded zip files
set keepzip=false

:: download m_map archive and unzip
if "%m_map%" == "true" (
    if EXIST "m_map" (
        echo "m_map" directory already exists
    ) else (
        echo Setting up m_map
        if NOT EXIST "m_map.zip" curl "https://www.eoas.ubc.ca/~rich/m_map1.4.zip" -o m_map.zip
        tar -xf m_map.zip 
        if NOT "%keepzip%" == "true" del m_map.zip > nul
    )
) else (
    echo Skipping m_map
)

:: move to datasets directory
if NOT EXIST "datasets" mkdir datasets

pushd datasets

:: download GSHHS_shp shoreline
if "%gshhs%" == "true" (
    if EXIST "GSHHS_shp" (
        echo "datasets/GSHHS_shp" directory ^(global shoreline^) already exists
    ) else (
        echo Setting up GSHHS shoreline
        if NOT EXIST "gshhs.zip" curl "http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip" -o gshhs.zip
        tar -xf gshhs.zip GSHHS_shp
        if NOT "%keepzip%" == "true" del gshhs.zip > nul
    )
) else (
    echo Skipping GSHHS shoreline
)

:: download SRTM15+ bathymetry
if "%srtm%" == "true" (
    if EXIST "SRTM15+.nc" (
        echo "datasets/SRTM15+.nc" global bathymetry file already exists
    ) else (
        echo Setting up SRTM15+ bathymetry
        curl "https://topex.ucsd.edu/pub/srtm15_plus/SRTM15_V2.4.nc" -o SRTM15+.nc
    )
) else (
    echo Skipping SRTM15+ bathymetry
)

:: download GEBCO_2020.nc bathymetry
if "%gebco%" == "true" (
    if EXIST "GEBCO_2020.nc" (
        echo "datasets/GEBCO_2020.nc" global bathymetry file already exists
    ) else (
        echo Setting up GEBCO_2020 bathymetry
        curl "https://www.bodc.ac.uk/data/open_download/gebco/gebco_2021/zip/" -o GEBCO_2020.nc
    )
) else (
    echo Skipping GEBCO_2020 bathymetry
)

popd 

pause
