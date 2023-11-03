#!/bin/bash
####################################################
# This script downloads m_map and some required    #
# datasets for the examples                        #
#                                                  #
# Requires wget and unzip                          #
####################################################

# Can toggle parameters if desired but the script
# checks for their existance too
# set to 'true' to download

m_map=true  # m_map mapping toolbox (reqd)
gshhs=true  # global shoreline
srtm=true   # SRTM15+V2.1 global bathymetry -> {user can select
gebco=false  # GEBCO_2020 global bathymetry  -> {desired source

keepzip=false # keep downloaded zip files

# download m_map archive and unzip
if $m_map; then
    if [ -d "m_map" ]; then
        echo "\"m_map\" directory already exists"
    else
        echo "Setting up m_map"
        if [ ! -f "m_map.zip" ]; then 
            wget --no-check-certificate --no-hsts "https://www.eoas.ubc.ca/~rich/m_map1.4.zip" -O m_map.zip
        fi
        unzip -q m_map.zip
        if ! $keepzip; then rm m_map.zip; fi
    fi
else
    echo "Skipping m_map setup"
fi

# move to datasets directory
if [ ! -d "datasets" ]; then
    mkdir datasets
fi
 
cd datasets/

# download m_map archive and unzip
if $gshhs; then
    if [ -d "GSHHS_shp" ]; then
        echo "\"datasets/GSHHS_shp\" directory (global shoreline) already exists"
    else
        echo "Setting up GSHHS shoreline"
        if [ ! -f "gshhs.zip" ]; then 
            wget --no-check-certificate --no-hsts "http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip" -O gshhs.zip
        fi
        unzip -q gshhs.zip "GSHHS_shp/*"
        if ! $keepzip; then rm gshhs.zip; fi
    fi
else
    echo "Skipping m_map"
fi

# download SRTM15+ bathymetry
if $srtm; then
    if [ -f "SRTM15+.nc" ]; then
        echo "\"datasets/SRTM15+.nc\" global bathymetry file already exists"
    else
        echo "Setting up SRTM15+ bathymetry"
        wget --no-check-certificate --no-hsts "https://topex.ucsd.edu/pub/srtm15_plus/SRTM15_V2.4.nc" -O SRTM15+.nc
    fi
else
    echo "Skipping SRTM15+ bathymetry"
fi

# download GEBCO_2020.nc bathymetry
if $gebco; then
    if [ -f "GEBCO_2020.nc" ]; then
        echo "\"datasets/GEBCO_2020.nc\" global bathymetry file already exists"
    else
        echo "Setting up GEBCO_2020 bathymetry"
        wget --no-check-certificate --no-hsts "https://www.bodc.ac.uk/data/open_download/gebco/gebco_2021/zip/" -O GEBCO_2020.nc
    fi
else
    echo "Skipping GEBCO_2020 bathymetry"
fi

read -p "Press [ENTER] to continue ..." answer
