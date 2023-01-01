#!/bin/bash
####################################################
# This script downloads m_map and some required    #
# datasets for the examples                        #
#                                                  #
# Requires wget and zip                            #
####################################################
# Can toggle parameters if desired but the script
# checks for their existance too
m_map=true  # m_map mapping toolbox (reqd)
gshhs=true  # global shoreline
srtm=true   # SRTM15+V2.1 global bathymetry -> {user can select
gebco=false # GEBCO_2020 global bathymetry  -> {desired source

if $m_map; then
  if [ -d "m_map" ]; then
     echo "m_map directory already exists"
  else
     # download m_map archive and unzip
     wget "http://www.eos.ubc.ca/%7Erich/m_map1.4.zip" -O m_map.zip
     unzip m_map.zip
     rm m_map.zip
  fi
fi

# move to datasets directory
cd datasets/

if $gshhs; then
  if [ -d "GSHHS_shp" ]; then
     echo "GSHHS_shp directory (global shoreline) already exists"
  else
     # download m_map archive and unzip
     wget "http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip" -O gshhs.zip
     unzip gshhs.zip "GSHHS_shp/*"
     rm gshhs.zip
  fi
fi

if $srtm; then
  if [ -f "SRTM15+.nc" ]; then
     echo "SRTM15+.nc global bathymetry file already exists"
  else
     # download SRTM15+ bathymetry
     wget "https://topex.ucsd.edu/pub/srtm15_plus/SRTM15_V2.4.nc"" -O SRTM15+.nc
  fi
fi

if $gebco; then
  if [ -f "GEBCO_2020.nc" ]; then
     echo "GEBCO_2020.nc global bathymetry file already exists"
  else
     # download GEBCO_2020.nc bathymetry
     wget "https://www.bodc.ac.uk/data/open_download/gebco/gebco_2021/zip/" -O GEBCO_2020.nc
  fi
fi
