#!/bin/bash

sJavaDir=/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/scripts/

java -cp $sJavaDir/SpeciesDistributionModeling.jar edu.ucsf.Rasters.AcsiiToNetcdf.AsciiToNetcdfLauncher \
--sVariableName=tempanomaly \
--sSource='http://data.giss.nasa.gov/gistemp/'' \
--sHistory='Created 2016-12-12 11:05:37 by SBBX_to_nc 2.0 - ILAND=1200, IOCEAN=NCDC/ER4, Base: 1951-1980' \
--sInstitution='NASA Goddard Institute for Space Studies' \
--sCellMethods='mean' \
--sReferences='http://data.giss.nasa.gov/gistemp/' \
--sUnits=K \
--sRasterDataPath='/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/tanomaly/tanom_rasterdata.csv' \
--sVariable=tempanomaly \
--sOutputPath='/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/tanomaly_compiled.nc' \
--sTitle='GISTEMP Surface Temperature Analysis' \
--sLongName='Surface temperature anomaly'
