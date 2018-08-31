#!/bin/bash

sJavaDir=/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/scripts/

java -cp $sJavaDir/SpeciesDistributionModeling.jar edu.ucsf.Rasters.AcsiiToNetcdf.AsciiToNetcdfLauncher \
--sVariableName=poc \
--sSource=MODIS \
--sHistory='smigen par=A20112442011273.L3m_MO_POC_poc_9km.nc.param' \
--sInstitution='NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean Biology Processing Group' \
--sCellMethods='Mean' \
--sReferences='http://oceandata.sci.gsfc.nasa.gov' \
--sUnits=mg m^-3 \
--sRasterDataPath='/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/poc/poc_rasterdata.csv' \
--sVariable=poc \
--sOutputPath='/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/poc_monthly_compiled.nc' \
--sTitle='MODIS Level-3 Standard Mapped Image' \
--sLongName='Particulate Organic Carbon, D. Stramski, 2007 (443/555 version)'
