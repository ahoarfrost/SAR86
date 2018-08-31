#!/bin/bash

sJavaDir=/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/scripts/

java -cp $sJavaDir/SpeciesDistributionModeling.jar edu.ucsf.Rasters.AcsiiToNetcdf.AsciiToNetcdfLauncher \
--sVariableName=sst \
--sSource=MODIS \
--sHistory='smigen par=A20141212014151.L3m_MO_SST_sst_9km.nc.param' \
--sInstitution='NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean Biology Processing Group' \
--sCellMethods='Mean' \
--sReferences='http://oceandata.sci.gsfc.nasa.gov' \
--sUnits=degree_C \
--sRasterDataPath='/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/sst/sst_rasterdata.csv' \
--sVariable=sst \
--sOutputPath='/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/sst_compiled.nc' \
--sTitle='MODIS Level-3 Standard Mapped Image' \
--sLongName='Sea Surface Temperature'
