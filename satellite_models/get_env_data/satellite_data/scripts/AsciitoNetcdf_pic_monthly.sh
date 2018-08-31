#!/bin/bash

sJavaDir=/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/scripts

java -cp $sJavaDir/SpeciesDistributionModeling.jar edu.ucsf.Rasters.AcsiiToNetcdf.AsciiToNetcdfLauncher \
--sVariableName=pic \
--sSource=MODIS \
--sHistory='smigen par=A20132442013273.L3m_MO_PIC_pic_9km.nc.param' \
--sInstitution='NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean Biology Processing Group' \
--sCellMethods='Mean' \
--sReferences='http://oceandata.sci.gsfc.nasa.gov' \
--sUnits=mol m^-3 \
--sRasterDataPath='/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/pic/pic_rasterdata.csv' \
--sVariable=pic \
--sOutputPath='/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/pic_monthly_compiled.nc' \
--sTitle='MODIS Level-3 Standard Mapped Image' \
--sLongName='Calcite Concentration, Balch and Gordon'
