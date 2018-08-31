#!/bin/bash

sJavaDir=/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/scripts

java -cp $sJavaDir/SpeciesDistributionModeling.jar edu.ucsf.Rasters.AcsiiToNetcdf.AsciiToNetcdfLauncher \
--sVariableName=chlor_a \
--sSource=MODIS \
--sHistory='smigen par=A20141212014151.L3m_MO_CHL_chlor_a_4km.nc.param' \
--sInstitution='NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean Biology Processing Group' \
--sCellMethods='Mean' \
--sReferences='http://oceandata.sci.gsfc.nasa.gov' \
--sUnits=mg m^-3 \
--sRasterDataPath='/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/chla/chla_rasterdata.csv' \
--sVariable=chlor_a \
--sOutputPath='/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/chla_monthly_compiled.nc' \
--sTitle='MODIS Level-3 Standard Mapped Image' \
--sLongName='Chlorophyll Concentration, OCI Algorithm'
