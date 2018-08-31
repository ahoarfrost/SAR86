#!/bin/bash

sJavaDir=/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/scripts

java -cp $sJavaDir/SpeciesDistributionModeling.jar edu.ucsf.Rasters.AcsiiToNetcdf.AsciiToNetcdfLauncher \
--sVariableName=par \
--sSource=MODIS \
--sHistory='smigen par=A20132442013273.L3m_MO_PAR_par_9km.nc.param' \
--sInstitution='NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean Biology Processing Group' \
--sCellMethods='Mean' \
--sReferences='http://oceandata.sci.gsfc.nasa.gov' \
--sUnits=einstein m^-2 day^-1 \
--sRasterDataPath='/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/par/par_rasterdata.csv' \
--sVariable=par \
--sOutputPath='/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/par_monthly_compiled.nc' \
--sTitle='MODIS Level-3 Standard Mapped Image' \
--sLongName='Photosynthetically Available Radiation, R. Frouin'
