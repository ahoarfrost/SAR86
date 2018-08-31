#!/bin/bash

sJavaDir=/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/scripts

java -cp $sJavaDir/SpeciesDistributionModeling.jar edu.ucsf.Rasters.AcsiiToNetcdf.AsciiToNetcdfLauncher \
--sVariableName=npp \
--sSource='MODIS, SeaWifs' \
--sHistory='not specified' \
--sInstitution='Oregon State' \
--sCellMethods='not specified' \
--sReferences='Ocean Productivity Group' \
--sUnits=mgC m-2 day-1 \
--sRasterDataPath='/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/npp/npp_rasterdata.csv' \
--sVariable=npp \
--sOutputPath='/Users/Adrienne/Projects/SAR86/satellite_models/get_env_data/satellite_data/npp_monthly_compiled.nc' \
--sTitle='Net Primary Productivity (VGPM model)' \
--sLongName='Net Primary Productivity'
