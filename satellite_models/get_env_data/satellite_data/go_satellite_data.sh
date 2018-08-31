mkdir sst
mkdir chla
mkdir par
mkdir poc
mkdir pic
mkdir npp
#mkdir tanomaly 

#download .nc files from repositories 
python scripts/get_nasa_data.py
#npp is in .gz format, unzip
gunzip npp/*.gz

#for each NEO .nc file in a folder for a single data type, convert to ESRI .asc file
#generate metadata.txt file w/ metadata need
#outputs metadata.txt file and .asc file for each .nc
python scripts/NEOtoASC.py
#python scripts/NEOtoASC_tanom.py

#make xxx_rasterdata.csv files for each data type with list asc files and time start/end
python scripts/make_rasterdatacsv.py

#make xxx_AsciitoNetcdf.sh by hand using metaata file
#requires SpeciesDistributionModeling.jar, xxx_rasterdata.csv, path to asc files
#outputs compiled .nc file with all the months together and unified metadata/data format
sh scripts/AsciitoNetcdf_sst_monthly.sh
sh scripts/AsciitoNetcdf_chla_monthly.sh
sh scripts/AsciitoNetcdf_par_monthly.sh
#sh scripts/AsciitoNetcdf_tanom.sh
sh scripts/AsciitoNetcdf_poc_monthly.sh
sh scripts/AsciitoNetcdf_pic_monthly.sh
sh scripts/AsciitoNetcdf_npp_monthly.sh


