#NEOtoASC.sh

#download nasa data
#needs get_nasa_data.sh, and get_nasa_data.py
sh get_nasa_data.sh

#process nasa data to esri asc, metadata txt files
#make sure have all approp nc files in folders 'sst','chla','par','poc','pic','npp'
python NEOtoASC.py

#make xxx_rasterdata.csv files for each data type
python make_rasterdata.py

#make xxx_AsciitoNetcdf.sh by hand using metaata file
#requires SpeciesDistributionModeling.jar, xxx_rasterdata.csv, path to asc files
sh sst_AsciitoNetcdf.sh
sh chla_AsciitoNetcdf.sh
sh par_AsciitoNetcdf.sh
sh tanom_AsciitoNetcdf.sh
sh poc_AsciitoNetcdf.sh
sh pic_AsciitoNetcdf.sh
sh npp_AsciitoNetcdf.sh
#sh currents_AsciitoNetcdf.sh
