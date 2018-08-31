mkdir sst
mkdir chla
mkdir par
mkdir poc
mkdir pic
mkdir npp

#download .nc files from repositories 
python get_nasa_data.py

#for each NEO .nc file in a folder for a single data type, convert to ESRI .asc file
#generate metadata.txt file w/ metadata need
python NEOtoASC.py
#generate xxx_RasterData.csv with list asc files and time start/end

#make xxx_rasterdata.csv files for each data type
python make_rasterdata.py


