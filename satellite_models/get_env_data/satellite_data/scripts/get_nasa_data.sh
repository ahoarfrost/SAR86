#download global t anomaly
#this has time bnds, monthly data from 1880-present. Need pull out months care about.
#mkdir globalTanomaly
#curl -o globalTanomaly/GlobalTAnomaly_1200ERSSTv4.nc.gz http://data.giss.nasa.gov/pub/gistemp/gistemp1200_ERSSTv4.nc.gz
#gunzip globalTanomaly/GlobalTAnomaly_1200ERSSTv4.nc.gz

mkdir sst
mkdir chla
mkdir par
mkdir poc
mkdir pic
mkdir npp
mkdir current

python get_nasa_data.py
