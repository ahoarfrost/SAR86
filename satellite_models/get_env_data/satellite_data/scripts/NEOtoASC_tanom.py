#convert globalTanomaly and currents to ESRI asc, which aren't in separate monthly files

import netCDF4 as nc
import numpy as np
import os

###globalTanomaly###
#okay new idea - subset .nc file to get np array that is the data you want to write;
#write header then data as asc file, bypass .nc files.

years = ['2008','2009','2010','2011','2012','2013']
month_names = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
year_mo = []
for year in years:
    for mo in month_names:
        y_mo = year+mo
        year_mo.append(y_mo)

gt = nc.Dataset('tanomaly/GlobalTAnomaly_1200ERSSTv4.nc')
metadata_info = dict([("title", ""),("references",""),("source",""),("institution",""),("cell_methods",""),("variable",""),("long_name",""),("units","")])
for time_ix,time in enumerate(range(1537,1609)):
    time_label = year_mo[time_ix]
    #subset tempanomaly for that month
    sub = gt.variables['tempanomaly'][time,:,:]
    sub_filled = sub.filled(fill_value=-9999)

    ncols = sub.shape[1]
    nrows = sub.shape[0]
    lon_center = round(-180+((float(360)/float(ncols))/2),8)
    lat_center = round(-90+((float(180)/float(nrows))/2),8)
    cell_size = round((float(180)/float(nrows)),8)
    fill = -9999

    asc_name = os.getcwd()+'/tanomaly/tanom_'+time_label+'..asc' #make sure this is going to be the same as your make_rasterdata.py output
    #write esri .asc file
    f = open(asc_name,'w')
    f.write("ncols "+str(ncols)+"\n")
    f.write("nrows "+str(nrows)+"\n")
    f.write("xllcenter "+str(lon_center)+"\n")
    f.write("yllcenter "+str(lat_center)+"\n")
    f.write("cellsize "+str(cell_size)+"\n")
    f.write("nodata_value "+str(fill)+"\n")
    f.write("\n".join(" ".join(map(str, x)) for x in sub_filled))
    f.close()

    title = "GISTEMP Surface Temperature Analysis"
    references = "http://data.giss.nasa.gov/gistemp/"
    source = "http://data.giss.nasa.gov/gistemp/"
    institution = "NASA Goddard Institute for Space Studies"
    cell_methods = "mean"
    variable = 'tempanomaly'
    long_name = "Surface temperature anomaly"
    units = "K"
    history = "Created 2016-12-12 11:05:37 by SBBX_to_nc 2.0 - ILAND=1200, IOCEAN=NCDC/ER4, Base: 1951-1980"

    #check if all but history entries of this metadata set are the same as previous - if not, replace with new; if so, keep
    metadata_temp = dict([("title", title),("references",references),("source",source),("institution",institution),("cell_methods",cell_methods),("variable",variable),("long_name",long_name),("units",units)])
    if metadata_temp!=metadata_info:
        metadata_info = metadata_temp

gt.close()


###currents###
