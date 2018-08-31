#for each NEO .nc file in a folder for a single data type, convert to ESRI .asc file
#generate metadata.txt file w/ metadata need
#generate xxx_RasterData.csv with list asc files and time start/end

import netCDF4 as nc
import numpy as np
import os

#folders to convert files in
data_types = ['sst','chla','par','poc','pic','npp'] #current and globalTanomaly still need reprocess, add later
vars = ['sst','chlor_a','par','poc','pic','npp']

wd = os.getcwd()

#for each folder of data,
for dt_ix,dt in enumerate(data_types):
    print "processing: %s" % dt
    #make list all the .nc files in a directory
    nc_files = []
    for file in os.listdir(dt):
        if file.endswith(".nc")|file.endswith(".hdf"):
            nc_files.append(file)
    print "-%s netcdf files to process" % len(nc_files)
    os.chdir(dt)

    #for each .nc file, generate .asc and metadata .txt
    metadata_info = dict([("title", ""),("references",""),("source",""),("institution",""),("cell_methods",""),("variable",""),("long_name",""),("units","")])
    for ncf in nc_files:
        asc_name = os.getcwd()+'/'+ncf[:-3]+'.asc' #make sure this is going to be the same as your make_rasterdata.py output
        #if asc file doesn't already exist, process
        if os.path.isfile(asc_name):
            print "--asc file exists, moving on; %s" % asc_name
        else:
            print "--processing file: %s" % ncf
            root = nc.Dataset(ncf)
            var = vars[dt_ix]
            variable = root.variables[var][:,:]
            try:
                var_filled = variable.filled(fill_value=-9999)
            except AttributeError:
                #npp isn't masked array so this throws an error, but the fill value is already -9999 so can skip 
                var_filled = variable
            ncols = variable.shape[1]
            nrows = variable.shape[0]
            lon_center = round(-180+((float(360)/float(ncols))/2),8)
            lat_center = round(-90+((float(180)/float(nrows))/2),8)
            cell_size = round((float(180)/float(nrows)),8)
            fill = -9999

            #write esri .asc file
            f = open(asc_name,'w')
            f.write("ncols "+str(ncols)+"\n")
            f.write("nrows "+str(nrows)+"\n")
            f.write("xllcenter "+str(lon_center)+"\n")
            f.write("yllcenter "+str(lat_center)+"\n")
            f.write("cellsize "+str(cell_size)+"\n")
            f.write("nodata_value "+str(fill)+"\n")
            f.write("\n".join(" ".join(map(str, x)) for x in var_filled))
            f.close()

            #define and write metadata .txt file
            if var=='npp':
                title = 'Monthly Net Primary Production'
                references = 'http://www.science.oregonstate.edu/ocean.productivity/index.php'
                source = 'VGPM model'
                institution = 'Ocean Productivity Group, Oregon State University'
                cell_methods = 'mean'
                variable = var
                long_name = 'Net Primary Productivity, VGPM model'
                units = 'mg C / m**2 / day'
                history = ''
            else:
                title = str(root.title)
                references = str(root.publisher_url)
                source = str(root.instrument)
                institution = str(root.institution)
                cell_methods = str(root.measure)
                variable = var
                long_name = str(root.variables[var].long_name)
                units = str(root.variables[var].units)
                history = str(root.history)

            #check if all but history entries of this metadata set are the same as previous - if not, replace with new; if so, keep
            metadata_temp = dict([("title", title),("references",references),("source",source),("institution",institution),("cell_methods",cell_methods),("variable",variable),("long_name",long_name),("units",units)])
            if metadata_temp!=metadata_info:
                metadata_info = metadata_temp

            root.close()
            
            #remove .nc ncf file to keep storage down 
            os.remove(ncf)

            #write metadata file for one data type
            text_name = os.getcwd()+'/metadata_'+dt+'.txt'
            g = open(text_name,'w')
            g.write("--sSource="+metadata_info['source']+" \ \n")
            g.write("--sHistory="+history+" \ \n")
            g.write("--sInstitution="+metadata_info['institution']+" \ \n")
            g.write("--sCellMethods="+metadata_info['cell_methods']+" \ \n")
            g.write("--sReferences="+metadata_info['references']+" \ \n")
            g.write("--Units="+metadata_info['units']+" \ \n")
            g.write("--sVariable="+metadata_info['variable']+" \ \n")
            g.write("--sTitle="+metadata_info['title']+" \ \n")
            g.write("--sLongName="+metadata_info['long_name']+" \ \n")
            g.close()

    os.chdir(wd)
