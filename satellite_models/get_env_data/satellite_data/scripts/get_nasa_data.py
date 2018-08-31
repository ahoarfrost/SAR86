import shutil
import urllib2
from contextlib import closing
import gzip
from ftplib import FTP
import os

#download the sst data
#ftp://podaac-ftp.jpl.nasa.gov/allData/modis/L3/aqua/11um/v2014.0/9km/monthly/[YYYY]/A[YYYY][sDDD][YYYY][eDDD].L3m_MO_SST_sst_9km.nc

years = ['2008','2009','2010','2011','2012','2013']
leap_id = [True,False,False,False,True,False]
leap_days = [['001','031'],['032','060'],['061','091'],['092','121'],['122','152'],['153','182'],['183','213'],['214','244'],['245','274'],['275','305'],['306','335'],['336','366']]
non_leap_days = [['001','031'],['032','059'],['060','090'],['091','120'],['121','151'],['152','181'],['182','212'],['213','243'],['244','273'],['274','304'],['305','334'],['335','365']]
month_names = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']

#download sst data, comes in month packets
print "gathering sst..."
list_files = []
for year_ix,year in enumerate(years):
    print "processing: %s" % year
    if leap_id[year_ix]==True:
        days = leap_days
    else:
        days = non_leap_days
    #download .nc file for every month in year
    dir = 'allData/modis/L3/aqua/11um/v2014.0/9km/monthly/'+year+'/'
    ftp = FTP('podaac-ftp.jpl.nasa.gov')
    ftp.login()
    ftp.cwd(dir)
    for month_ix,month in enumerate(month_names):
        #check if file exists
        filein = 'A'+year+days[month_ix][0]+year+days[month_ix][1]+'.L3m_MO_SST_sst_9km.nc'
        fileout = 'sst/sst_'+year+month_names[month_ix]+'.L3m_MO_SST_sst_9km.nc'
        if os.path.isfile(fileout):
            print "--%s exists, moving on" % month
        else:
            print "--%s" % month
            with open(fileout, "wb") as f:
                ftp.retrbinary('RETR ' + filein, f.write)
            f.close()
            list_files.append(fileout)
    ftp.quit()


#download chla data - similar to sst ftp
#ftp://podaac-ftp.jpl.nasa.gov/allData/modis/L3/aqua/chlA/v2014.0/4km/monthly/
#[YYYY]/A[YYYY][sDDD][YYYY][eDDD].L3m_MO_CHL_chlor_a_4km.nc
print "gathering chla..."
list_files = []
for year_ix,year in enumerate(years):
    print "processing: %s" % year
    if leap_id[year_ix]==True:
        days = leap_days
    else:
        days = non_leap_days
    #download .nc file for every month in year
    dir = 'allData/modis/L3/aqua/chlA/v2014.0/4km/monthly/'+year+'/'
    ftp = FTP('podaac-ftp.jpl.nasa.gov',timeout=60)
    ftp.login()
    ftp.cwd(dir)
    for month_ix,month in enumerate(month_names):
        filein = 'A'+year+days[month_ix][0]+year+days[month_ix][1]+'.L3m_MO_CHL_chlor_a_4km.nc'
        fileout = 'chla/chla_'+year+month_names[month_ix]+'.L3m_MO_CHL_chlor_a_4km.nc'
        if os.path.isfile(fileout):
            print "--%s exists, moving on" % month
        else:
            print "--%s" % month
            with open(fileout, "wb") as f:
                ftp.retrbinary('RETR ' + filein, f.write)
            f.close()
            list_files.append(fileout)
    ftp.quit()

#download PAR data - on ocean color data, similar to sst but http
#http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A[YYYY][sDDD][YYYY][eDDD].L3m_MO_PAR_par_9km.nc
print "gathering par..."
list_files = []
for year_ix,year in enumerate(years):
    print "processing: %s" % year
    if leap_id[year_ix]==True:
        days = leap_days
    else:
        days = non_leap_days
    #download .nc file for every month in year
    for month_ix,month in enumerate(month_names):
        url = 'https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A'+year+days[month_ix][0]+year+days[month_ix][1]+'.L3m_MO_PAR_par_9km.nc'
        fileout = 'par/par_'+year+month_names[month_ix]+'.L3m_MO_PAR_par_9km.nc'
        if os.path.isfile(fileout):
            print "--%s exists, moving on" % month
        else:
            print "--%s" % month
            with closing(urllib2.urlopen(url)) as r:
                with open(fileout, 'wb') as f:
                    shutil.copyfileobj(r, f)
            list_files.append(fileout)


#POC data - ocean data, similar to par
#http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A20083062008335.L3m_MO_POC_poc_9km.nc
print "gathering poc..."
list_files = []
for year_ix,year in enumerate(years):
    print "processing: %s" % year
    if leap_id[year_ix]==True:
        days = leap_days
    else:
        days = non_leap_days
    #download .nc file for every month in year
    for month_ix,month in enumerate(month_names):
        url = 'https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A'+year+days[month_ix][0]+year+days[month_ix][1]+'.L3m_MO_POC_poc_9km.nc'
        fileout = 'poc/poc_'+year+month_names[month_ix]+'.L3m_MO_POC_poc_9km.nc'
        if os.path.isfile(fileout):
            print "--%s exists, moving on" % month
        else:
            print "--%s" % month
            with closing(urllib2.urlopen(url)) as r:
                with open(fileout, 'wb') as f:
                    shutil.copyfileobj(r, f)
            list_files.append(fileout)


#PIC data - ocean data, similar to par
#http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A20030012003031.L3m_MO_PIC_pic_9km.nc
print "gathering pic..."
list_files = []
for year_ix,year in enumerate(years):
    print "processing: %s" % year
    if leap_id[year_ix]==True:
        days = leap_days
    else:
        days = non_leap_days
    #download .nc file for every month in year
    for month_ix,month in enumerate(month_names):
        url = 'https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A'+year+days[month_ix][0]+year+days[month_ix][1]+'.L3m_MO_PIC_pic_9km.nc'
        fileout = 'pic/pic_'+year+month_names[month_ix]+'.L3m_MO_PIC_pic_9km.nc'
        if os.path.isfile(fileout):
            print "--%s exists, moving on" % month
        else:
            print "--%s" % month
            with closing(urllib2.urlopen(url)) as r:
                with open(fileout, 'wb') as f:
                    shutil.copyfileobj(r, f)
            list_files.append(fileout)

#NPP data - ocean productivity Group
#http://orca.science.oregonstate.edu/data/2x4/monthly/vgpm.r2014.m.chl.m.sst/hdf/vgpm.2008032.hdf.gz
#vgpm.[YYYY][sDDD].hdf.gz
#vgpm = net primary production (units of mg C / m**2 / day) based on the standard vgpm algorithm
print "gathering npp..."
list_files = []
for year_ix,year in enumerate(years):
    print "processing: %s" % year
    if leap_id[year_ix]==True:
        days = leap_days
    else:
        days = non_leap_days
    #download .nc file for every month in year
    for month_ix,month in enumerate(month_names):
        url = 'http://orca.science.oregonstate.edu/data/2x4/monthly/vgpm.r2014.m.chl.m.sst/hdf/vgpm.'+year+days[month_ix][0]+'.hdf.gz'
        fileout = 'npp/npp_'+year+month_names[month_ix]+'.hdf.gz'
        if os.path.isfile(fileout):
            print "--%s exists, moving on" % month
        else:
            print "--%s" % month
            with closing(urllib2.urlopen(url)) as r:
                with open(fileout, 'wb') as f:
                    shutil.copyfileobj(r, f)
            list_files.append(fileout)
            gzip.open(fileout,'rb')

#check if netcdf python module can handle hdf

#OSCAR current speed/direction data - podaac
#ftp://podaac-ftp.jpl.nasa.gov/allData/oscar/preview/L4/oscar_third_deg/oscar_vel[YYYY].nc.gz
#by year, each file has time arrays; 5d res only; 1/3 deg res; need separate out all the 5d slots (?)
'''
print "gathering current speed & direction..."
list_files = []
for year_ix,year in enumerate(years):
    print "processing: %s" % year
    #download .nc file for every month in year
    dir = 'allData/oscar/preview/L4/oscar_third_deg/'
    ftp = FTP('podaac-ftp.jpl.nasa.gov')
    ftp.login()
    ftp.cwd(dir)
    filein = 'oscar_vel'+year+'.nc.gz'
    fileout = 'current/current_'+year+'.nc.gz'
    with open(fileout, "wb") as f:
        ftp.retrbinary('RETR ' + filein, f.write)
    f.close()
    list_files.append(fileout)
    ftp.quit()
    gzip.open(fileout,'rb')
'''
