'''
# # Compiling env data for satellite models
#
# ## Two types of satellite data: historical, and contemporary
#
# I ran the commands in go_satellite_data.sh to get sst_, chla_, poc_, pic_, npp_, and par_ compiled.nc files, which contain the 9km (or 4km for chla) resolution monthly average data for 2009-2013.
#
# I also got a bunch of historical averages for a lot of data (from JLadau), contained in the historical_data folder; a summary of what these fields are is found in historical_data_guide.xlsx. Some of these have depth values, some are just from surface; some are monthly averages, some are annual; some are 50-yr averages, some are decadal; etc.
#
# ### Need to go through each data source, and pull the data point from the correct latitude, longitude, and depth (nearest possible point, won't be exact) for each TARA site
'''

# ## Create TARA dataframe
#
# --latitude
# --longitude
# --year sampled (from sampling date)
# --month sampled (from sampling date)
# --sampling depth
# --filter_range (from filter lower and filter upper)
#

# ## Get general info from samples (lat, lon, depth, etc)

import pandas as pd
import numpy as np
import datetime as dt
import netCDF4

samples = pd.read_csv("TARASampleDescriptionTable.csv")
samples = samples[['Sample label [TARA_station#_environmental-feature_size-fraction]', 'INSDC run accession number(s)',
         'Date/Time [yyyy-mm-ddThh:mm]', 'Latitude [degrees North]', 'Longitude [degrees East]', 'Sampling depth [m]',
         'Size fraction lower threshold [micrometre]', 'Size fraction upper threshold [micrometre]']]
samples.head()

env = pd.DataFrame({'TARA_sample_label':samples['Sample label [TARA_station#_environmental-feature_size-fraction]']})
env['run_ids'] = samples['INSDC run accession number(s)']
env['filter_range'] = samples['Size fraction lower threshold [micrometre]']+'-'+samples['Size fraction upper threshold [micrometre]'].astype('string')
env['latitude'] = samples['Latitude [degrees North]']
env['longitude'] = samples['Longitude [degrees East]']
env['depth_sampled'] = samples['Sampling depth [m]']
env['year_sampled'] = samples['Date/Time [yyyy-mm-ddThh:mm]'].apply(lambda x: dt.datetime.strptime(x, '%Y-%m-%dT%H:%M').year)
env['month_sampled'] = samples['Date/Time [yyyy-mm-ddThh:mm]'].apply(lambda x: dt.datetime.strptime(x, '%Y-%m-%dT%H:%M').month)


# ## Add contemporary satellite data
#
# * sst
# * chla
# * npp
# * par
# * pic
# * poc

#function to get gregorian days since e.g. 1970-01-01
def get_yrmo_gregorian_timedelta(year, month, datetime_since):
    leap_days_in_month = {1:31, 2:29, 3:31, 4:30, 5:31, 6:30, 7:31, 8:31, 9:30, 10:31, 11:30, 12:31}
    common_days_in_month = {1:31, 2:28, 3:31, 4:30, 5:31, 6:30, 7:31, 8:31, 9:30, 10:31, 11:30, 12:31}
    leap_years = {2008:True, 2009:False, 2010:False, 2011:False, 2012:True, 2013:False}

    if leap_years[year]==True:
        day = leap_days_in_month[month]
    else:
        day = common_days_in_month[month]

    return (dt.datetime(year, month, day) - datetime_since).days

#function to pull satellite value from lat-lon array for correct time

def pull_satellite_data(nc,year,month,z_variable,lat_target,lon_target,dep_target, lat="lat",lon="lon",time="time", vert="vert", index_order=['time', 'lat', 'lon']):
    import numpy as np
    indices = {}

    #get lat and lon indexes in netcdf
    lats = nc.variables[lat][:]
    lons = nc.variables[lon][:]
    indices['lat'] = np.argmin(np.abs(lats-lat_target))
    indices['lon'] = np.argmin(np.abs(lons-lon_target))

    #get time index in netcdf; time is usually in gregorian days since 1970-01-01
    if time in nc.variables and len(nc.variables[time][:])>1:
        if nc.variables[time].units=='days since 1970-01-01':
            if len(nc.variables[time])==12:
                indices['time'] = month-1
            else:
                days_since = get_yrmo_gregorian_timedelta(year,month, dt.datetime(1970,1,1))
                indices['time'] = int(np.argwhere(nc.variables[time][:]==days_since))
        elif nc.variables[time].units=='Month' or nc.variables[time].units=='Months':
            #this works if the time is 1:12 for the months of the year
            indices['time'] = month-1
        else:
            raise AttributeError("time units in netcdf are not 'days since 1970-01-01' or 'Month'")
    elif time in nc.variables and len(nc.variables[time][:])==1:
        indices['time'] = 0

    #get vert index
    if vert in nc.variables and len(nc.variables[vert][:])==1:
        indices['vert'] = 0
    elif vert in nc.variables and len(nc.variables[vert][:])>1:
        deps = nc.variables['vert'][:]
        if min(nc.variables['vert'])<0:
            dep_target=dep_target*-1
        indices['vert'] = np.argmin(np.abs(deps-dep_target))

    #get data at correct index
    ix_order = [indices[k] for k in index_order]
    subset = nc.variables[z_variable]
    for ix in ix_order:
        subset = subset[ix]
    #if is masked value, and int of data is -9999, change to NaN
    if np.ma.is_masked(subset):
        if int(subset.data)==-9999:
            subset = np.nan
    print "retrieving index %s, month is %s: %s" % (ix_order, month, subset)
    return subset

#function to create new satellite column
def create_satellite_column(env, nc, z_var, year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled',
                            latitude='latitude', longitude='longitude', lat="lat",lon="lon",time="time", vert="vert",
                            index_order=['time', 'lat', 'lon']):
    column = []
    for ix, row in env.iterrows():
        year = row[year_sampled]
        month = row[month_sampled]
        lat_target = row[latitude]
        lon_target = row[longitude]
        dep_target = row[depth_sampled]

        answer = pull_satellite_data(nc, year=year, month=month, z_variable=z_var, lat_target=lat_target, lon_target=lon_target, dep_target=dep_target,
                                    lat=lat,lon=lon,time=time, vert=vert, index_order=index_order)
        column.append(answer)
    return column

#get chl satellite data
chl_nc = netCDF4.Dataset("satellite_data/chla_monthly_compiled.nc")
chl = create_satellite_column(env=env, nc=chl_nc, z_var="chlor_a")
env['chl_satellite'] = chl

#sst
print "processing sst..."
sst_nc = netCDF4.Dataset("satellite_data/sst_compiled.nc")
sst = create_satellite_column(env=env, nc=sst_nc, z_var="sst")
env['sst_satellite'] = sst
#par
print "processing par..."
par_nc = netCDF4.Dataset("satellite_data/par_monthly_compiled.nc")
par = create_satellite_column(env=env, nc=par_nc, z_var="par")
env['par_satellite'] = par
#pic
print "processing pic..."
pic_nc = netCDF4.Dataset("satellite_data/pic_monthly_compiled.nc")
pic = create_satellite_column(env=env, nc=pic_nc, z_var="pic")
env['pic_satellite'] = pic
#poc
print "processing poc..."
poc_nc = netCDF4.Dataset("satellite_data/poc_monthly_compiled.nc")
poc = create_satellite_column(env=env, nc=poc_nc, z_var="poc")
env['poc_satellite'] = poc
#npp
print "processing npp..."
npp_nc = netCDF4.Dataset("satellite_data/npp_monthly_compiled.nc")
npp = create_satellite_column(env=env, nc=npp_nc, z_var="npp")
env['npp_satellite'] = npp


# ## Add historical satellite data
print "adding calcite..."
calcite_nc = netCDF4.Dataset("historical_data/calciteAnmeanBiooracle.nc")
calcite = create_satellite_column(env=env, nc=calcite_nc, z_var='calcite', year_sampled='year_sampled', month_sampled='month_sampled',
                            latitude='latitude', longitude='longitude', lat="lat",lon="lon",time="time", vert="vert",
                            index_order=['time', 'lat', 'lon'])
env['calcite_historical'] = calcite

print "adding chlorophyll monthly historical data..."
chlormo_nc = netCDF4.Dataset("historical_data/chloMomeanNASA.nc")
chlormo = create_satellite_column(env=env, nc=chlormo_nc, z_var='Chlorophyll_Concentration',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon', 'vert', 'time'])
env['chla_monthly_historical'] = chlormo

print "adding chlorophyll annual historical data..."
chloran_nc = netCDF4.Dataset("historical_data/chlorAnmeanBiooracle.nc")
chloran = create_satellite_column(env=env, nc=chloran_nc, z_var='chlor',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','lat', 'lon'])
env['chla_annual_historical'] = chloran

print "adding chlorophyll annual range historical data..."
chlrange_nc = netCDF4.Dataset("historical_data/chlorAnrangeBiooracle.nc")
chlrange = create_satellite_column(env=env, nc=chlrange_nc, z_var='chlorrange',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','lat', 'lon'])
env['chla_annualrange_historical'] = chlrange

print "adding cloud fraction annual historical data..."
cloud_nc = netCDF4.Dataset("historical_data/cldAnmeanBiooracle.nc")
cloudan = create_satellite_column(env=env, nc=cloud_nc, z_var='cld',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','lat', 'lon'])
env['cloudfraction_annual_historical'] = cloudan

print "add cloud fraction monthly historical data..."
cloudmo_nc = netCDF4.Dataset("historical_data/cloudfracMomeanNASA.nc")
cloudmo = create_satellite_column(env=env, nc=cloudmo_nc, z_var='cloud_fraction',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat','lon','vert','time'])
env['cloudfraction_monthly_historical'] = cloudmo

print "adding cloud fraction std dev historical data..."
cloudsd_nc = netCDF4.Dataset("historical_data/cloudfracStdevNASA.nc")
cloudsd = create_satellite_column(env=env, nc=cloudsd_nc, z_var='AnnualStdev_Cloud_Fraction',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon'])
env['cloudfraction_annualstdev_historical'] = cloudsd

print "diffuse attenuation historical data..."
da_nc = netCDF4.Dataset("historical_data/daAnmeanBiooracle.nc")
da = create_satellite_column(env=env, nc=da_nc, z_var='Diffuse_attenuation_coefficient',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon'])
env['diffuseattenuation_annual_historical'] = da

print "adding day length monthly historical..."
day_nc = netCDF4.Dataset("historical_data/daylengthMomeanEarthtools.nc")
day = create_satellite_column(env=env, nc=day_nc, z_var='Day_Length_on_15th_Day_of_Month',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon', 'vert', 'time'])
env['daylength_monthly_historical'] = day

print "dust flux annual historical data..."
dustan_nc = netCDF4.Dataset("historical_data/dustAnmeanJickells.nc")
dustan = create_satellite_column(env=env, nc=dustan_nc, z_var='Dust_Deposition',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon'])
env['dustflux_annual_historical'] = dustan

print "dust flux monthly historical..."
dustmo_nc = netCDF4.Dataset("historical_data/dustMomeanJickells.nc")
dustmo = create_satellite_column(env=env, nc=dustmo_nc, z_var='Dust_Deposition',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon','vert','time'])
env['dustflux_monthly_historical'] = dustmo

print "dust std dev historical data..."
dustsd_nc = netCDF4.Dataset("historical_data/dustStdevJickells.nc")
dustsd = create_satellite_column(env=env, nc=dustsd_nc, z_var='AnnualStdev_Dust_Deposition',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon'])
env['dustflux_annualstdev_historical'] = dustsd

print "solar insolation annual historical data..."
insolan_nc = netCDF4.Dataset("historical_data/insolationAnmeanBiooracle.nc")
insolan = create_satellite_column(env=env, nc=insolan_nc, z_var='Solar_Insolation',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon'])
env['solarinsolation_annual_historical'] = insolan

print "solar insolation monthly historical data..."
insolmo_nc = netCDF4.Dataset("historical_data/insolationMomeanNASA.nc")
insolmo = create_satellite_column(env=env, nc=insolmo_nc, z_var='solar_insolation',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon', 'vert', 'time'])
env['solarinsolation_monthly_historical'] = insolmo

print "solar insolation std dev historical data..."
insolsd_nc = netCDF4.Dataset("historical_data/insolationStdevNASA.nc")
insolsd = create_satellite_column(env=env, nc=insolsd_nc, z_var='AnnualStdev_Solar_Insolation',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon'])
env['solarinsolation_annualstdev_historical'] = insolsd

print "distance from land historical data..."
land_nc = netCDF4.Dataset("historical_data/landdistAnmeanReady.nc")
land = create_satellite_column(env=env, nc=land_nc, z_var='LandDist',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon'])
env['distfromland_annual_historical'] = land

print "pycnocline depth historical..."
pyc_nc = netCDF4.Dataset("historical_data/mixedlayerdensityMomeanMontegut.nc")
pyc = create_satellite_column(env=env, nc=pyc_nc, z_var='Mixed_Layer_Depth',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon', 'vert', 'time'])
env['pycnoclinedepth_monthly_historical'] = pyc

print "pycnocline depth std dev historical..."
pycsd_nc = netCDF4.Dataset("historical_data/mixedlayerdensityStdevMontegut.nc")
pycsd = create_satellite_column(env=env, nc=pycsd_nc, z_var='AnnualStdev_Mixed_Layer_Depth_02',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon'])
env['pycnoclinedepth_annualstdev_historical'] = pycsd

print "thermocline depth monthly historical..."
therm_nc = netCDF4.Dataset("historical_data/mixedlayertempMomeanMontegut.nc")
therm = create_satellite_column(env=env, nc=therm_nc, z_var='Mixed_Layer_Depth',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon', 'vert', 'time'])
env['thermoclinedepth_monthly_historical'] = therm

print "thermocline depth std dev historical data..."
thermsd_nc = netCDF4.Dataset("historical_data/mixedlayertempStdevMontegut.nc")
thermsd = create_satellite_column(env=env, nc=thermsd_nc, z_var='AnnualStdev_Mixed_Layer_Depth_11',
                                  year_sampled='year_sampled', month_sampled='month_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon'])
env['thermoclinedepth_annualstdev_historical'] = thermsd

print "nitrate annual historical..."
nitrate_nc = netCDF4.Dataset("historical_data/nitrateAnmeanWOA.nc")
nitrate = create_satellite_column(env=env, nc=nitrate_nc, z_var='nitrate',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','vert','lat', 'lon'])
env['nitrate_annual_historical'] = nitrate

print "nitrate monthly historical..."
nitratemo_nc = netCDF4.Dataset("historical_data/nitrateMomeanWOA.nc")
nitratemo = create_satellite_column(env=env, nc=nitratemo_nc, z_var='nitrate', year_sampled='year_sampled', month_sampled='month_sampled',
                                  depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','vert','lat', 'lon'])
env['nitrate_monthly_historical'] = nitratemo

print "np ratio monthly historical..."
np_nc = netCDF4.Dataset("historical_data/npratioMomeanWOA.nc")
np = create_satellite_column(env=env, nc=np_nc, z_var='Nitrate_Phosphate_Ratio',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon', 'vert', 'time'])
env['npratio_monthly_historical'] = np

print "ocean depth historical..."
oceandep_nc = netCDF4.Dataset("historical_data/oceandepthAnmeanNASA.nc")
oceandep = create_satellite_column(env=env, nc=oceandep_nc, z_var='DepthMean',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon'])
env['oceandepth_historical'] = oceandep

print "dissolved oxygen annual historical..."
oxyan_nc = netCDF4.Dataset("historical_data/oxygendissolvedAnmeanWOA.nc")
oxyan = create_satellite_column(env=env, nc=oxyan_nc, z_var='oxygendissolved',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','vert','lat', 'lon'])
env['oxygendissolved_annual_historical'] = oxyan

print "dissolved oxygen monthly historical..."
oxymo_nc = netCDF4.Dataset("historical_data/oxygendissolvedMomeanWOA.nc")
oxymo = create_satellite_column(env=env, nc=oxymo_nc, z_var='oxygendissolved',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','vert','lat', 'lon'])
env['oxygendissolved_monthly_historical'] = oxymo

print "oxygen saturation annual historical..."
oxysatan_nc = netCDF4.Dataset("historical_data/oxygensaturationAnmeanWOA.nc")
oxysatan = create_satellite_column(env=env, nc=oxysatan_nc, z_var='oxygensaturation',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','vert','lat', 'lon'])
env['oxygensaturation_annual_historical'] = oxysatan

print "oxygen saturation monthly historical..."
oxysatmo_nc = netCDF4.Dataset("historical_data/oxygensaturationMomeanWOA.nc")
oxysatmo = create_satellite_column(env=env, nc=oxysatmo_nc, z_var='oxygensaturation',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','vert','lat', 'lon'])
env['oxygensaturation_monthly_historical'] = oxysatmo

print "apparent oxygen utilization annual historical..."
aouan_nc = netCDF4.Dataset("historical_data/oxygenutilizationAnmeanWOA.nc")
aouan = create_satellite_column(env=env, nc=aouan_nc, z_var='oxygenutilization',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','vert','lat', 'lon'])
env['oxygenutilization_annual_historical'] = aouan

print "apparent oxygen utilization monthly historical..."
aoumo_nc = netCDF4.Dataset("historical_data/oxygenutilizationMomeanWOA.nc")
aoumo = create_satellite_column(env=env, nc=aoumo_nc, z_var='oxygenutilization',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','vert','lat', 'lon'])
env['oxygenutilization_monthly_historical'] = aoumo

print "PAR annual historical..."
par_nc = netCDF4.Dataset("historical_data/parAnmeanBiooracle.nc")
par = create_satellite_column(env=env, nc=par_nc, z_var='par',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','lat', 'lon'])
env['par_annual_historical'] = par

print "pH annual historical..."
ph_nc = netCDF4.Dataset("historical_data/phAnmeanBiooracle.nc")
ph = create_satellite_column(env=env, nc=ph_nc, z_var='pH',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon'])
env['ph_annual_historical'] = ph

print "phosphate annual historical..."
phosan_nc = netCDF4.Dataset("historical_data/phosphateAnmeanWOA.nc")
phosan = create_satellite_column(env=env, nc=phosan_nc, z_var='phosphate',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','vert','lat', 'lon'])
env['phosphate_annual_historical'] = phosan

print "phosphate monthly historical..."
phosmo_nc = netCDF4.Dataset("historical_data/phosphateMomeanWOA.nc")
phosmo = create_satellite_column(env=env, nc=phosmo_nc, z_var='phosphate',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','vert','lat', 'lon'])
env['phosphate_monthly_historical'] = phosmo

print "salinity annual historical..."
salan_nc = netCDF4.Dataset("historical_data/salinityAnmeanWOA.nc")
salan = create_satellite_column(env=env, nc=salan_nc, z_var='salinity',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','vert','lat', 'lon'])
env['salinity_annual_historical'] = salan

print "salinity monthly historical..."
salmo_nc = netCDF4.Dataset("historical_data/salinityMomeanWOA.nc")
salmo = create_satellite_column(env=env, nc=salmo_nc, z_var='salinity',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','vert','lat', 'lon'])
env['salinity_monthly_historical'] = salmo

print "silicate annual historical..."
silicatean_nc = netCDF4.Dataset("historical_data/silicateAnmeanWOA.nc")
silicatean = create_satellite_column(env=env, nc=silicatean_nc, z_var='silicate',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','vert','lat', 'lon'])
env['silicate_annual_historical'] = silicatean

print "silicate monthly historical..."
silicatemo_nc = netCDF4.Dataset("historical_data/silicateMomeanWOA.nc")
silicatemo = create_satellite_column(env=env, nc=silicatemo_nc, z_var='silicate',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','vert','lat', 'lon'])
env['silicate_monthly_historical'] = silicatemo

print "SST annual historical..."
sstan_nc = netCDF4.Dataset("historical_data/sstAnmeanBiooracle.nc")
sstan = create_satellite_column(env=env, nc=sstan_nc, z_var='Sea_surface_temperature_mean',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['lat', 'lon'])
env['sst_annual_historical'] = sstan

print "ocean temperature monthly historical..."
tempmo_nc = netCDF4.Dataset("historical_data/temperatureMomeanWOA.nc")
tempmo = create_satellite_column(env=env, nc=tempmo_nc, z_var='temperature',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','vert','lat', 'lon'])
env['oceantemp_monthly_historical'] = tempmo

print "ocean temperature annual historical..."
tempan_nc = netCDF4.Dataset("historical_data/watertempAnmeanWOA.nc")
tempan = create_satellite_column(env=env, nc=tempan_nc, z_var='watertemp',
                                  year_sampled='year_sampled', month_sampled='month_sampled', depth_sampled='depth_sampled', latitude='latitude', longitude='longitude',
                                  lat="lat",lon="lon",time="time", vert="vert",
                                  index_order=['time','vert','lat', 'lon'])
env['oceantemp_annual_historical'] = tempan


#save file
print "saving remote data!"
print env.head()
env.to_csv("env_remote_data_TARA.csv", index=False)
