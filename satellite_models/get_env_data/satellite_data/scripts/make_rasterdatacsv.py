#make rasterdata.csv files
import csv
import os

#define function, given asc_name, that will calculate time start and time end and make raster data csv rows
def generate_csv(prefix='/sst/sst_',suffix='.L3m_MO_SST_sst_9km.asc'):
    years = ['2008','2009','2010','2011','2012','2013']
    leap_id = [True,False,False,False,True,False]
    leap_days = [['001','031'],['032','060'],['061','091'],['092','121'],['122','152'],['153','182'],['183','213'],['214','244'],['245','274'],['275','305'],['306','335'],['336','366']]
    non_leap_days = [['001','031'],['032','059'],['060','090'],['091','120'],['121','151'],['152','181'],['182','212'],['213','243'],['244','273'],['274','304'],['305','334'],['335','365']]
    month_names = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
    month_nums = ['01','02','03','04','05','06','07','08','09','10','11','12']

    asc_list = [['ASCII_RASTER_PATH','CLIMATOLOGY_START_DATE','CLIMATOLOGY_END_DATE']]
    for year_ix,year in enumerate(years):
        if leap_id[year_ix]==True:
            days = leap_days
        else:
            days = non_leap_days
        for month_ix,month in enumerate(month_names):
            day_start = '01'
            day_end = int(days[month_ix][1])-int(days[month_ix][0])+1
            mo = month_nums[month_ix]
            asc_name = os.getcwd()+prefix+year+month_names[month_ix]+suffix
            time_start = year+'-'+mo+'-'+day_start
            time_end = year+'-'+mo+'-'+str(day_end)

            asc = [asc_name,time_start,time_end]
            asc_list.append(asc)

    return asc_list

###sst###
sst_list = generate_csv(prefix='/sst/'+'sst_',suffix='.L3m_MO_SST_sst_9km.asc')
#save sst file
with open('sst/sst_rasterdata.csv', "wb") as f:
    writer = csv.writer(f)
    writer.writerows(sst_list)
f.close()

###chl###
chl_list = generate_csv(prefix='/chla/chla_',suffix='.L3m_MO_CHL_chlor_a_4km.asc')
#save file
with open('chla/chla_rasterdata.csv', "wb") as f:
    writer = csv.writer(f)
    writer.writerows(chl_list)
f.close()

###par###
par_list = generate_csv(prefix='/par/par_',suffix='.L3m_MO_PAR_par_9km.asc')
#save file
with open('par/par_rasterdata.csv', "wb") as f:
    writer = csv.writer(f)
    writer.writerows(par_list)
f.close()

###poc###
poc_list = generate_csv(prefix='/poc/poc_',suffix='.L3m_MO_POC_poc_9km.asc')
#save file
with open('poc/poc_rasterdata.csv', "wb") as f:
    writer = csv.writer(f)
    writer.writerows(poc_list)
f.close()

###pic###
pic_list = generate_csv(prefix='/pic/pic_',suffix='.L3m_MO_PIC_pic_9km.asc')
#save file
with open('pic/pic_rasterdata.csv', "wb") as f:
    writer = csv.writer(f)
    writer.writerows(pic_list)
f.close()

###npp###
npp_list = generate_csv(prefix='/npp/npp_',suffix='..asc')
#save file
with open('npp/npp_rasterdata.csv', "wb") as f:
    writer = csv.writer(f)
    writer.writerows(npp_list)
f.close()

###globalTanomaly###
'''
tanom_list = generate_csv(prefix='/globalTanomaly/tanom_',suffix='.asc')
#save file
with open('tanomaly/tanom_rasterdata.csv', "wb") as f:
    writer = csv.writer(f)
    writer.writerows(tanom_list)
f.close()
'''

###currents####
