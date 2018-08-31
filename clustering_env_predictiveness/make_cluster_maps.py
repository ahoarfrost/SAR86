from MapProjectionsFns import *
import netCDF4 as nc

projections = nc.Dataset("envforprojections_2009.nc", "r+", format="NETCDF4")
contemp = projections['contemporary']
monthly = projections['monthly']
annual = projections['annual']

#reorder into expected order 
general_feats = ['latitude', 'longitude', 'depth_sampled']
satellite_feats = ['chl_satellite', 'sst_satellite', 'par_satellite', 'pic_satellite', 'poc_satellite', 'npp_satellite']
monthly_feats = ['chla_monthly_historical', 'cloudfraction_monthly_historical', 'daylength_monthly_historical', 
                 'dustflux_monthly_historical', 'solarinsolation_monthly_historical', 'pycnoclinedepth_monthly_historical',
                 'thermoclinedepth_monthly_historical', 'nitrate_monthly_historical', 'npratio_monthly_historical',
                 'oxygendissolved_monthly_historical', 'oxygensaturation_monthly_historical', 'oxygenutilization_monthly_historical',
                 'phosphate_monthly_historical', 'salinity_monthly_historical', 'silicate_monthly_historical', 'oceantemp_monthly_historical'
                ]
annual_feats = ['chla_annual_historical', 'chla_annualrange_historical', 'cloudfraction_annual_historical', 'cloudfraction_annualstdev_historical',
                'diffuseattenuation_annual_historical', 'par_annual_historical', 'salinity_annual_historical',
                'thermoclinedepth_annualstdev_historical', 'nitrate_annual_historical', 'solarinsolation_annual_historical',
                'distfromland_annual_historical', 'oxygendissolved_annual_historical', 'sst_annual_historical', 
                'pycnoclinedepth_annualstdev_historical', 'solarinsolation_annualstdev_historical', 'oceandepth_historical', 
                'dustflux_annual_historical', 'oxygensaturation_annual_historical', 'dustflux_annualstdev_historical', 
                'oxygenutilization_annual_historical', 'phosphate_annual_historical', 'silicate_annual_historical', 
                'calcite_annual_historical', 'oceantemp_annual_historical', 'ph_annual_historical'
                ]

feat_order = general_feats+satellite_feats+monthly_feats+annual_feats+['intercept']

# create map projections
projnc = nc.Dataset("predictiveness_clusters_projections_k5.nc", "w", format="NETCDF4")

def add_projection(contemp, monthly, annual, 
                   group_name, centroids, 
                   depth_sampled=-0.556543241, satellite_month=0, monthly_month=0, 
                   satellite_feats=satellite_feats, monthly_feats=monthly_feats, annual_feats=annual_feats):

    #create k6 group in nc
    print("creating netcdf group",group_name,"...")
    kGroup = projnc.createGroup(group_name)

    #add dimensions - 
    #contemporary dims should be lat, lon, and time which will be unlimited (but for now 4 layers with quarterly months (Jan 2009, Apr 2009, July 2009, Oct 2009))
    lat = kGroup.createDimension("lat", 2160)
    lon = kGroup.createDimension("lon", 4320)

    #create lat and lon and time variables as 64-bit floats and ints
    latitudes = kGroup.createVariable("lat","f8",("lat",))
    longitudes = kGroup.createVariable("lon","f8",("lon",))

    #assign lat/lon values at 9km res to the lat/lon variables (same as contemp)
    kGroup['lat'][:] = contemp['lat'][:]
    kGroup['lon'][:] = contemp['lon'][:]

    #extract predictions and add to nc
    print("creating feature matrices...")
    feat_matrix = get_feat_matrices(contemp, monthly, annual, 
                                    depth_sampled=depth_sampled, satellite_month=satellite_month, 
                                    monthly_month=monthly_month, satellite_feats=satellite_feats, 
                                    monthly_feats=monthly_feats, annual_feats=annual_feats)
    
    for clust in centroids.columns:
        print('---processing', clust, '---')
        print("scoring the matrix...")
        testgene_scores = get_score_matrix(feat_matrix, weights=centroids[clust])
        print("squishing through the sigmoid...")
        sigmoids = get_sigmoid_matrix(testgene_scores)
        print("Eureka!")
        
        #create a variable in our nc file for the sig if it doesn't exist, else just add data
        if str(clust) in kGroup.variables:
            kGroup[str(clust)][:] = sigmoids
        else:
            kGroup.createVariable(str(clust), "f8", ("lat", "lon"))
            kGroup[str(clust)][:] = sigmoids

            
## projections from k5 centroids, all surface
centroids5 = pd.read_csv("centroidsk5.csv", index_col="feature")

#Jan 2009
add_projection(contemp, monthly, annual, 
               group_name="Jan2009_k5", centroids=centroids5, 
               depth_sampled=-0.556543241, satellite_month=0, monthly_month=0, 
               satellite_feats=satellite_feats, monthly_feats=monthly_feats, annual_feats=annual_feats)

#Apr 2009
add_projection(contemp, monthly, annual, 
               group_name="Apr2009_k5", centroids=centroids5, 
               depth_sampled=-0.556543241, satellite_month=1, monthly_month=3, 
               satellite_feats=satellite_feats, monthly_feats=monthly_feats, annual_feats=annual_feats)

#July 2009
add_projection(contemp, monthly, annual, 
               group_name="Jul2009_k5", centroids=centroids5, 
               depth_sampled=-0.556543241, satellite_month=2, monthly_month=6, 
               satellite_feats=satellite_feats, monthly_feats=monthly_feats, annual_feats=annual_feats)

#Oct 2009
add_projection(contemp, monthly, annual, 
               group_name="Oct2009_k5", centroids=centroids5, 
               depth_sampled=-0.556543241, satellite_month=3, monthly_month=9, 
               satellite_feats=satellite_feats, monthly_feats=monthly_feats, annual_feats=annual_feats)

projnc.close()

'''
#to tune the number of clusters, I also made projections for k4, k6, k7, k9 and visually inspected the most important features to identify which ones make oceanographic sense and are not redundant with each other. k5 looks good doing this. Code for how to do this below.
## k6
centroids6 = pd.read_csv("centroidsk6.csv", index_col="feature")

add_projection(contemp, monthly, annual, 
               group_name="Jan2009_k6", centroids=centroids6, 
               depth_sampled=-0.556543241, satellite_month=0, monthly_month=0, 
               satellite_feats=satellite_feats, monthly_feats=monthly_feats, annual_feats=annual_feats)


## k4
centroids4 = pd.read_csv("centroidsk4.csv", index_col="feature")

add_projection(contemp, monthly, annual, 
               group_name="Jan2009_k4", centroids=centroids4, 
               depth_sampled=-0.556543241, satellite_month=0, monthly_month=0, 
               satellite_feats=satellite_feats, monthly_feats=monthly_feats, annual_feats=annual_feats)

## k7
centroids7 = pd.read_csv("centroidsk7.csv", index_col="feature")

add_projection(contemp, monthly, annual, 
               group_name="Jan2009_k7", centroids=centroids7, 
               depth_sampled=-0.556543241, satellite_month=0, monthly_month=0, 
               satellite_feats=satellite_feats, monthly_feats=monthly_feats, annual_feats=annual_feats)

## k9
centroids9 = pd.read_csv("centroidsk9.csv", index_col="feature")

add_projection(contemp, monthly, annual, 
               group_name="Jan2009_k9", centroids=centroids9, 
               depth_sampled=-0.556543241, satellite_month=0, monthly_month=0, 
               satellite_feats=satellite_feats, monthly_feats=monthly_feats, annual_feats=annual_feats)

#look at top 5 most impt env features
for clust in centroids6.columns:
    print(clust)
    centroids6['abs'] = np.abs(centroids6[clust])
    print(centroids6.sort_values('abs', ascending=False)[clust][0:5])

for clust in centroids5.columns:
    print(clust)
    centroids5['abs'] = np.abs(centroids5[clust])
    print(centroids5.sort_values('abs', ascending=False)[clust][0:5])

for clust in centroids7.columns:
    print(clust)
    centroids7['abs'] = np.abs(centroids7[clust])
    print(centroids7.sort_values('abs', ascending=False)[clust][0:5])

for clust in centroids4.columns:
    print(clust)
    centroids4['abs'] = np.abs(centroids4[clust])
    print(centroids4.sort_values('abs', ascending=False)[clust][0:5])

'''


