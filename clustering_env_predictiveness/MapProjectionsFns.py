
#functions for making map projections

#lat/lon at 9km resolution - grab it from sst
import netCDF4 as nc
import datetime as dt
import numpy as np
import pandas as pd

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


def insert_contemp(var, var_target, source_nc, target_nc):
    #create variable and attributes
    if var_target in target_nc.variables:
        var_target = target_nc[var_target]
    else:
        var_target = target_nc.createVariable(var_target,"f8",("lat","lon","time"))

    #set attributes as same as original
    for attr in source_nc[var].ncattrs():
        var_target.setncattr(attr, str(getattr(source_nc[var], attr)))
    #insert into projections
    if source_nc[var].shape == (72,2160,4320):
        for day_ix, day in enumerate(target_nc['time'][:]):
            subset = np.array(source_nc[var][int(np.argwhere(source_nc['time'][:]==target_nc['time'][day_ix])),:,:])
            if np.ma.is_masked(subset):
                print("filling masked")
                subset = subset.filled(fill_value=np.nan)
            #insert into projections
            var_target[:,:,day_ix] = subset
    else:
        print('shape not as expected')



def insert_monthly(var, var_target, source_nc, target_nc):
    from mpl_toolkits import basemap
    #create variable and attributes
    if var_target in target_nc.variables:
        var_target = target_nc[var_target]
    else:
        var_target = target_nc.createVariable(var_target,"f8",("lat","lon","time"))

    #set attributes as same as original
    for attr in source_nc[var].ncattrs():
        var_target.setncattr(attr, str(getattr(source_nc[var], attr)))
    #insert into projections

    print("processing ", var)
    #change resolution
    lats_source = source_nc['lat'][:]
    lons_source = source_nc['lon'][:]
    lats_fine = target_nc['lat'][:]
    lons_fine = target_nc['lon'][:]
    lons_sub, lats_sub = np.meshgrid(lons_fine, lats_fine)

    if source_nc[var].shape in [(360,720,1,12), (360, 720, 14, 12)]:
        for mo_ix, month in enumerate(target_nc['time'][:]):
            var_source = source_nc[var][:,:,0,mo_ix]
            var_fine = basemap.interp(var_source, lons_source, lats_source, lons_sub, lats_sub, order=1)
            if np.ma.is_masked(var_fine):
                var_fine = var_fine.filled(fill_value=np.nan)
            #insert into projections
            var_target[:,:,mo_ix] = var_fine

    elif source_nc[var].shape in [(12, 37, 180, 360), (12, 57, 180, 360)]:
        #need to use len(vert)-2 layer; 36 or 56
        vert_layer = len(source_nc['vert'][:])-2
        print('using vertical layer ', vert_layer, ' for ', var)
        for mo_ix, month in enumerate(target_nc['time'][:]):
            var_source = source_nc[var][mo_ix,vert_layer,:,:]
            var_fine = basemap.interp(var_source, lons_source, lats_source, lons_sub, lats_sub, order=1)
            if np.ma.is_masked(var_fine):
                var_fine = var_fine.filled(fill_value=np.nan)
            #insert into projections
            var_target[:,:,mo_ix] = var_fine

    else:
        print('shape not as expected')


def insert_annual(var, var_target, source_nc, target_nc):
    from mpl_toolkits import basemap
    #create variable and attributes
    if var_target in target_nc.variables:
        var_target = target_nc[var_target]
    else:
        var_target = target_nc.createVariable(var_target,"f8",("lat","lon"))

    #set attributes as same as original
    for attr in source_nc[var].ncattrs():
        var_target.setncattr(attr, str(getattr(source_nc[var], attr)))
    #insert into projections

    print("processing ", var)
    #change resolution
    lats_source = source_nc['lat'][:]
    lons_source = source_nc['lon'][:]
    lats_fine = target_nc['lat'][:]
    lons_fine = target_nc['lon'][:]
    lons_sub, lats_sub = np.meshgrid(lons_fine, lats_fine)

    if source_nc[var].shape in [(1, 2160, 4320)]:
        var_source = source_nc[var][0,:,:]
        var_fine = var_source

    elif source_nc[var].shape in [(1, 102, 180, 360)]:
        print('using level 100 for ', var)
        #gotta use len(vert)-2 level to get surface data (centered around 5m)
        var_source = source_nc[var][0,100,:,:]
        var_fine = basemap.interp(var_source, lons_source, lats_source, lons_sub, lats_sub, order=1)

    elif source_nc[var].shape in [(360, 720)]:
        var_source = source_nc[var][:,:]
        var_fine = basemap.interp(var_source, lons_source, lats_source, lons_sub, lats_sub, order=1)

    else:
        print('shape not as expected')

    if np.ma.is_masked(var_fine):
        print("filling masked")
        var_fine = var_fine.filled(fill_value=np.nan)
        #insert into projections
    var_target[:,:] = var_fine

def reorder_features(df, feat_order):
    df_reordered = df.reindex(feat_order)
    return df_reordered


def get_feat_matrices(contemp, monthly, annual, depth_sampled=-0.556543241, satellite_month=0, monthly_month=0, satellite_feats=satellite_feats, monthly_feats=monthly_feats, annual_feats=annual_feats):

    #this also requires contemp, monthly, annual,
    #satellite_feats, monthly_feats, annual_feats

    feat_matrices = []
    #extract 2160x4320 matrix for each feature in same order as weights, append to list

    ###general descriptors
    print('getting general feats')
    lats = contemp['lat'][:].reshape(2160,1)
    lat_matrix = np.repeat(lats, 4320, axis=1)
    #normalize, subtract mean and divide by std
    lat_matrix = (lat_matrix - np.nanmean(lat_matrix))/float(np.nanstd(lat_matrix))

    lons = contemp['lon'][:].reshape(1,4320)
    lon_matrix = np.repeat(lons, 2160, axis=0)
    lon_matrix = (lon_matrix - np.nanmean(lon_matrix))/float(np.nanstd(lon_matrix))

    #depth_sampled is all surface for now, stealing that value from TARA env_remote_data_preprocessed
    depth_matrix = np.broadcast_to(depth_sampled, (2160,4320))

    feat_matrices.extend([lat_matrix, lon_matrix, depth_matrix])

    ###satellite data
    print('getting satellite feats')
    #satellite_month 0 for jan, 1 for apr, 2 for july, 3 for oct
    for sat in satellite_feats:
        mat = contemp[sat][:,:,satellite_month]
        #fill if it's a masked array
        if np.ma.is_masked(mat):
            mat = mat.filled(fill_value=np.nan)
        #normalize, subtract mean and divide by std
        mat = (mat - np.nanmean(mat))/float(np.nanstd(mat))
        #add to feat_matrices
        feat_matrices.extend([mat])

    ###monthly data
    print('getting monthly historical feats')
    #monthly_month 0 for jan, 3 for apr, 6 for july, 9 for oct
    for mo in monthly_feats:
        mat = monthly[mo][:,:,monthly_month]
        #normalize, subtract mean and divide by std
        mat = (mat - np.nanmean(mat))/float(np.nanstd(mat))
        #add to feat_matrices
        feat_matrices.extend([mat])

    ###annual data
    print('getting annual historical feats')
    for ann in annual_feats:
        mat = annual[ann][:,:]
        #normalize, subtract mean and divide by std
        mat = (mat - np.nanmean(mat))/float(np.nanstd(mat))
        #add to feat_matrices
        feat_matrices.extend([mat])

    ###intercept
    print('adding intercept')
    int_mat = np.broadcast_to(1, shape=(2160,4320))
    feat_matrices.extend([int_mat])

    print('DONE!')
    return np.array(feat_matrices)



def get_score_matrix(feat_matrix, weights):
    #broadcast weights to same size each feat's matrix so can multiply them together
    broadcast_weights = np.array([np.broadcast_to(w, shape=(2160,4320)) for w in weights])
    assert broadcast_weights.shape==feat_matrix.shape, "feat and weight matrices are not same shape"

    #multiply each value in each feat matrix by its coefficient weight
    weighted_matrix = feat_matrix*broadcast_weights
    assert weighted_matrix.shape==feat_matrix.shape, "weighted feat and feat matrices are not same shape"

    #sum all 51 matrices to get score matrix
    scores = weighted_matrix.sum(axis=0)

    return scores


def get_sigmoid_matrix(score_matrix):
    #apply the sigmoid function to any scalar or array of any size
    sigmoid = 1.0/(1+np.exp(-score_matrix))
    return sigmoid


def predict_pres_abs(sigmoid_matrix):
    prediction_matrix = sigmoid_matrix.copy()
    prediction_matrix[prediction_matrix>=0.5] = 1
    prediction_matrix[prediction_matrix<0.5] = 0
    return prediction_matrix


def create_prediction_matrix(contemp, monthly, annual, gene, depth_sampled=-0.556543241, satellite_month=0, monthly_month=0, satellite_feats=satellite_feats, monthly_feats=monthly_feats, annual_feats=annual_feats):
    print("creating feature matrices...")
    feat_matrix = get_feat_matrices(contemp, monthly, annual, depth_sampled=-0.556543241, satellite_month=0, monthly_month=0, satellite_feats=satellite_feats, monthly_feats=monthly_feats, annual_feats=annual_feats)
    print("scoring the matrix...")
    testgene_scores = get_score_matrix(feat_matrix, weights=gene)
    print("squishing through the sigmoid...")
    testsig_matrix = get_sigmoid_matrix(testgene_scores)
    print("making predictions...")
    testpredictions = predict_pres_abs(testsig_matrix)
    print("Eureka!")
    return testsig_matrix, testpredictions
