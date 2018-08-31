import pandas as pd
import numpy as np
from sklearn.cluster import k_means
#import plotly as py
#import plotly.graph_objs as go
#py.offline.init_notebook_mode()
#import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')


## read in data
print("reading in model data...")
acc2080 = pd.read_csv("../satellite_models/Balanced2080Genes_LogRegAccuracies_L1.csv")
coef2080 = pd.read_csv("../satellite_models/Balanced2080Genes_LogRegCoefficients_L1.csv", index_col="feature")
#need features to be cols, genes want to cluster to be rows for sklearn, so transform
Tcoef = coef2080.transpose()

#define feat order so can reorder into expected order 
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

# train final clusters with k=6
#train and save centroids and cluster labels for diff k
print("training clusters for k =", 5)
clusters = k_means(Tcoef,n_clusters=5, random_state=1423)
cluster_names = ['cluster'+str(num) for num in range(1,6)]

centroids = pd.DataFrame(clusters[0],columns=Tcoef.columns,index=cluster_names)
#need to transpose so each column is a vector of feats
centroids = centroids.transpose()
centroids = centroids.reindex(feat_order)

labels = pd.DataFrame({"gene":Tcoef.index, "cluster":clusters[1]+1})

print("saving cluster centroids and labels...")
centroids.to_csv("centroidsk5.csv")
labels.to_csv("clusterlabelsk5.csv", index=False)
print("...done!")


'''
# choose number clusters - going with k=6
#cluster sites;
n_clusters = []
inertias = []
dists = []
for cluster_candidate in range(1,20):
    print("training for k =",cluster_candidate)
    n_clusters.append(cluster_candidate)
    clusters = k_means(Tcoef,n_clusters=cluster_candidate, random_state=1423)
    inertias.append(clusters[2])
    
    from scipy.spatial.distance import squareform, pdist
    centroids = pd.DataFrame(clusters[0],columns=Tcoef.columns)
    centroids = centroids.transpose()
    dist = pd.DataFrame(squareform(pdist(centroids.T)), columns=centroids.columns, index=centroids.columns)
    dists.append(np.mean(np.mean(dist)))

trace1 = go.Scatter(
    x = n_clusters,
    y = inertias,
    mode = 'lines+markers'
    )
py.offline.iplot([trace1], filename='inertias')

trace1 = go.Scatter(
    x = n_clusters,
    y = dists,
    mode = 'lines+markers'
    )
py.offline.iplot([trace1], filename='inertias')

# look at k=4, k=6, k=7, k=9 more closely 
## train clusters at k=4,5,6,7,9

#train and save centroids and cluster labels for diff k
for k in [4,5,6,7,9]:
    print("training clusters for k =", k)
    clusters = k_means(Tcoef,n_clusters=k, random_state=1423)
    cluster_names = ['cluster'+str(num) for num in range(1,k+1)]
    
    centroids = pd.DataFrame(clusters[0],columns=Tcoef.columns,index=cluster_names)
    #need to transpose so each column is a vector of feats
    centroids = centroids.transpose()
    centroids = centroids.reindex(feat_order)
    
    labels = pd.DataFrame({"gene":Tcoef.index, "cluster":clusters[1]+1})
    
    centroids.to_csv("centroidsk"+str(k)+".csv")
    labels.to_csv("clusterlabelsk"+str(k)+".csv", index=False)

#read these in separately

centroids4 = pd.read_csv("centroidsk4.csv", index_col="feature")
centroids5 = pd.read_csv("centroidsk5.csv", index_col="feature")
centroids6 = pd.read_csv("centroidsk6.csv", index_col="feature")
centroids7 = pd.read_csv("centroidsk7.csv", index_col="feature")
centroids9 = pd.read_csv("centroidsk9.csv", index_col="feature")

labels4 = pd.read_csv("clusterlabelsk4.csv")
labels5 = pd.read_csv("clusterlabelsk5.csv")
labels6 = pd.read_csv("clusterlabelsk6.csv")
labels7 = pd.read_csv("clusterlabelsk7.csv")
labels9 = pd.read_csv("clusterlabelsk9.csv")

## cluster evenness
from collections import Counter
print(Counter(labels['cluster']))
plt.hist(labels4['cluster'], bins=4)
plt.figure()
plt.hist(labels6['cluster'], bins=6)
plt.figure()
plt.hist(labels7['cluster'], bins=7)
plt.figure()
plt.hist(labels9['cluster'], bins=9)
plt.show()

## heatmap of centroids
#k4
centroids4_heatmap = centroids4.reindex(feat_order)

trace = go.Heatmap(z=[list(centroids4_heatmap[column][1:]) for column in centroids4_heatmap.columns],
                   x=centroids4_heatmap.index[1:],
                   y=centroids4_heatmap.columns,
                   colorscale=[[0.0, 'rgb(40,0,255)'],[0.4, 'rgb(255,255,255)'], [0.6, 'rgb(255,255,255)'], [1.0, 'rgb(255,0,40)']],
                   zmin=-2,
                   zmax=2
                  )
data=[trace]
py.offline.iplot(data, filename='cluster4-heatmap')

#k5
centroids5_heatmap = centroids5.reindex(feat_order)

trace = go.Heatmap(z=[list(centroids5_heatmap[column][1:]) for column in centroids5_heatmap.columns],
                   x=centroids5_heatmap.index[1:],
                   y=centroids5_heatmap.columns,
                   colorscale=[[0.0, 'rgb(40,0,255)'],[0.4, 'rgb(255,255,255)'], [0.6, 'rgb(255,255,255)'], [1.0, 'rgb(255,0,40)']],
                   zmin=-2,
                   zmax=2
                  )
data=[trace]
py.offline.iplot(data, filename='cluster4-heatmap')

#k6
centroids6_heatmap = centroids6.reindex(feat_order)

trace = go.Heatmap(z=[list(centroids6_heatmap[column][1:]) for column in centroids6_heatmap.columns],
                   x=centroids6_heatmap.index[1:],
                   y=centroids6_heatmap.columns,
                   colorscale=[[0.0, 'rgb(40,0,255)'],[0.4, 'rgb(255,255,255)'], [0.6, 'rgb(255,255,255)'], [1.0, 'rgb(255,0,40)']],
                   zmin=-2,
                   zmax=2
                  )
data=[trace]
py.offline.iplot(data, filename='cluster6-heatmap')

#k7
centroids7_heatmap = centroids7.reindex(feat_order)

trace = go.Heatmap(z=[list(centroids7_heatmap[column][1:]) for column in centroids7_heatmap.columns],
                   x=centroids7_heatmap.index[1:],
                   y=centroids7_heatmap.columns,
                   colorscale=[[0.0, 'rgb(40,0,255)'],[0.4, 'rgb(255,255,255)'], [0.6, 'rgb(255,255,255)'], [1.0, 'rgb(255,0,40)']],
                   zmin=-2,
                   zmax=2
                  )
data=[trace]
py.offline.iplot(data, filename='cluster7-heatmap')

#k9
centroids9_heatmap = centroids9.reindex(feat_order)

trace = go.Heatmap(z=[list(centroids9_heatmap[column][1:]) for column in centroids9_heatmap.columns],
                   x=centroids9_heatmap.index[1:],
                   y=centroids9_heatmap.columns,
                   colorscale=[[0.0, 'rgb(40,0,255)'],[0.4, 'rgb(255,255,255)'], [0.6, 'rgb(255,255,255)'], [1.0, 'rgb(255,0,40)']],
                   zmin=-2,
                   zmax=2
                  )
data=[trace]
py.offline.iplot(data, filename='cluster9-heatmap')

## compare cluster similarities of clusters across diff k

#k4 vs k7
from scipy.spatial import distance
clustdist47 = pd.DataFrame()
for col4 in centroids4.columns:
    for col7 in centroids7.columns:
        dist = distance.euclidean(centroids4[col4],centroids7[col7])
        clustdist47.loc[col4,col7] = dist
clustdist47.columns = centroids7.columns
clustdist47.index = centroids4.columns

trace = go.Heatmap(z=[list(clustdist47[column]) for column in clustdist47.columns],
                   x=clustdist47.index,
                   y=clustdist47.columns,
                   colorscale=[[0.0, 'rgb(255,255,255)'],[1.0, 'rgb(255,0,40)']],
                   zmin=0,
                   zmax=4
                  )
data=[trace]
py.offline.iplot(data, filename='cluster-heatmap')

#k4 vs k9
from scipy.spatial import distance
clustdist49 = pd.DataFrame()
for col4 in centroids4.columns:
    for col9 in centroids9.columns:
        dist = distance.euclidean(centroids4[col4],centroids9[col9])
        clustdist49.loc[col4,col9] = dist
clustdist49.columns = centroids9.columns
clustdist49.index = centroids4.columns

trace = go.Heatmap(z=[list(clustdist49[column]) for column in clustdist49.columns],
                   x=clustdist49.index,
                   y=clustdist49.columns,
                   colorscale=[[0.0, 'rgb(255,255,255)'],[1.0, 'rgb(255,0,40)']],
                   zmin=0,
                   zmax=4
                  )
data=[trace]
py.offline.iplot(data, filename='cluster-heatmap')

#k7 vs k9
from scipy.spatial import distance
clustdist79 = pd.DataFrame()
for col7 in centroids7.columns:
    for col9 in centroids9.columns:
        dist = distance.euclidean(centroids7[col7],centroids9[col9])
        clustdist79.loc[col7,col9] = dist
clustdist79.columns = centroids9.columns
clustdist79.index = centroids7.columns

trace = go.Heatmap(z=[list(clustdist79[column]) for column in clustdist79.columns],
                   x=clustdist79.index,
                   y=clustdist79.columns,
                   colorscale=[[0.0, 'rgb(255,255,255)'],[1.0, 'rgb(255,0,40)']],
                   zmin=0,
                   zmax=4
                  )
data=[trace]
py.offline.iplot(data, filename='cluster-heatmap')

#k5 vs k7
from scipy.spatial import distance
clustdist75 = pd.DataFrame()
for col7 in centroids7.columns:
    for col5 in centroids5.columns:
        dist = distance.euclidean(centroids7[col7],centroids5[col5])
        clustdist75.loc[col7,col5] = dist
clustdist75.columns = centroids5.columns
clustdist75.index = centroids7.columns

trace = go.Heatmap(z=[list(clustdist75[column]) for column in clustdist75.columns],
                   x=clustdist75.index,
                   y=clustdist75.columns,
                   colorscale=[[0.0, 'rgb(255,255,255)'],[1.0, 'rgb(255,0,40)']],
                   zmin=0,
                   zmax=4
                  )
data=[trace]
py.offline.iplot(data, filename='cluster-heatmap')

#k5 vs k6
from scipy.spatial import distance
clustdist65 = pd.DataFrame()
for col6 in centroids6.columns:
    for col5 in centroids5.columns:
        dist = distance.euclidean(centroids6[col6],centroids5[col5])
        clustdist65.loc[col6,col5] = dist
clustdist65.columns = centroids5.columns
clustdist65.index = centroids6.columns

trace = go.Heatmap(z=[list(clustdist65[column]) for column in clustdist65.columns],
                   x=clustdist65.index,
                   y=clustdist65.columns,
                   colorscale=[[0.0, 'rgb(255,255,255)'],[1.0, 'rgb(255,0,40)']],
                   zmin=0,
                   zmax=4
                  )
data=[trace]
py.offline.iplot(data, filename='cluster-heatmap')

# five and six clusters look pretty different still. two k5 clusters being split in ~two

#k7 vs k6
from scipy.spatial import distance
clustdist67 = pd.DataFrame()
for col6 in centroids6.columns:
    for col7 in centroids7.columns:
        dist = distance.euclidean(centroids6[col6],centroids7[col7])
        clustdist67.loc[col6,col7] = dist
clustdist67.columns = centroids7.columns
clustdist67.index = centroids6.columns

trace = go.Heatmap(z=[list(clustdist67[column]) for column in clustdist67.columns],
                   x=clustdist67.index,
                   y=clustdist67.columns,
                   colorscale=[[0.0, 'rgb(255,255,255)'],[1.0, 'rgb(255,0,40)']],
                   zmin=0,
                   zmax=4
                  )
data=[trace]
py.offline.iplot(data, filename='cluster-heatmap')


# difference between 6 and 7 clusters very little, really just cluster 3 is being split up and only a little

#k4 vs k6
from scipy.spatial import distance
clustdist64 = pd.DataFrame()
for col6 in centroids6.columns:
    for col4 in centroids4.columns:
        dist = distance.euclidean(centroids6[col6],centroids4[col4])
        clustdist64.loc[col6,col4] = dist
clustdist64.columns = centroids4.columns
clustdist64.index = centroids6.columns

trace = go.Heatmap(z=[list(clustdist64[column]) for column in clustdist64.columns],
                   x=clustdist64.index,
                   y=clustdist64.columns,
                   colorscale=[[0.0, 'rgb(255,255,255)'],[1.0, 'rgb(255,0,40)']],
                   zmin=0,
                   zmax=4
                  )
data=[trace]
py.offline.iplot(data, filename='cluster-heatmap')

## compare clusters to themselves
#k4 to k4
from scipy.spatial import distance
clustdist44 = pd.DataFrame()
for col6 in centroids4.columns:
    for col4 in centroids4.columns:
        dist = distance.euclidean(centroids4[col6],centroids4[col4])
        clustdist44.loc[col6,col4] = dist
clustdist44.columns = centroids4.columns
clustdist44.index = centroids4.columns

trace = go.Heatmap(z=[list(clustdist44[column]) for column in clustdist44.columns],
                   x=clustdist44.index,
                   y=clustdist44.columns,
                   colorscale=[[0.0, 'rgb(255,255,255)'],[1.0, 'rgb(255,0,40)']],
                   zmin=0,
                   zmax=4
                  )
data=[trace]
py.offline.iplot(data, filename='cluster-heatmap')

#k5 to k5
from scipy.spatial import distance
clustdist99 = pd.DataFrame()
for col6 in centroids5.columns:
    for col4 in centroids5.columns:
        dist = distance.euclidean(centroids5[col6],centroids5[col4])
        clustdist99.loc[col6,col4] = dist
clustdist99.columns = centroids5.columns
clustdist99.index = centroids5.columns

trace = go.Heatmap(z=[list(clustdist99[column]) for column in clustdist99.columns],
                   x=clustdist99.index,
                   y=clustdist99.columns,
                   colorscale=[[0.0, 'rgb(255,255,255)'],[1.0, 'rgb(255,0,40)']],
                   zmin=0,
                   zmax=4
                  )
data=[trace]
py.offline.iplot(data, filename='cluster-heatmap')

#k6 to k6
from scipy.spatial import distance
clustdist99 = pd.DataFrame()
for col6 in centroids6.columns:
    for col4 in centroids6.columns:
        dist = distance.euclidean(centroids6[col6],centroids6[col4])
        clustdist99.loc[col6,col4] = dist
clustdist99.columns = centroids6.columns
clustdist99.index = centroids6.columns

trace = go.Heatmap(z=[list(clustdist99[column]) for column in clustdist99.columns],
                   x=clustdist99.index,
                   y=clustdist99.columns,
                   colorscale=[[0.0, 'rgb(255,255,255)'],[1.0, 'rgb(255,0,40)']],
                   zmin=0,
                   zmax=4
                  )
data=[trace]
py.offline.iplot(data, filename='cluster-heatmap')
'''
