import pandas as pd
import numpy as np

import plotly as py
import plotly.graph_objs as go
#py.offline.init_notebook_mode(connected=True)
import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')

from collections import Counter


# # Taxonomic Enrichment in Clusters
print("doing taxonomic enrichment analysis...")
# ## read in files
labels = pd.read_csv("../clustering_env_predictiveness/clusterlabelsk5.csv")
annotations = pd.read_csv("../annotation_SAR86.tab", sep="\t")
labeled_annotations = labels.merge(annotations, left_on="gene",right_on="orf_id")
genes = pd.read_csv("../satellite_models/get_gene_data/GeneCoverageMatrix_Normalized037_Binary.csv")

#pull genes from known genomes out
figs = [x for x in labeled_annotations['gene'] if 'fig' in x]
genomes = [x.split('.peg')[0] for x in figs]
genomeset = list(set(genomes))

genome_map = {'fig|1208365.4':'SAR86E',
              'fig|1123867.3':'SAR86B',
              'fig|1123866.3':'SAR86A',
              'fig|1123865.4':'SAR86D',
              'fig|1123864.3':'SAR86C'
             }

clusters = list(set(labeled_annotations['cluster']))

#find percent genes in each cluster belong to each of our five genomes
genestrains = pd.DataFrame({"cluster":clusters})
for genome in genomeset:
    col = []
    #subset genes corresponding to cluster
    for clust in clusters:
        subset = labels[labels['cluster']==clust]

        genes_in_genome = [x for x in subset['gene'] if genome in x]
        #percent genes in that cluster from that genome
        col.append(len(genes_in_genome)/len(subset) * 100)
    genestrains[genome] = col

#find whether genomes enriched in single clusters (SAR86A and E in particular of interest, other three have very few genes total)
genomeenrich = pd.DataFrame({"cluster":clusters})
for genome in genomeset:
    print("---",genome)
    total_genes_in_genome = [x for x in labels['gene'] if genome in x]
    expected_genome_in_cluster = len(total_genes_in_genome)/6.0
    print(len(total_genes_in_genome),"genes in genome")
    col = []
    #subset genes corresponding to cluster
    for clust in clusters:
        subset = labels[labels['cluster']==clust]

        genes_in_genome = [x for x in subset['gene'] if genome in x]
        #percent genes in that genome in that cluster
        col.append((len(genes_in_genome) - expected_genome_in_cluster)/expected_genome_in_cluster)
    genomeenrich[genome] = col

'''
#heatmaps for genestrains and genomeenrich

trace = go.Heatmap(z=[row[1:] for ix,row in genestrains.iterrows()],
                   x=[genome_map[c] for c in genestrains.columns[1:]],
                   y=genestrains['cluster'],
                   colorscale=[[0.0, 'rgb(255,255,255)'], [1.0, 'rgb(255,0,40)']],
                   zmin=0.2,
                   zmax=6
                  )
data=[trace]
py.offline.iplot(data, filename='cluster-heatmap')

trace = go.Heatmap(z=[row[1:] for ix,row in genomeenrich.iterrows()],
                   x=[genome_map[c] for c in genomeenrich.columns[1:]],
                   y=genomeenrich['cluster'],
                   colorscale=[[-1.0, 'rgb(40,0,255)'], [1.0, 'rgb(255,0,40)']],
                   zmin=-1,
                   zmax=1
                  )
data=[trace]
py.offline.iplot(data, filename='cluster-heatmap')
'''

## contig enrichment in clusters
scf = [x for x in labeled_annotations['gene'] if 'scf' in x]
contigs = [s.split('_')[0] for s in scf]

contigset = [x for x, count in Counter(contigs).most_common()]

contigenrich = pd.DataFrame({"cluster":clusters})
contig_counts = pd.DataFrame()
for contig in contigset:
    col = []
    counts = []

    total_genes_on_contig = [x for x in labels['gene'] if contig in x]
    expected_genes_per_cluster = len(total_genes_on_contig)/5.0

    overall = len(total_genes_on_contig)/len(labels)

    #subset genes corresponding to cluster
    for clust in clusters:
        subset = labels[labels['cluster']==clust]

        genes_in_genome = [x for x in subset['gene'] if contig in x]
        #percent genes in that cluster from that contig
        col.append((len(genes_in_genome) - expected_genes_per_cluster)/expected_genes_per_cluster)

        counts.append(len(genes_in_genome)/len(subset))
    counts.append(overall)

    contig_counts[contig] = counts
    contigenrich[contig] = col

contig_counts.index = ['c1','c2','c3','c4','c5','overall']


trace = go.Heatmap(z=[row[1:100] for ix,row in contigenrich.iterrows()],
                   x=contigenrich.columns[1:100],
                   y=contigenrich['cluster'],
                   colorscale=[[-1.0, 'rgb(40,0,255)'], [1.0, 'rgb(255,0,40)']],
                   zmin=-1,
                   zmax=1
                  )
data=[trace]
py.offline.iplot(data, filename='cluster-heatmap')


#piechart
contigbin = contigenrich>0
contigsums = Counter(contigbin.sum(axis=0)[1:])
clabels = [str(x[0])+" : "+str(x[1])+" contigs" for x in contigsums.items()]
ccounts = [count for label,count in contigsums.items()]

fig1, ax1 = plt.subplots()
ax1.pie(ccounts, labels=clabels, shadow=True, startangle=90, colors=["#6666FF","#66FFCC","#CC66FF","#FFFC79","#73FDFF"])
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.savefig("figures/contig_enriched_piechart.png")
plt.show()


## t test
# are contig enrichment in clusters significantly different than null? Answer: yes definitely
from scipy.stats import ttest_1samp

for cluster, row in contigenrich.iterrows():
    print(cluster+1)
    print(ttest_1samp(row, popmean=0))


## what's cluster membership of TARA sites?
print("analyzing cluster membership in TARA sites...")

genes.index = genes['index']
Tgenes = genes.T
Tgenes['gene'] = Tgenes.index
merged = labels.merge(Tgenes, on="gene")

cluster_membership = pd.DataFrame()
for gene in merged.columns[2:]:
    total_genes = np.sum(merged[gene])
    for clust in set(merged['cluster']):
        subset = merged[merged['cluster']==clust]
        count_percent = np.sum(subset[gene])/total_genes
        cluster_membership.loc[clust,gene] = count_percent
cluster_membership = cluster_membership.T

run_keys = pd.read_csv('../satellite_models/get_gene_data/runs_samples_keys.csv')
sampleinfo = pd.read_csv("../satellite_models/get_env_data/env_remote_data_TARA.csv")
env = run_keys.merge(sampleinfo,on='tara_label')

new = cluster_membership
cluster_membership['run_id'] = cluster_membership.index
new = env[['run_id','longitude']].merge(cluster_membership, on='run_id')
new = new.sort_values('longitude')

#TARA cluster membership sorted by longitude
sites = new['run_id']
lons = new['longitude']
colors = ["#6666FF", "#66FFCC", "#FFCC66", "#CC66FF", "#FFA5CA", "#07A6DE", "#FF6666"]
traces = []
for ix,clust in enumerate(new.columns[2:]):
    trace = go.Bar(
        x=sites,
        y=new[clust],
        name="c"+str(clust),
        marker=dict(color=colors[ix])
    )
    traces.append(trace)

data = traces
layout = go.Layout(
    barmode='stack'
)

fig = go.Figure(data=data, layout=layout)
py.offline.iplot(fig, filename='stacked-bar')


## does proportion c1/c5 and c3/c4 reflect true relative abundance SAR86E/A at TARA sites?
relative = pd.read_csv("relativeAbundanceSAR86.csv", index_col="species_id")
relative = relative.T

miniclusters = merged[merged['cluster'].isin([1,3,4,5])]
mini_membership = pd.DataFrame()
for site in miniclusters.columns[2:]:
    total_genes = np.sum(miniclusters[site])
    for clust in set(miniclusters['cluster']):
        subset = miniclusters[miniclusters['cluster']==clust]
        count_percent = np.sum(subset[site])/total_genes
        mini_membership.loc[clust,site] = count_percent

mini_membership = mini_membership.T

relative['site'] = relative.index
relative['SAR86Apercent'] = relative['SAR86A']/(relative['SAR86A']+relative['SAR86E']) * 100
relative['SAR86Epercent'] = relative['SAR86E']/(relative['SAR86A']+relative['SAR86E']) * 100

mini_membership['site'] = mini_membership.index
relmerge = relative.merge(mini_membership, on='site').dropna()
relmerge.head()

#scatterplot relative abundance SAR86A vs cluster proportion c3+c4 (as opposed to SAR86E only, so adds up to one and SAR86E is mirror of SAR86A, why only plot SAR86A below)
plt.figure(figsize=(8,8))
#plt.scatter(x=relmerge['SAR86Epercent'],y=relmerge[1]+relmerge[5], c='#FF0028', alpha=0.5)
plt.scatter(x=relmerge['SAR86Apercent'],y=relmerge[3]+relmerge[4], c='#2800FF', alpha=0.5)
plt.xlabel('relative abundance SAR86')
plt.ylabel('cluster proportion SAR86')
#plt.legend(labels=['SAR86E','SAR86A'])
plt.savefig("figures/genome_RelabunVsClusterproportion.png")
plt.show()


from scipy.stats.stats import pearsonr
print("relative abundance vs cluster proportion pearsonR:", pearsonr(x=relmerge['SAR86Apercent'],y=relmerge[3]+relmerge[4]))


# TARA station map

mapdata = pd.read_csv("mapinfo.csv")

colors=['green','yellow','red','turquoise','purple']
depths = ['DCM','SRF','MES.OMZ','MIX','MES']

traces = []
for ix, dep in enumerate(depths):
    #ix,subset in enumerate(mapdata.groupby('depth.id')):
    #ix = int(subset[0])
    data = pd.DataFrame(mapdata[mapdata['depth.id']==dep])
    traces.append(
            go.Scattergeo(
            lon = data['LONG'],
            lat = data['LAT'],
            name = list(data['depth.id'])[0],
            hovertext=data['sample_label'],
            mode = 'markers',
            marker = dict(
                size = 8,
                opacity = 0.4,
                symbol = 'circle',
                color = colors[ix])
            )
    )
data = go.Data(traces)

layout = go.Layout(
        title = 'TARA stations by depth',
        geo = dict(
            scope="world",
            showland = True,
            landcolor = "rgb(250, 250, 250)",
            subunitcolor = "rgb(217, 217, 217)",
            countrycolor = "rgb(217, 217, 217)",
            countrywidth = 0.5,
            subunitwidth = 0.5
        ),
    )

fig = go.Figure(data=data, layout=layout)
py.offline.iplot( fig, validate=False, filename='clustermap' )


# Function Enrichment in Clusters
print("doing functional enrichment in clusters...")
split = [str(x).split("||") for x in labeled_annotations["PFams"]]
id = []
for s in split:
    temp = []
    for br in s:
        if 'PF' in br:
            temp.append(br)
    if len(temp)>0:
        id.append(temp[0])
    else:
        id.append(np.nan)

labeled_annotations['parsed_PFams'] = id

sorted_pfams = [x for x, count in Counter(labeled_annotations['parsed_PFams']).most_common()][1:]

pfam_enriched = pd.DataFrame()
pfam_counts = pd.DataFrame()
for pfam in sorted_pfams:
    col = []
    counts = []

    #total percent genes annotated to pfam
    total_genes_in_pfam = labeled_annotations[labeled_annotations['parsed_PFams']==pfam]

    overall = len(total_genes_in_pfam)/len(labeled_annotations)

    #expected_pfam_in_clusters = len(total_genes_in_pfam)/5.0
    #pfam_enriched.loc[pfam,'overall'] = overall
    for clust in range(1,6):
        total_genes_in_clust = labeled_annotations[labeled_annotations['cluster']==clust]
        #expected number genes with this annotation in this cluster is proportional to total number of genes in cluster
        expected_pfam_in_clusters = (len(total_genes_in_clust)/len(labeled_annotations)) * len(total_genes_in_pfam)

        pfam_in_cluster = total_genes_in_clust[total_genes_in_clust['parsed_PFams']==pfam]
        col.append((len(pfam_in_cluster) - expected_pfam_in_clusters)/expected_pfam_in_clusters)

        counts.append(len(pfam_in_cluster)/len(total_genes_in_clust))

    counts.append(overall)

    pfam_enriched[pfam] = col
    pfam_counts[pfam] = counts

pfam_counts.index = ['c1','c2','c3','c4','c5','overall']

#pfam enriched heatmap
trace = go.Heatmap(z=[row[0:500] for ix,row in pfam_enriched.iterrows()],
                   x=pfam_enriched.columns[0:500],
                   y=['c1','c2','c3','c4','c5'],
                   colorscale=[[-1.0, 'rgb(40,0,255)'], [1.0, 'rgb(255,0,40)']],
                   zmin=-1,
                   zmax=1
                  )
data=[trace]
py.offline.iplot(data, filename='cluster-heatmap')

## t test
#
# are clusters' functional enrichments significantly diff than overall pangenome? Answer: no, except maybe c2 and maaaaybe c4
#
# note c2 is the cluster that predicts genes to be present at negative longitudes - aka Pacific basin. This is also the only locations where we have mesopelagic samples; given that we know functional capabilities stratify with depth, this seems the likely culprit to this one cluster being functionally different than the others; you also see that TARA sites with particularly high proportions of cluster 2 genes encompass most of the mesopelagic TARA samples, although this isn't uniformly true

from scipy.stats import ttest_1samp
for cluster, row in pfam_enriched.iterrows():
    print("c",cluster+1, "ttest")
    print(ttest_1samp(row, popmean=0))

#to see depths sampled in tara samples with highest proportions cluster 2
#new.sort_values(2, ascending=False).merge(env[['run_id','depth_sampled', 'tara_label']], on='run_id')
