#run get_gene_data to get normalized gene presence/absence data for all the TARA sites
#outputs GeneCoverageMatrix_Normalized037_Binary.csv and GeneCoverageMatrix_Normalized.csv
cd satellite_models/get_gene_data
python get_gene_data.py

cd ..

#train logistic regression models for genes present at between 20-80% of TARA sites (and for all sites for comparison)
#outputs BalancedGenes_LogRegAccuracies_L1.csv and BalancedGenes_LogRegAccuracies_L1.csv
python run_logreg.py

cd ../clustering_env_predictiveness

#kmeans cluster model coefficients
python run_kmeans.py

#make map projections clusters
python make_cluster_maps.py

cd ../enrichment_analysis

#enrichment analysis
python enrichment_analysis.py
