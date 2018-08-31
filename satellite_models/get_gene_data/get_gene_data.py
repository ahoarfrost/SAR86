# # get gene coverage for all the TARA runs
import pandas as pd

#list run ids
run_keys = pd.read_csv("runs_samples_keys.csv")

coverage = pd.DataFrame()

for run_id in run_keys['run_id']:
    print("processing run", run_id, "...")
    #read in genes file
    genes_path = 'sar86_genes/'
    file_name = run_id+'.SAR86.genes'
    try:
        cnv = pd.read_csv(genes_path+file_name,delimiter='\t')
        #drop copy_number column 
        cnv = cnv.drop('copy_number', axis=1)
        #pivot dataframe
        cnv['index'] = run_id 
        pivoted = cnv.pivot(index='index',columns='gene_id',values='coverage')
        #append to copy_number dataframe
        if len(coverage)==0:
            coverage = pivoted
        else:
            coverage = coverage.append(pivoted)
    except FileNotFoundError:
        print("file",run_id,";",run_keys[run_keys['run_id']==run_id]['tara_label'].to_string(index=False),"not found, skipping")

print("saving coverage matrix...")
coverage.to_csv('GeneCoverageMatrix.csv',index_label="run_id")


# # normalize coverage to scale to 1 for each TARA site 
# (divide each row by the row's mean)

def scale_rows_to_one(df):
    #mean of each row if axis is 1, mean of each column if axis is 0
    rowmeans = df.mean(axis=1, skipna=True)
    #divide each row by its mean
    norm = df.div(rowmeans, axis=0)
    return norm

print("normalizing gene coverage data...")
normalized = scale_rows_to_one(df=coverage)
print("saving normalized coverage matrix...")
normalized.to_csv("GeneCoverageMatrix_Normalized.csv")


'''
histograms were inspected manually to choose the appropriate presence/absence gene coverage threshold. See below for examples of how to do this:
# ### make some histograms of per-sample gene depth

import matplotlib.pyplot as plt

plt.hist(normalized.loc['ERR598993',:], bins=50, range=(1.5,20))
plt.hist(normalized.loc['ERR598993',:], bins=50, range=(1,20))
plt.hist(normalized.loc['ERR315859',:], bins=50, range=(2,20))
plt.hist(normalized.loc['ERR315859',:], bins=50, range=(1,20))
plt.hist(normalized.loc['ERR315858',:], bins=50, range=(2,20))
plt.hist(normalized.loc['ERR315861',:], bins=50, range=(1,20))
plt.hist(normalized.loc['ERR598951',:], bins=50, range=(1,20))
plt.hist(normalized.loc['ERR599043',:], bins=50, range=(0.6,20))
plt.hist(normalized.loc['ERR318620',:], bins=50, range=(1,20))
plt.hist(normalized.loc['ERR318620',:], bins=50, range=(0.37,20))
plt.hist(normalized.loc['ERR315862',:], bins=50, range=(0.37,20))
plt.hist(normalized.loc['ERR315862',:], bins=50, range=(1,20))

#plot normalized histograms, see a peak below e^-1 that makes it easier to see what looks like noise and actual low coverage gene data

import numpy as np
plt.hist(np.log(normalized.loc['ERR315862',:]+0.00000000000001), bins=50, range=(-10,10))
plt.hist(np.log(normalized.loc['ERR315862',:]+0.00000000000001), bins=50, range=(-10,10))
plt.hist(np.log(normalized.loc['ERR599043',:]+0.00000000000001), bins=50, range=(-5,5))

# It looks to me like particularly where there's a really clear modal/bimodal/(trimodal?) distribution, a cutoff of e^-1 looks pretty good, like it more or less bisects the peak in the very small numbers (e^-3, e^-4...) and the definitely-signal peak in the higher numbers. But, sometimes it looks like it's closer to e^0, or could be a smaller number...
# 
# I'm going to go with a cutoff of e^-1 (=0.37).
'''

# ## convert to presence/absence
# 
# use 0.37 as cutoff (which is e^-1)

def binarize_matrix(df, threshold=0.37):
    new_df = df.copy()
    new_df[new_df > threshold] = 1
    new_df[new_df <= threshold] = 0
    return new_df

normalized_37 = binarize_matrix(normalized, threshold=0.37)
normalized_37.to_csv("GeneCoverageMatrix_Normalized037_Binary.csv")

