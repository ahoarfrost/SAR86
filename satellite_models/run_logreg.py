# coding: utf-8
import pandas as pd
import numpy as np
#import plotly as py
#import plotly.graph_objs as go
#py.offline.init_notebook_mode(connected=True)

## read in gene data
print("reading in gene and env data...")
genes = pd.read_csv("get_gene_data/GeneCoverageMatrix_Normalized037_Binary.csv")
genes.rename(columns={'index': 'run_id'}, inplace=True)

#get list of variables going to predict
target_ids = genes.columns[1:]

#add tara sample id label
run_keys = pd.read_csv('get_gene_data/runs_samples_keys.csv')
genes = run_keys.merge(genes,on='run_id')

## get env feature data
env = pd.read_csv("get_env_data/env_remote_data_TARA_preprocessed.csv", index_col=False)

#merge run_keys with env so that each run id has the env data associated with it
env = run_keys.merge(env,on='tara_label')
env_features = ['latitude', 'longitude', 'depth_sampled', 'chl_satellite',
       'sst_satellite', 'par_satellite', 'pic_satellite', 'poc_satellite',
       'npp_satellite', 'calcite_annual_historical', 'chla_monthly_historical',
       'chla_annual_historical', 'chla_annualrange_historical',
       'cloudfraction_annual_historical', 'cloudfraction_monthly_historical',
       'cloudfraction_annualstdev_historical',
       'diffuseattenuation_annual_historical', 'daylength_monthly_historical',
       'dustflux_annual_historical', 'dustflux_monthly_historical',
       'dustflux_annualstdev_historical', 'solarinsolation_annual_historical',
       'solarinsolation_monthly_historical',
       'solarinsolation_annualstdev_historical',
       'distfromland_annual_historical', 'pycnoclinedepth_monthly_historical',
       'pycnoclinedepth_annualstdev_historical',
       'thermoclinedepth_monthly_historical',
       'thermoclinedepth_annualstdev_historical', 'nitrate_annual_historical',
       'nitrate_monthly_historical', 'npratio_monthly_historical',
       'oceandepth_historical', 'oxygendissolved_annual_historical',
       'oxygendissolved_monthly_historical',
       'oxygensaturation_annual_historical',
       'oxygensaturation_monthly_historical',
       'oxygenutilization_annual_historical',
       'oxygenutilization_monthly_historical', 'par_annual_historical',
       'ph_annual_historical', 'phosphate_annual_historical',
       'phosphate_monthly_historical', 'salinity_annual_historical',
       'salinity_monthly_historical', 'silicate_annual_historical',
       'silicate_monthly_historical', 'sst_annual_historical',
       'oceantemp_monthly_historical', 'oceantemp_annual_historical',
       'intercept']


## join gene and env data into one df

#join cnv_with_label with env to get features in right order
data = env.merge(genes,on=['run_id','tara_label'])

## remove sites with no/low abundance SAR86
print("removing sites with no/low abundance SAR86 or missing env feature data...")
coverage = pd.read_csv("get_gene_data/GeneCoverageMatrix_Normalized.csv", index_col='index')

five = [] #runs where >1000 genes have >5x coverage
ten = [] #runs where >500 genes have >10x coverage
none = [] #runs with <=1000 genes with >5x coverage
noneten = [] #runs with <=500 genes with >10x coverage
for ix, row in coverage.iterrows():
    #print min(row[target_ids]), max(row[target_ids]), row[target_ids].mean()
    if len(row[row>5])>1000:
        five.append(ix)
    else:
        none.append(ix)
    if len(row[row>10])>500:
        ten.append(ix)
    else:
        noneten.append(ix)

#going to go with the five array above; looking through none, those all appear to definitely be low-SAR86-abundance samples. 
#I am possibly including some samples that are high coverage of many of the genes but don't have the signature 
#peak at 5 or 10ish coverage we see where a particular strain is definitely present

data_short = data[data['run_id'].isin(five)]
data_short = data_short.dropna() #also drop any samples for which don't have env data

data_none = data[data['run_id'].isin(none)] #tara sites where sar86 not present/low abundance - these are mostly mesopelagic samples
print("...number of TARA sites used for model building:", len(data_short))

## split into train/test/validation sets
# 72% train, 8% valid, 20% test
print("splitting TARA sites into 72% training, 8% validation, and 20% test sets...")

from sklearn.model_selection import train_test_split
train, test = train_test_split(data_short, test_size = 0.2, random_state=4448)
train, valid = train_test_split(train, test_size=0.1, random_state=833)

## define targets and logreg fn (with L1)
balanced_4060_targets = []
balanced_2080_targets = []
nonzero_targets = []

for column in target_ids:
    yes_no = data_short[column]>0
    if sum(yes_no)>(len(data_short)*0.4) and sum(yes_no)<(len(data_short)*0.6): 
        balanced_4060_targets.append(column)
    if sum(yes_no)>(len(data_short)*0.2) and sum(yes_no)<(len(data_short)*0.8):  
        balanced_2080_targets.append(column)
    if sum(yes_no)>0: #gene found in at least one site
        nonzero_targets.append(column)
print("number of genes present at between 40-60% of TARA sites:",len(balanced_4060_targets))
print ("number of genes present at between 20-80% of TARA sites:", len(balanced_2080_targets))
print("number of genes found in at least one site:",len(nonzero_targets))


print("building models for",len(balanced_2080_targets),"genes with 20-80% balanced classes...")

def run_logregs(target_genes, feature_cols, train_data, test_data,penalty='l1',C=0.01,max_iter=500,verbose=0, random_state=1423):
    from sklearn.linear_model import LogisticRegression
    train_features = train_data[feature_cols]
    test_features = test_data[feature_cols]
    
    train_targets = train_data[target_ids]
    test_targets = test_data[target_ids]
    
    coefficients = pd.DataFrame()
    coefficients['feature'] = feature_cols

    accuracies = pd.DataFrame()
    accuracies['gene_id'] = target_genes
    accuracies['test_accuracy']=""
    accuracies['improvement_accuracy']=""
    accuracies['train_accuracy'] = ""
    accuracies['train_improvement_accuracy'] = ""
    accuracies.set_index('gene_id',drop=True,inplace=True)

    train_in = train_features.as_matrix()
    test_in = test_features.as_matrix()

    for target,gene_id in enumerate(target_genes):
        if verbose==1:
            if target <= 1 or (target <= 100 and target % 10 == 0) or (target <= 1000 and target % 100 == 0)             or (target <= 10000 and target % 1000 == 0) or target % 10000 == 0:
                print("processing gene_id %s: %s" % (target,gene_id))

        train_targets_in = np.array(train_targets[target_genes[target]])>0
        test_targets_in = np.array(test_targets[target_genes[target]])>0
        
        log = LogisticRegression(penalty=penalty, C=C, max_iter=max_iter, verbose=0)
        #if all one class in training set will throw error; just skip in this case
        try:
            log.fit(X=train_in,y=train_targets_in)
        except ValueError:
            print("problem with fit for target ", target)
            continue
        
        coef = log.coef_[0]
        test_accuracy = log.score(X=test_in,y=test_targets_in)
        train_accuracy = log.score(X=train_in,y=train_targets_in)
        
        #if gene is present in >50% of sites, majority class accuracy is positive_class_frequency; if below, it's 1-positive class (absence is majority class)
        positive_class_frequency = float(sum(test_targets_in))/len(test_targets_in)
        if positive_class_frequency>0.5:
            majority_class_accuracy = positive_class_frequency
        else:
            majority_class_accuracy = 1-positive_class_frequency

        improvement_accuracy = test_accuracy-majority_class_accuracy
        train_improvement_accuracy = train_accuracy-majority_class_accuracy

        coefficients[gene_id] = coef
        accuracies.loc[gene_id]['test_accuracy'] = test_accuracy
        accuracies.loc[gene_id]['improvement_accuracy'] = improvement_accuracy
        accuracies.loc[gene_id]['train_accuracy'] = train_accuracy
        accuracies.loc[gene_id]['train_improvement_accuracy'] = train_improvement_accuracy
        

    return coefficients, accuracies

'''
## tune penalty parameter
#C=0.7 looks like a reasonable balance of not too much regularization and not too much of a hit to validation accuracies, and was what was used here. Code for how this was decided below (or see accompanying notebook)

def hist_improvements(data_to_plot, names):
    traces = []
    for ix, data in enumerate(data_to_plot):
        trace = go.Histogram(
            x=data['improvement_accuracy'],
            opacity=0.5,
            name=names[ix]
        )
        traces.append(trace)


    data = traces
    layout = go.Layout()
    fig = go.Figure(data=data, layout=layout)

    py.offline.iplot(fig, filename='overlaid histogram')

Cs = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
accuracies = []
penalties = []
coefficients = []
for c in Cs:
    print("---training on penalty", c, "...---")
    coef, acc = run_logregs(target_genes=balanced_2080_targets, feature_cols=env_features, train_data=train, test_data=valid, 
                            penalty='l1',C=c,max_iter=100,verbose=1)
    accuracies.append(acc)
    coefficients.append(coef)
    penalties.append(c)

means = []
for c, df in enumerate(accuracies):
    means.append(np.mean(accuracies[c]['train_improvement_accuracy']))

traces = []
for c, df in enumerate(accuracies):  
    trace = go.Histogram(
        x = df['improvement_accuracy'],
        opacity=0.2,
        name=penalties[c]
    )
    traces.append(trace)

layout = go.Layout(barmode='overlay')
fig = go.Figure(data=traces, layout=layout)

py.offline.iplot(fig)

for c, df in enumerate(accuracies):
    print("-----",penalties[c],"------")
    print(df[['improvement_accuracy', 'train_improvement_accuracy']].sort_values('improvement_accuracy', ascending=False).head())

#how many env variables with nonzero coefficients?
for c,coef in enumerate(coefficients):
    print("------C = %s------" % penalties[c])
    for col in coef.columns[1:15]:
        print(sum(coef[col]!=0))

trace = go.Scatter(
    x = penalties,
    y = means
)

data = [trace]

py.offline.iplot(data)

### it looks like above C=0.5, the impact on accuracies isn't huge; let's test overfitting/underfitting at a few of the upper values
coef_06, acc_06 = run_logregs(target_genes=balanced_2080_targets, feature_cols=env_features, train_data=train, test_data=valid, penalty='l1',C=0.6,max_iter=200,verbose=0)
print("penalty c=0.6 percent overfit as expected:", len(acc_06[acc_06['test_accuracy']<acc_06['train_accuracy']])/float(len(acc_06)))
print("mean validation accuracy: ", acc_06['test_accuracy'].mean())

coef_07, acc_07 = run_logregs(target_genes=balanced_2080_targets, feature_cols=env_features, train_data=train, test_data=valid, penalty='l1',C=0.7,max_iter=200,verbose=0)
print("penalty c=0.7 percent overfit as expected: ", len(acc_07[acc_07['test_accuracy']<acc_07['train_accuracy']])/float(len(acc_07)))
print("mean validation accuracy: ", acc_07['test_accuracy'].mean())

coef_08, acc_08 = run_logregs(target_genes=balanced_2080_targets, feature_cols=env_features, train_data=train, test_data=valid, penalty='l1',C=0.8,max_iter=200,verbose=0)
print("penalty c=0.8 percent overfit as expected: ", len(acc_08[acc_08['test_accuracy']<acc_08['train_accuracy']])/float(len(acc_08)))
print("mean validation accuracy: ", acc_08['test_accuracy'].mean())

coef_09, acc_09 = run_logregs(target_genes=balanced_2080_targets, feature_cols=env_features, train_data=train, test_data=valid, penalty='l1',C=0.9,max_iter=200,verbose=0)
print("penalty c=0.9 percent overfit as expected: ", len(acc_09[acc_09['test_accuracy']<acc_09['train_accuracy']])/float(len(acc_09)))
print("mean validation accuracy: ", acc_09['test_accuracy'].mean())

coef_10, acc_10 = run_logregs(target_genes=balanced_2080_targets, feature_cols=env_features, train_data=train, test_data=valid, penalty='l1',C=1.0,max_iter=200,verbose=0)
print("penalty c=1.0 percent overfit as expected: ", len(acc_10[acc_10['test_accuracy']<acc_10['train_accuracy']])/float(len(acc_10)))
print("mean validation accuracy: ", acc_10['test_accuracy'].mean())
'''

'''
## Run multiple L1 models at C=0.7 with no random_state set, compare which features used most often and how consistent that is
#### (which features are nonzero are super consistent. see example code below.)

#define as function, take how big of rank genes want to take (top ten, top twenty), list of feature names, list of gene column names to look at, and dataframe of data
def plot_most_predictive(features, genes, data, penalty='l1', rank=None, plot=True):
    #create a dict that counts the number of times each environmental feature is in the top X most predictive genes
    counts = {}
    for feature in features:
        counts[feature] = 0
        
    num_feats = {}
    
    for gene in genes:
        data['sort'] = data[gene].abs()
        
        if penalty=='l1':
            #find nonzero features
            nonzero = data[data['sort'] > 0]
        
            if rank is not None:
                #find x most predictive coefficients for that gene
                most_predictive_feats = list(nonzero.sort_values('sort',ascending=False)['feature'][0:rank])

                #add one to counts for each feat in most_predictive_feats
                for feat in most_predictive_feats:
                    counts[feat] = counts[feat]+1    
            else:
                #add one to counts for each feat in nonzero 
                for feat in nonzero['feature']:
                    counts[feat] = counts[feat]+1
        
            #add one to num_feats for number of nonzero
            if len(data[data['sort'] > 0]) in num_feats.keys():
                num_feats[len(data[data['sort'] > 0])] = num_feats[len(data[data['sort'] > 0])] + 1
            else:
                num_feats[len(data[data['sort'] > 0])] = 1
                
        if penalty=='l2':
            #for each feature, add the absolute value of the weight to the count in counts
            for feat in counts.keys():
                counts[feat] = counts[feat] + float(data[data['feature']==feat]['sort'])
                
    if penalty=='l2':
        #divide each key in counts by number of genes to get average weight, rather than sum
        for key in counts.keys():
            counts[key] = counts[key]/float(len(genes))
    
    
    countsdf = pd.DataFrame.from_dict(counts,orient='index')
    countsdf.columns = ['count']
    countsdf_sorted = countsdf.sort_values('count',ascending=False)
    
    if penalty=='l1':
        featshist = pd.DataFrame.from_dict(num_feats,orient='index')
        featshist.columns = ['count']
        featshist_sorted = featshist.sort_values('count',ascending=False)

    #plot
    if plot==True:
        data = [go.Bar(
                    x=list(countsdf_sorted.index),
                    y=list(countsdf_sorted['count']),
                    marker=dict(
                    color='lavender'
                )
            )]

        layout = go.Layout(
                yaxis = dict(
                      domain = [0.2, 1]
                    ),
                title="Number of Models with a Feature as Nonzero"
            )

        fig = go.Figure(data=data, layout=layout)

        py.offline.iplot(fig)

    if penalty=='l1':
        return countsdf_sorted, featshist_sorted
    if penalty=='l2':
        return countsdf_sorted

temp_feats = sorted(env_features)
featcounts_07 = pd.DataFrame({'feature':temp_feats})
for x in range(10):
    print("processing run",x,"...")
    #train a logreg model, will have slightly different results each time
    coef, acc = run_logregs(target_genes=balanced_2080_targets, feature_cols=env_features, train_data=train, test_data=valid, 
                        penalty='l1',C=0.7,max_iter=500,verbose=0, random_state=None)

    countsdf, num_feats = plot_most_predictive(features=list(coef['feature']),genes=list(acc.index),data=coef, plot=False)
    
    if list(featcounts_07['feature'])==list(countsdf.sort_index().index):
        featcounts_07['count_'+str(x)] = list(countsdf.sort_index()['count'])
    else:
        raise ValueError("feature names not in same order")

featcounts_07.sort_values('count_0', ascending=False)
print(featcounts_07)
'''

'''
## compare L1 to L2 regularization
#### the mean improvement accuracy between l2 and l1 (at C=0.7) are similar, the amount of overfitting and the skew are similar; I'm going to go with L1 regularization at C=0.7, which has the added advantage of feature selection of the most important features

coef_l2, acc_l2 = run_logregs(target_genes=balanced_2080_targets, feature_cols=env_features, 
                              train_data=train, test_data=valid, 
                              penalty='l2',C=0.7, max_iter=1000,verbose=1)


print("percent models overfit as expected:"," len(acc_l2[acc_l2['test_accuracy']<acc_l2['train_accuracy']])/float(len(acc_l2)))

print("mean L2 improvement accuracy:", acc_l2['improvement_accuracy'].mean(), "median L2 improvement accuracy:", acc_l2['improvement_accuracy'].median())
print("mean L1 improvement accuracy:", acc_07['improvement_accuracy'].mean(), "median L1 improvement accuracy:", acc_07['improvement_accuracy'].median())

countsdf_l2, num_feats_l2 = plot_most_predictive(features=list(coef_l2['feature']),genes=list(acc_l2.index),
                                                 data=coef_l2, plot=False)

counts_merged = pd.DataFrame({'l1':list(countsdf['count']),
                              'l1_feat':countsdf.index,
                              'l2':list(countsdf_l2['count']),
                             'l2_feat':countsdf_l2.index})
print ("percent of top-5 l1 feats also in l2 feats:")
print (len([x for x in counts_merged['l1_feat'][0:5] if x in list(counts_merged['l2_feat'][0:5])])/5.0)
print ("percent of top-10 l1 feats also in l2 feats:")
print (len([x for x in counts_merged['l1_feat'][0:10] if x in list(counts_merged['l2_feat'][0:10])])/10.0)
print ("percent of top-20 l1 feats also in l2 feats:")
print (len([x for x in counts_merged['l1_feat'][0:20] if x in list(counts_merged['l2_feat'][0:20])])/20.0)
print ("percent of top-30 l1 feats also in l2 feats:")
print (len([x for x in counts_merged['l1_feat'][0:30] if x in list(counts_merged['l2_feat'][0:30])])/30.0)
print ("percent of top-40 l1 feats also in l2 feats:")
print (len([x for x in counts_merged['l1_feat'][0:40] if x in list(counts_merged['l2_feat'][0:40])])/40.0)
print ("percent of top-50 l1 feats also in l2 feats:")
print (len([x for x in counts_merged['l1_feat'][0:50] if x in list(counts_merged['l2_feat'][0:50])])/50.0)
'''

# So let's run one final time (using test data instead of validation set for test accuracy this time), and save the final model

coef, acc = run_logregs(target_genes=balanced_2080_targets, feature_cols=env_features, 
                        train_data=train, test_data=test, 
                        penalty='l1',C=0.7,max_iter=1500,verbose=0)

acc_sorted = acc.sort_values('improvement_accuracy', ascending=False)
print("mean test accuracy: ", acc_sorted['test_accuracy'].mean())
print("median test accuracy: ", acc_sorted['test_accuracy'].median())
print("std dev test accuracy: ", acc_sorted['test_accuracy'].std())

#save the accuracies and coefficient info
print("saving results...")
coef.to_csv('Balanced2080Genes_LogRegCoefficients_L1.csv', index=False)
acc_sorted.to_csv('Balanced2080Genes_LogRegAccuracies_L1.csv')

'''
# making plots, take a look at characteristics of the model
countsdf, num_feats = plot_most_predictive(features=list(coef['feature']),genes=list(acc.index),data=coef)

#histogram train vs test accuracies
trace0 = go.Histogram(
    x=acc['test_accuracy'],
    opacity=0.5,
    xbins=dict(start=acc['test_accuracy'].min(), size=0.035, end=acc['test_accuracy'].max()),
    name="test accuracy"
)
trace1 = go.Histogram(
    x=acc['train_accuracy'],
    opacity=0.5,
    xbins=dict(start=acc['train_accuracy'].min(), size=0.035, end=acc['train_accuracy'].max()),
    name="train accuracy"
)

data = [trace0, trace1]

layout = go.Layout(barmode='overlay')

fig = go.Figure(data=data, layout=layout)
py.offline.iplot(fig)

#histogram test improvement accuracies vs train improvement accuracies
trace0 = go.Histogram(
    x=acc['improvement_accuracy'],
    opacity=0.5,
    xbins=dict(start=acc['improvement_accuracy'].min(), size=0.035, end=acc['improvement_accuracy'].max()),
    name="improvement accuracy"
)
trace1 = go.Histogram(
    x=acc['train_improvement_accuracy'],
    opacity=0.5,
    xbins=dict(start=acc['train_improvement_accuracy'].min(), size=0.035, end=acc['train_improvement_accuracy'].max()),
    name="train improvement accuracy"
)

data = [trace0, trace1]

layout = go.Layout(barmode='overlay')
fig = go.Figure(data=data, layout=layout)
py.offline.iplot(fig)

#histogram number of nonzero features per model
data = [go.Bar(
                x=list(num_feats.index),
                y=list(num_feats['count']),
                marker=dict(
                color='pink'
            )
        )]

layout = go.Layout(
        title="Number of Nonzero Features per Model"
    )
fig = go.Figure(data=data, layout=layout)

py.offline.iplot(fig)

#improvement accuracy plot
genes = acc_sorted.index
x_vals = np.arange(len(genes))
improvements = acc_sorted['improvement_accuracy']

trace0 = go.Bar(
    x = x_vals,
    y = improvements,
    text= genes,
    marker=dict(
        color='rgb(130,180,200)',
        line=dict(
            color='rgb(8,48,107)'
        )
    ),
    opacity=0.9
)

plot_data = [trace0]
layout = go.Layout(
    title='Improvement accuracy for all models (balanced genes only)',
)

fig = go.Figure(data=plot_data, layout=layout)
py.offline.iplot(fig)

#heatmap of coefficients for each env feature for each gene model
countsdf['order'] = range(0,len(countsdf))
countsdf['feature'] = countsdf.index
acc['gene_id'] = acc.index
coef_heatmap = countsdf[['feature','order']].merge(coef,on="feature")
coef_heatmap = coef_heatmap[(coef_heatmap.T != 0).any()]

trace = go.Heatmap(z=[list(row[acc['gene_id']]) for idx, row in coef_heatmap.iterrows()],
                   x=acc['gene_id'],
                   y=coef_heatmap['feature'],
                   zmin=-2,
                   zmax=2)

layout = go.Layout(
            xaxis = dict(
                  domain = [0.1, 1]
                ),
            yaxis = dict(
                  domain = [0.2, 1]
                )
        )
    
fig = go.Figure(data=[trace], layout=layout)

py.offline.iplot(fig)


# This is really interesting - the features are ordered, from bottom to top, by how many models contained it as a nonzero feature. It really looks like there are some clusters of genes that were predicted by different sets of feats than others. This should be interesting to look into more with the annotations. 
'''

# Run logreg on *all* the genes, even those with less balanced classes
print("training models on *all* the genes, including those with very imbalanced classes...")
print("(note there are a few genes only present at one site in training and none in test set, skipping these genes)")
#there are several genes in nonzero_targets that are absent from the training data (e.g. only one present and it got shuffled into test data)
#so remove those and use the new set
train_targets = train[nonzero_targets]
nonzero = [x for ix, x in enumerate(nonzero_targets) if sum(np.array(train_targets[nonzero_targets[ix]])>0)]

coef_all, acc_all = run_logregs(target_genes=nonzero, feature_cols=env_features, train_data=train, test_data=test, penalty='l1',C=0.7,max_iter=1000,verbose=0)

acc_all.replace('', np.nan, inplace=True) #replace empty strings in those rows with problem with fit with np.nan so can dropna
acc_all_sorted = acc_all.dropna().sort_values('improvement_accuracy', ascending=False)
print("mean test accuracy: ", acc_all_sorted['test_accuracy'][5:].mean())
print("median test accuracy: ", acc_all_sorted['test_accuracy'][5:].median())
print("std dev test accuracy: ", acc_all_sorted['test_accuracy'][5:].std())
acc_all_sorted = acc_all_sorted[5:]

for gene in acc_all_sorted.index:
    yes_no = data_short[gene]>0
    proportion = sum(yes_no)/float(len(data_short))
    acc_all_sorted.loc[gene,"proportion_sites"] = proportion
acc_all_sorted.tail()

#save the accuracies and coefficient info
print("saving results...")
coef_all.to_csv('AllGenes_LogRegCoefficients_L1.csv')
acc_all_sorted.to_csv('AllGenes_LogRegAccuracies_L1.csv')


'''
#some plots for the all-genes models

#improvement accuracy vs proportion sites gene is present at (worse accuracies for imbalanced genes)
trace = go.Scatter(
    x = acc_all_sorted['proportion_sites'],
    y = acc_all_sorted['improvement_accuracy'],
    mode = 'markers',
    marker= dict(
                opacity= 0.05
                )
)

data = [trace]

py.offline.iplot(data)
### so there is some effect on accuracy with the more imbalanced class genes, as expected

#improvement accuracy plot
genes = acc_all_sorted.index
x_vals = np.arange(len(genes))
improvements = acc_all_sorted['improvement_accuracy']

trace0 = go.Bar(
    x = x_vals,
    y = improvements,
    text= genes,
    marker=dict(
        color='rgb(158,202,225)',
        line=dict(
            color='rgb(8,48,107)'
        )
    ),
    opacity=0.6
)

plot_data = [trace0]
layout = go.Layout(
    title='Improvement accuracy for all models (balanced genes only)',
)

fig = go.Figure(data=plot_data, layout=layout)
py.offline.iplot(fig)

#num models env features most predictive in 
countsdf_all, num_feats_all = plot_most_predictive(features=list(coef_all['feature']),genes=list(acc_all_sorted.index),data=coef_all)

countsdf_all['feature'] = countsdf_all.index
countsdf_merged = countsdf_all.merge(countsdf[['count','feature', 'order']], on="feature")

print("percent of top-5 all-gene feats also in balanced-gene feats:")
print(len([x for x in countsdf['feature'][0:5] if x in list(countsdf_all.index[0:5])])/5.0)
print("percent of top-10 all-gene feats also in balanced-gene feats:")
print(len([x for x in countsdf['feature'][0:10] if x in list(countsdf_all.index[0:10])])/10.0)
print("percent of top-20 all-gene feats also in balanced-gene feats:")
print(len([x for x in countsdf['feature'][0:20] if x in list(countsdf_all.index[0:20])])/20.0)
print("percent of top-30 all-gene feats also in balanced-gene feats:")
print(len([x for x in countsdf['feature'][0:30] if x in list(countsdf_all.index[0:30])])/30.0)
'''

print("...done!")