# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 15:11:50 2022

@author: Administrator
"""
import time
start = time.time()
import sys
import dill
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
from sortedcontainers import SortedList
# from matplotlib import cm
# import pickle as pkl
import argparse
from collections import Counter
from multiprocessing import Pool

parser = argparse.ArgumentParser()
parser.add_argument('--train', type = str, required = True, help = 'training file')
parser.add_argument('--test', type = str, required = True, help = 'test file')
parser.add_argument('--proteins', type = str, required = True, help = '.dill file of protein object list')
parser.add_argument('--terms', type = str, required = True, help = '.dill file of GO term object list')
parser.add_argument('-p','--prefix', type = str, required = True, help = 'output prefix, must include directory')
parser.add_argument('-r','--restart', 
                    type = str, 
                    required = False, 
                    help = 'load model from file instead of fitting', 
                    default = 'false')
args = parser.parse_args()

#########
depth = 30
n_trees = 140
use_weights = False
cores = 32
rng = np.random.default_rng(0)
seed = 0
fdr_targ = 0.1
#########

#[tp,fp,tn,fn]
def TPR(rates):
    tp,fp,tn,fn = rates
    divisor = tp + fn
    if divisor == 0:
        return 0
    else:
        return tp/divisor

def FPR(rates):
    tp,fp,tn,fn = rates
    divisor = fp + tn
    if divisor == 0:
        return 0
    else:
        return fp/divisor

def FDR(rates):
    tp,fp,tn,fn = rates
    divisor = tp + fp
    if divisor == 0:
        return 0
    else:
        return fp/divisor

#inputs
with open(args.proteins,'rb') as pickle:
    proteins = dill.load(pickle)
with open(args.terms,'rb') as pickle:
    terms = dill.load(pickle)
go_dict = {t.id:t for t in terms}
prot_dict = {p.pp_id:p for p in proteins}

#data processing

train_dataset = pd.read_csv(args.train, sep = '\t')
train_dataset = train_dataset.dropna()

train_dataset['depth_long'] = [d if d > 1 else 2 for d in train_dataset['depth_long']]
train_dataset['depth_long'] = [d if d < 10 else 10 for d in train_dataset['depth_long']]

train_dataset = train_dataset.append([train_dataset[train_dataset['depth_long'] > 6]] * 2, ignore_index=True)
train_dataset = train_dataset.append([train_dataset[train_dataset['depth_long'] > 8]] * 4, ignore_index=True)

pairs = train_dataset[['protein1','protein2']].copy()
train_dataset = train_dataset.drop(['protein1','protein2'], axis = 1)

train_features = train_dataset.copy()
train_labels = train_features.pop('depth_long')

#Model fitting

title = 'randomForest_final_depth-{d}_trees-{n}_weights-{w}'.format(d = depth,
                                                                    n = n_trees,
                                                                    w = use_weights)

if args.restart == 'false':
    sys.stdout.write(title + ' started at ' + time.asctime(time.localtime(time.time())) + '\n')
    sys.stdout.flush()
    
    model = RandomForestRegressor(n_jobs = cores, max_depth=depth, n_estimators=n_trees, verbose=1)
    model.fit(train_features, train_labels)
    
    sys.stdout.write(title + ' model fit at ' + time.asctime(time.localtime(time.time())) + '\n')
    sys.stdout.flush()
    
    with open(args.prefix + title + '-model.dill', 'wb') as outfile:
        dill.dump(model, outfile)
else:
    with open(args.restart,'rb') as dill_file:
        model = dill.load(dill_file)

#FDR control
test_prots = set(list(pairs['protein1']) + list(pairs['protein2']))
pseudo_test = set(rng.choice(list(test_prots),300))
pairs['predictions'] = model.predict(train_features)
pairs['depth_long'] = train_labels 
pair_hits = pairs[pairs.apply(lambda x: (x['protein1'] in pseudo_test) != (x['protein2'] in pseudo_test), axis = 1)]

def test_cut(cut):
    def FDR(rates):
        tp,fp,tn,fn = rates
        divisor = tp + fp
        if divisor == 0:
            return 0
        else:
            return fp/divisor
    
    hits = pair_hits[pair_hits['predictions'] > float(cut)].reset_index(drop = True)
    informative_queries = [p for p in set(list(hits['protein1']) + list(hits['protein2'])) if p in pseudo_test]
    true_vals = []
    false_vals = []
    for query in informative_queries:
        true_go = prot_dict[query].annotation
        term_hits = []
        query_hits = hits[hits.apply(lambda x: x['protein1'] == query or x['protein2'] == query, axis = 1)]
        for donor in query_hits.apply(lambda x: x['protein1'] if x['protein1'] != query else x['protein2'], axis = 1):
            term_hits.extend(prot_dict[donor].annotation)
        counts = Counter(term_hits)
        true_vals.extend([counts[t] for t in true_go])
        false_vals.extend([counts[t] for t in counts.keys() if t not in true_go])

    true_vals = SortedList(true_vals)
    false_vals = SortedList(false_vals)
    if not true_vals:
        return None
    cutoffs = range(-1,max(true_vals))
    confusions = []
    for cutoff in cutoffs:
        tp = len(true_vals) - true_vals.bisect_right(cutoff)
        fp = len(false_vals) - false_vals.bisect_right(cutoff)
        tn = false_vals.bisect_right(cutoff)
        fn = true_vals.bisect_right(cutoff)
        confusions.append([tp,fp,tn,fn])
    
    fdr_targ = 0.1
    min_delta = min([abs(FDR(c) - fdr_targ) for c in confusions])
    fdr_idx = next(i for i,c in enumerate(confusions) if abs(FDR(c) - fdr_targ) == min_delta)
    fdr_confusion = confusions[fdr_idx]
    fdr_cutoff = cutoffs[fdr_idx]
    return (fdr_confusion, fdr_cutoff, cut)

depth_cuts = np.linspace(5,7)
with Pool(cores) as p:
    results = [r for r in p.map(test_cut, depth_cuts) if r is not None]


best_result = max([r[0][0] for r in results])
fdr_cutoff = next(r[1] for r in results if r[0][0] == best_result)
depth_cutoff = next(r[2] for r in results if r[0][0] == best_result)
fdr_confusion = next(r[0] for r in results if r[0][0] == best_result)

  
#Training set analysis figures
binned = list(zip([round(p,1) for p in pairs['predictions']], pairs['depth_long']))
categories = (sorted(list(set([round(p,1) for p in pairs['predictions']]))), sorted(list(set(pairs['depth_long']))))
counts = Counter(binned)
heatmap = np.array([[counts[(i,j)] for j in categories[1]] for i in categories[0]])
heatmap = np.log(heatmap + 1)

colormap = plt.cm.get_cmap('plasma')
sm = plt.cm.ScalarMappable(cmap=colormap)
sm.set_clim(vmin=min(heatmap.flatten()), vmax=max(heatmap.flatten()))

fig, ax = plt.subplots()
im = ax.imshow(heatmap,cmap = 'plasma',  extent = [0, 10, 0, 10], origin = 'lower')
ax.set_yticks(categories[0])
ax.set_xticks(categories[1])
ax.set_yticklabels([str(c) for c in categories[0]])
ax.set_xticklabels([str(c) for c in categories[1]])
ax.set_ylabel('Predicted Depth')
ax.set_xlabel('Shared Term Depth')
clb = plt.colorbar(sm)
clb.ax.set_title('Log Counts',fontsize=8)
fig.savefig(args.prefix + 'trainDataHeatmap.png')
plt.close('all')

hits = pair_hits[pair_hits['predictions'] > float(fdr_cutoff)].reset_index(drop = True)
print(hits.shape)
informative_queries = [p for p in set(list(hits['protein1']) + list(hits['protein2'])) if p in pseudo_test]
print(len(informative_queries))
print(fdr_cutoff)
true_vals = []
false_vals = []
for query in informative_queries:
    true_go = prot_dict[query].annotation
    term_hits = []
    query_hits = hits[hits.apply(lambda x: x['protein1'] == query or x['protein2'] == query, axis = 1)]
    for donor in hits.apply(lambda x: x['protein1'] if x['protein1'] != query else x['protein2'], axis = 1):
        term_hits.extend(prot_dict[donor].annotation)
    counts = Counter(term_hits)
    true_vals.extend([counts[t] for t in true_go])
    false_vals.extend([counts[t] for t in counts.keys() if t not in true_go])

true_vals = SortedList(true_vals)
false_vals = SortedList(false_vals)
cutoffs = range(-1,max(true_vals))
confusions = []
for cutoff in cutoffs:
    tp = len(true_vals) - true_vals.bisect_right(cutoff)
    fp = len(false_vals) - false_vals.bisect_right(cutoff)
    tn = false_vals.bisect_right(cutoff)
    fn = true_vals.bisect_right(cutoff)
    confusions.append([tp,fp,tn,fn])


fig, ax = plt.subplots()
ax.plot([FPR(c) for c in confusions],[TPR(c) for c in confusions],'-k',linewidth = 1)
ax.plot([0,1],[0,1],'--r',linewidth = 0.5)
ax.set_aspect('equal')
ax.set_ylim((0,1))
ax.set_xlim((0,1))
ax.set_ylabel('True Positive Rate')
ax.set_xlabel('False Positive Rate')
ax.set_title('Training Set ROC')
ax.scatter(FPR(fdr_confusion),TPR(fdr_confusion), s= 3, c = 'r', marker = '.')
fig.savefig(args.prefix + 'trainDataROC.png')
plt.close('all')

fig, ax = plt.subplot()
ax.plot([r[2] for r in results], [r[0][0] for r in results], '-k', linewidth = 1)
ax.set_ylabel('Number of True Positives Identified')
ax.set_xlabel('Predicted Depth Cutoff')
fig.savefig(args.prefix + 'depthCutoffSensitivityProfile.png')
plt.close('all')

#Test Dataset
test_dataset = pd.read_csv(args.test, sep = '\t')
test_dataset = test_dataset.dropna()
pairs = test_dataset[['protein1','protein2']].copy()
test_dataset = test_dataset.drop(['protein1','protein2'], axis = 1)


test_dataset['depth_long'] = [d if d > 1 else 2 for d in test_dataset['depth_long']]
test_dataset['depth_long'] = [d if d < 10 else 10 for d in test_dataset['depth_long']]

pairs['depth_long'] = test_dataset.pop('depth_long')
pairs['predictions'] = model.predict(test_dataset)
pair_hits = pairs[pairs.apply(lambda x: (x['protein1'] in pseudo_test) != (x['protein2'] in pseudo_test), axis = 1)]

#Test set analysis figures
binned = list(zip([round(p,1) for p in pairs['predictions']], pairs['depth_long']))
categories = (sorted(list(set([round(p,1) for p in pairs['predictions']]))), sorted(list(set(pairs['depth_long']))))
counts = Counter(binned)
heatmap = np.array([[counts[(i,j)] for j in categories[1]] for i in categories[0]])
heatmap = np.log(heatmap + 1)

colormap = plt.cm.get_cmap('plasma')
sm = plt.cm.ScalarMappable(cmap=colormap)
sm.set_clim(vmin=min(heatmap.flatten()), vmax=max(heatmap.flatten()))

fig, ax = plt.subplots()
im = ax.imshow(heatmap,cmap = 'plasma')
ax.set_yticks(categories[0])
ax.set_xticks(categories[1])
ax.set_yticklabels([str(c) for c in categories[0]])
ax.set_xticklabels([str(c) for c in categories[1]])
ax.set_ylabel('Predicted Depth')
ax.set_xlabel('Shared Term Depth')
clb = plt.colorbar(sm)
clb.ax.set_title('Log Counts',fontsize=8)
fig.savefig(args.prefix + 'testDataHeatmap.png')
plt.close('all')

hits = pair_hits[pair_hits['predictions'] > float(fdr_cutoff)].reset_index(drop = True)
informative_queries = [p for p in set(list(hits['protein1']) + list(hits['protein2'])) if p in pseudo_test]
true_vals = []
false_vals = []
for query in informative_queries:
    true_go = prot_dict[query].annotation
    term_hits = []
    query_hits = hits[hits.apply(lambda x: x['protein1'] == query or x['protein2'] == query, axis = 1)]
    for donor in hits.apply(lambda x: x['protein1'] if x['protein1'] != query else x['protein2'], axis = 1):
        term_hits.extend(prot_dict[donor].annotation)
    counts = Counter(term_hits)
    true_vals.extend([counts[t] for t in true_go])
    false_vals.extend([counts[t] for t in counts.keys() if t not in true_go])

true_vals = SortedList(true_vals)
false_vals = SortedList(false_vals)
cutoffs = range(-1,max(true_vals))
confusions = []
for cutoff in cutoffs:
    tp = len(true_vals) - true_vals.bisect_right(cutoff)
    fp = len(false_vals) - false_vals.bisect_right(cutoff)
    tn = false_vals.bisect_right(cutoff)
    fn = true_vals.bisect_right(cutoff)
    confusions.append([tp,fp,tn,fn])


test_fdr_confusion = confusions[cutoffs.index(fdr_cutoff)]
fig, ax = plt.subplots()
ax.plot([FPR(c) for c in confusions],[TPR(c) for c in confusions],'-k',linewidth = 1)
ax.plot([0,1],[0,1],'--r',linewidth = 0.5)
ax.set_aspect('equal')
ax.set_ylim((0,1))
ax.set_xlim((0,1))
ax.set_ylabel('True Positive Rate')
ax.set_xlabel('False Positive Rate')
ax.set_title('Test Set ROC')
ax.scatter(FPR(confusions[cutoffs.index(fdr_cutoff)]),
           TPR(confusions[cutoffs.index(fdr_cutoff)]),
           s= 3, c = 'r', marker = '.')
ax.text(FPR(test_fdr_confusion) + .02,
        TPR(test_fdr_confusion) - .01,
        'Achieved FDR = {r}'.format(r=str(round(FDR(test_fdr_confusion)*100,1))),
        fontsize = 6)
fig.savefig(args.prefix + 'testDataROC.png')
plt.close('all')

sys.stdout.write(title + ' finished at ' + time.asctime(time.localtime(time.time())) + '\n')
sys.stdout.flush()

print('this took: ' + str((time.time() - start)/60) + ' minutes\n')


