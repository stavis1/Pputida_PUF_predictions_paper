# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 17:02:20 2022

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
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--edgelist', type = str, required = True, help = 'training file')
parser.add_argument('--proteins', type = str, required = True, help = '.dill file of protein object list')
parser.add_argument('--terms', type = str, required = True, help = '.dill file of GO term object list')
parser.add_argument('-p','--prefix', type = str, required = True, help = 'output prefix, must include directory')
parser.add_argument('-r','--restart', 
                    type = str, 
                    required = False, 
                    help = 'load model from file instead of fitting', 
                    default = 'false')
args = parser.parse_args()
out_dir = os.path.abspath(args.prefix) + '/'

#########
depth = 30
n_trees = 140
use_weights = False
cores = 32
rng = np.random.default_rng(0)
sim_cut = 6
#########

#inputs
with open(args.proteins,'rb') as pickle:
    proteins = dill.load(pickle)
with open(args.terms,'rb') as pickle:
    terms = dill.load(pickle)
go_dict = {t.id:t for t in terms}
prot_dict = {p.pp_id:p for p in proteins}

test_set = set(rng.choice([p.pp_id for p in proteins if not p.unannotated], 300, replace = False))
with open(f'{out_dir}test_set.txt','w') as test_file:
    test_file.write('\n'.join(test_set))

unannotated = set([p.pp_id for p in proteins if p.unannotated])

def partition(p1,p2):
    if p1 in unannotated or p2 in unannotated:
        return 'Q'
    elif p1 in test_set or p2 in test_set:
        return 'Te'
    else:
        return 'Tr'

#data processing

dataset = pd.read_csv(args.edgelist, sep = '\t')
dataset = dataset.dropna()
dataset['partition'] = [partition(p1,p2) for p1,p2 in zip(dataset['protein1'],dataset['protein2'])]


train_dataset = dataset[[p == 'Tr' for p in dataset['partition']]]

train_dataset['depth_long'] = [d if d > 1 else 2 for d in train_dataset['depth_long']]
train_dataset['depth_long'] = [d if d < 10 else 10 for d in train_dataset['depth_long']]

train_dataset = train_dataset.append([train_dataset[train_dataset['depth_long'] > 6]] * 2, ignore_index=True)
train_dataset = train_dataset.append([train_dataset[train_dataset['depth_long'] > 8]] * 4, ignore_index=True)

pairs = train_dataset[['protein1','protein2']].copy()
train_dataset = train_dataset.drop(['protein1','protein2'], axis = 1)

train_features = train_dataset.copy()
train_features.pop('partition')
train_labels = train_features.pop('depth_long')

#Model fitting

title = 'randomForest_final_depth-{d}_trees-{n}_weights-{w}'.format(d = depth,
                                                                    n = n_trees,
                                                                    w = use_weights)

sys.stdout.write(title + ' started at ' + time.asctime(time.localtime(time.time())) + '\n')
sys.stdout.flush()

model = RandomForestRegressor(n_jobs = cores, max_depth=depth, n_estimators=n_trees, verbose=1)
model.fit(train_features, train_labels)

sys.stdout.write(title + ' model fit at ' + time.asctime(time.localtime(time.time())) + '\n')
sys.stdout.flush()

with open(args.prefix + title + '-model.dill', 'wb') as outfile:
    dill.dump(model, outfile)

#Predict 

dataset['predicted_depth'] = model.predict(dataset[[c for c in dataset.columns if c not in ['depth_long','protein1','protein2','partition']]])

dataset[['protein1', 'protein2', 'partition', 'depth_long', 'predicted_depth']].to_csv(f'{out_dir}SLT01_proteinSimilarityModelResults.tsv', 
                                                                                       sep = '\t',
                                                                                       index = False)

dataset[dataset['predicted_depth'] > sim_cut][[c for c in dataset.columns if c != 'depth_long']].to_csv(f'{out_dir}SLT01_filteredEdgelist.tsv',
                                                                                                        sep = '\t',
                                                                                                        index = False)

sys.stdout.write(title + ' finished at ' + time.asctime(time.localtime(time.time())) + '\n')
sys.stdout.flush()

print('this took: ' + str((time.time() - start)/60) + ' minutes\n')


