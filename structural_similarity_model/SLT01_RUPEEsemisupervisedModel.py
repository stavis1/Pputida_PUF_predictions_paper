# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 12:42:24 2022

@author: Administrator
"""

import os
import sys
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.semi_supervised import SelfTrainingClassifier
from copy import deepcopy
import dill
from sortedcontainers import SortedList

data_file = os.path.abspath(sys.argv[1])
output_dir = os.path.abspath(sys.argv[2]) + '/'
if len(sys.argv) > 3:
    cores = int(sys.argv[3])
else:
    cores = os.cpu_count()

######
n_trees = 140
depth = 40
self_train_threshold = 0.95
fdr_targ = 0.01
######

def FDR_control(true_probs, false_probs):
    all_probs = sorted(list(true_probs) + list(false_probs))
    true_probs = SortedList(true_probs)
    false_probs = SortedList(false_probs)
    for cut in all_probs:
        fp = len(false_probs) - false_probs.bisect_right(cut)
        tp = len(true_probs) - true_probs.bisect_right(cut)
        fdr = fp/(fp+tp)
        if fdr < fdr_targ:
            return cut

input_data = pd.read_csv(data_file,sep = '\t')
features = input_data[['tm_max','tm_sum','rmsd_min','rmsd_mean',
                       'seqid_max','seqid_sum','notgap_max','notgap_sum',
                       'eukaryote_any','eukaryote_all','phylo_min','phylo_mean','n_hits']]
labels = deepcopy(input_data['label'])
labels = [l if p != 'Te' else -1 for l,p in zip(labels,input_data['partition'])]

base_model = RandomForestClassifier(n_estimators = n_trees,
                                    max_depth = depth,
                                    n_jobs = cores,
                                    random_state = 0,
                                    verbose = 1)

calibrated_model = CalibratedClassifierCV(base_model,
                                          n_jobs = cores)

self_train_model = SelfTrainingClassifier(calibrated_model,
                                          threshold = self_train_threshold,
                                          verbose = True)

self_train_model.fit(features, labels)
predictions = self_train_model.predict_proba(features)

true_probs = [p for p,l in zip(predictions[:,1],labels) if l == 1]
false_probs = [p for p,l in zip(predictions[:,1],labels) if l == 0]
cutoff = FDR_control(true_probs, false_probs)
print(f'cutoff: {cutoff}')

true_probs = [p for p,l,i in zip(predictions[:,1],input_data['label'],input_data['partition']) if l == 1 and i == 'Te']
false_probs = [p for p,l,i in zip(predictions[:,1],input_data['label'],input_data['partition']) if l == 0 and i == 'Te']
tp = len([p for p in true_probs if p > cutoff])
fp = len([p for p in false_probs if p > cutoff])
tn = len([p for p in false_probs if p < cutoff])
fn = len([p for p in true_probs if p < cutoff])

print(f'confusion matrix: tp {tp}, fp {fp}, tn {tn}, fn {fn}')

out_data = pd.DataFrame({'protein':input_data['protein'],
                         'term':input_data['term'],
                         'partition':input_data['partition'],
                         'label':input_data['label'],
                         'prediction':predictions[:,1],
                         'predicted':[p > cutoff for p in predictions[:,1]]})
out_data.to_csv(f'{output_dir}predictions.tsv', sep = '\t', index = False)

with open(f'{output_dir}SLT01_RUPEEselfTrainedModel.dill','wb') as model_file:
    dill.dump(self_train_model, model_file)






















