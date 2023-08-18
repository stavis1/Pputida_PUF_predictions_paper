# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 15:13:10 2022

@author: Administrator
"""

import sys
import os
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.semi_supervised import SelfTrainingClassifier
from sortedcontainers import SortedList
import dill

data_file = sys.argv[1]
out_dir = os.path.abspath(sys.argv[2]) + '/'
if len(sys.argv) > 3:
    cores = int(sys.argv[3])
else:
    cores = os.cpu_count()

#####
n_trees = 100
depth = 60
self_train_threshold = 0.9
fdr_targ = 0.01
#####

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

data = pd.read_csv(data_file, sep = '\t')
labels = [l if p != 'Te' else -1 for l,p in zip(data['label'],data['partition'])]
features = data[[c for c in data.columns if c not in ['label','protein','term','partition']]]


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

test_true = [p for p,l,pa in zip(predictions[:,1],data['label'],data['partition']) if l == 1 and pa == 'Te']
test_false = [p for p,l,pa in zip(predictions[:,1],data['label'],data['partition']) if l == 0 and pa == 'Te']
tp = len([p for p in test_true if p > cutoff])
fp = len([p for p in test_false if p > cutoff])
tn = len([p for p in test_false if p < cutoff])
fn = len([p for p in test_true if p < cutoff])
print(f'confusion matrix: tp {tp}, fp {fp}, tn {tn}, fn {fn}')

out_data = pd.DataFrame({'protein':data['protein'],
                         'term':data['term'],
                         'partition':data['partition'],
                         'label':data['label'],
                         'prediction':predictions[:,1],
                         'predicted':[p > cutoff for p in predictions[:,1]]})
out_data.to_csv(f'{out_dir}predictions.tsv', sep = '\t', index = False)
with open(f'{out_dir}SLT01_termCentricSelfTrainedModel.dill','wb') as model_file:
    dill.dump(self_train_model, model_file)

