# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 14:21:01 2022

@author: Administrator
"""

import subprocess
import os
import json
import numpy as np

json_files = [f for f in os.listdir() if f.endswith('.json.gz')]

with open('SLT01_alphafoldPredsSummaryMetrics.txt','w') as out_file:
    for j_file in json_files:
        subprocess.run(f'gunzip {j_file}', shell = True)
        with open(j_file[:-3],'rb') as json_file:
            scores = json.load(json_file)
        mean = np.mean(scores['confidenceScore'])
        subprocess.run(f'rm {j_file[:-3]}', shell = True)
        out_file.write(f'{j_file}\t{mean}\n')
    

