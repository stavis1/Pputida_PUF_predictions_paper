# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 14:47:57 2022

@author: Administrator
"""

import sys
import os
import pandas as pd
import re
from sklearn.neighbors import KernelDensity
import numpy as np
# from sortedcontainers import SortedList
import subprocess

def get_peaks(vals):
    #peaks are defined as any value that is larger that two values on either side
    peaks = []
    for i in range(3,len(vals) - 3):
        if not vals[i] == 0:
            if all(j < vals[i] for j in [vals[i-2],vals[i-1],vals[i+1],vals[i+2]]):
                peaks.append(i)
    return peaks

def process_delta_mz(delta):
    if abs(delta) > 5:
        for peak in bins:
            if delta > peak[0] and delta < peak[2]:
                return delta - peak[1]
    return delta

def is_exp_file(experiment, file):
    if 'PXD013011' in file:
        return experiment in file
    else:
        return experiment in raw_new_dict[file.split('.')[0]]

def expect_mass(sequence):
    base_seq = re.sub(r'n?\[\d+\]','',sequence)
    mod_masses = [float(m) for m in re.findall(r'\[(\d+)\]',sequence)]
    return np.sum([aa_mass[aa] for aa in base_seq]) + np.sum(mod_masses) + 18.01

in_dir = os.path.abspath(sys.argv[1]) + '/'
out_dir = os.path.abspath(sys.argv[2]) + '/'
annsolo_dir = in_dir + 'SLT01_coEx-annSoLoOut/'
config_dir = in_dir + 'SLT01_annSoLo-config/'

aa_mass = {
    "A":71.03711, "C":103.00919, "D":115.02694, "E":129.04259, "F":147.06841, "G":57.02146, "H":137.05891, 
    "I":113.08406, "K":128.09496, "L":113.08406, "M":131.04049, "N":114.04293, "P":97.05276, "Q":128.05858,
    "R":156.10111, "S": 87.03203, "T":101.04768, "V": 99.06841, "W":186.07931, "Y":163.06333
}

with open(config_dir + 'valid_defects.txt','r') as defect_file:
    valid_defects = [float(d.strip()) for d in defect_file.readlines()]
with open(config_dir + 'names.txt','r') as name_file:
    lines = name_file.readlines()
raw_new_dict = {n.split('\t')[1].split('.')[0]:n.split('\t')[0].split('.')[0] for n in lines}
# new_raw_dict = {n.split('\t')[0].split('.')[0]:n.split('\t')[1].split('.')[0] for n in lines}

#information from the sepctral library file
lid_gene_dict = {} #library ID -> gene name
with open(config_dir + 'SLT01_coEx-concat.sptxt','r') as sptxt:
    for line in sptxt:
        if line.startswith('LibID:'):
            lib_id = re.search(r'LibID: (\d+)', line).group(1)
        elif line.startswith('Comment:'):
            gene = re.search(r'Protein=([^ ]+)',line).group(1)
            lid_gene_dict[lib_id] = gene

#identify modifications from mass defects
mztabs = [annsolo_dir + f for f in os.listdir(annsolo_dir) if f.endswith('.mztab')][:5]
deltas = []
for mztab in mztabs:
    results = pd.read_csv(mztab,sep ='\t',skiprows=(42))
    deltas.extend([(o_mz - e_mz) * z for e_mz,o_mz,z in zip(results['exp_mass_to_charge'],results['calc_mass_to_charge'],results['charge'])])

#kernel density model of mass deviations to find groups of values that look like real PTMs
big_deltas = np.array([d for d in deltas if abs(d) > 5])
x_d = np.linspace(-500, 500, 100000)
kde = KernelDensity(bandwidth=.3, kernel='gaussian')
kde.fit(big_deltas[:, None])
logprob = kde.score_samples(x_d[:, None])
density = np.exp(logprob)

#wide bandwith kernel density model of mass deviations for local noise floor estimation
kde_wide = KernelDensity(bandwidth=20, kernel='gaussian')
kde_wide.fit(big_deltas[:, None])
logprob_wide = kde_wide.score_samples(x_d[:, None])
baseline = np.exp(logprob_wide) + 0.001

#denoise 
density = [0 if d < b else d for d,b in zip(density,baseline)]

#get peaks and check that their mass defect corresponds with a valid mod
peaks = get_peaks(density)
good_peaks = [p for p in peaks if any(abs(p - m) <= 0.02 for m in valid_defects)]

#bin delta mass values around their peak
peaks.sort()
bins = []
for i,peak in enumerate(peaks):
    if peak in good_peaks:
        #if two adjacent peaks are not baseline resolved the breakpoint is set as the minimum value between them
        if i > 0 and min(density[peaks[i-1]:peak]) > 0:
            bin_lower = next(peaks[i-1] + j for j,v in enumerate(density[peaks[i-1]:peak]) if v == min(density[peaks[i-1]:peak]))
        #otherwise the edges of a peak are set to where it hits 0
        else:
            bin_lower = next(j for j in range(peak,0,-1) if density[j] == 0)
        if i < (len(peaks) - 1) and min(density[peak:peaks[i+1]]) > 0:
            bin_upper = next(peak + j for j,v in enumerate(density[peak:peaks[i+1]]) if v == min(density[peak:peaks[i+1]]))
            assert bin_upper != 42
        else:
            bin_upper = next(j for j in range(peak,len(density)) if density[j] == 0)
        bins.append((x_d[bin_lower],next(m for m in valid_defects if abs(peak - m) <= 0.02),x_d[bin_upper]))
assert all(b[0] < x_d[b[1]] and x_d[b[1]] < b[2] for b in bins)

#build experiment specific input file for flashLFQ
experiments = ['PXD013011','PXD016028','PXD016114']
for experiment in experiments:
    mztabs = [annsolo_dir + f for f in os.listdir(annsolo_dir) if f.endswith('.mztab') and is_exp_file(experiment, f)]
    data_out = out_dir + 'SLT01_coEx-' + experiment + '-flashLFQdata.tsv'
    with open(data_out,'w') as out_file:
        out_file.write('\t'.join(['File Name','Base Sequence','Full Sequence','Peptide Monoisotopic Mass',
                                  'Scan Retention Time','Precursor Charge','Protein Accession\n']))
    for mztab in mztabs:
        results = pd.read_csv(mztab,sep ='\t',skiprows=(42))
        results['retention_time'] = [r/60 for r in results['retention_time']]
        results['protein'] = [';'.join(re.findall(r'PP_\d{4}',lid_gene_dict[str(i)])) for i in results['opt_ms_run[1]_cv_MS:1003062_spectrum_index']]
        results['delta_mz'] = [o_mz - e_mz for e_mz,o_mz in zip(results['exp_mass_to_charge'],results['calc_mass_to_charge'])]
        results['delta_mass'] = [process_delta_mz(mz*z) for mz,z in zip(results['delta_mz'],results['charge'])]
        results['modified'] = [mz*z != d for mz,z,d in zip(results['delta_mz'],results['charge'],results['delta_mass'])]
        # results['calc_mass'] = [mz*z if not m else (mz*z) - d for mz,z,m,d in zip(results['calc_mass_to_charge'],results['charge'],
        #                                                                         results['modified'],results['delta_mass'])]
        results['expect_mass'] = [expect_mass(s) for s in results['sequence']]
        results['base_seq'] = [re.sub(r'n?\[\d+\]','',s) for s in results['sequence']]
        results['sequence'] = [s if not m else s + 'modified' for s,m in zip(results['sequence'],results['modified'])]
        if not experiment == 'PXD013011':
            results['filename'] = [raw_new_dict[mztab.split('/')[-1].split('.')[0]] + '.mzML'] * results.shape[0]
        else:
            results['filename'] = [mztab.split('/')[-1].split('.')[0] + '.mzML'] * results.shape[0]
        results = results.loc[[bool(p) for p in results['protein']]]
        results.to_csv(data_out, sep = '\t', mode = 'a', index = False, header = False, columns = ['filename','base_seq','sequence','expect_mass',
                                                                                        'retention_time','charge','protein'])
subprocess.run('sudo docker run --rm -v /home/cades/SLT02/data:/mnt/data smithchemwisc/flashlfq:1.0.3 --idt /mnt/data/SLT01_coEx-flashLFQ/SLT01_coEx-PXD013011-flashLFQdata.tsv --rep /mnt/data/mzml --out /mnt/data/PXD013011 --mbr --sha', 
               shell = True)

subprocess.run('sudo docker run --rm -v /home/cades/SLT02/data:/mnt/data smithchemwisc/flashlfq:1.0.3 --idt /mnt/data/SLT01_coEx-flashLFQ/SLT01_coEx-PXD016028-flashLFQdata.tsv --rep /mnt/data/mzml --out /mnt/data/PXD016028 --mbr --sha',
               shell = True)

subprocess.run('sudo docker run --rm -v /home/cades/SLT02/data:/mnt/data smithchemwisc/flashlfq:1.0.3 --idt /mnt/data/SLT01_coEx-flashLFQ/SLT01_coEx-PXD016114-flashLFQdata.tsv --rep /mnt/data/mzml --out /mnt/data/PXD016114 --mbr --sha',
               shell = True)

