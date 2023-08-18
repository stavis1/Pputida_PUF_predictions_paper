# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 16:08:50 2022

@author: Administrator
"""
import Bio.PDB as bpdb
import json
import sys
import dill
import os
from multiprocessing import Pool

in_dir = os.path.abspath(sys.argv[1]) + '/'
out_dir = os.path.abspath(sys.argv[2]) + '/'
protein_dill = sys.argv[3]
if len(sys.argv) > 4:
    cores = int(sys.argv[4])
else:
    cores = os.cpu_count()

#######
min_score = 70
#######

def trim_file(structure_file, scores_json, out_file):
    class ResSelect(bpdb.Select):
        def accept_residue(self, res):
            if res.id[1] > good_start and res.id[1] < good_end:
                return True
            else:
                return False
    
    with open(scores_json,'rb') as json_file:
        scores = json.load(json_file)
    good_start = next(i for i,s in enumerate(scores['confidenceScore']) if s > min_score)
    good_end = len(scores['confidenceScore']) - next(i for i,s in enumerate(reversed(scores['confidenceScore'])) if s > min_score)
    
    if good_end - good_start < 30:
        return
    
    structure = bpdb.MMCIFParser().get_structure('structure', structure_file)
    io = bpdb.PDBIO()
    io.set_structure(structure)
    io.save(out_file, ResSelect())
    return


with open(protein_dill, 'rb') as dill_file:
    proteins = dill.load(dill_file)

prot_dict = {p.uni_id:p for p in proteins}
uni_pp = {p.uni_id:p.pp_id for p in proteins}
good_uni = set(uni_pp.keys())

cifs = [f for f in os.listdir(in_dir) if f.endswith('.cif')]
jsons = [f for f in os.listdir(in_dir) if f.endswith('.json')]
jobs = []
for uni in good_uni:
    try:
        cif = in_dir + next(c for c in cifs if uni in c)
    except StopIteration:
        continue
    scores_file = in_dir + next(j for j in jsons if uni in j)
    if prot_dict[uni].unannotated:
        out_file = out_dir + uni_pp[uni] + '-trim-U.pdb'
    else:
        out_file = out_dir + uni_pp[uni] + '-trim-A.pdb'
    jobs.append((cif,scores_file,out_file))

with Pool(cores) as p:
    p.starmap(trim_file, jobs)




