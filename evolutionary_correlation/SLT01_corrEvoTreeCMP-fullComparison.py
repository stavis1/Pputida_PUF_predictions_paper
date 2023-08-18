# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 12:05:37 2021

@author: Administrator
"""
import sys
import os
import subprocess
from multiprocessing import cpu_count
import itertools
import time
import math

def treecmp_within(infile):
    outfile = infile[0] + '.w.out'
    sys.stdout.write(r'java -jar ./TreeCmp_v2.0-b76/bin/treeCmp.jar' +
                                ' -m ' + #matrix mode
                                ' -d ' + metric + #what metrics to use
                                ' -i ' + tmpdir + infile[0] + #input file has the reference newick trees
                                ' -o ' + outdir + outfile + 
                                ' -P -W')
    sys.stdout.flush()
    comparison = subprocess.Popen(r'java -jar ./TreeCmp_v2.0-b76/bin/treeCmp.jar' +
                                ' -m ' + #matrix mode
                                ' -d ' + metric + #what metrics to use
                                ' -i ' + tmpdir + infile[0] + #input file has the reference newick trees
                                ' -o ' + outdir + outfile + 
                                ' -P -W',#P prunes unpaired species and W allows for 0 distance branches
                                shell = True)
    return (comparison, outfile)

def treecmp_between(files):
    outfile = '__'.join(files) + '.b.out'
    sys.stdout.write(r'java -jar ./TreeCmp_v2.0-b76/bin/treeCmp.jar' +
                                ' -r ' + tmpdir + files[0] + #reference mode makes all pairwise comparisons between files 1 and 2
                                ' -d ' + metric +
                                ' -i ' + tmpdir + files[1] +
                                ' -o ' + outdir + outfile +
                                ' -P -W')
    sys.stdout.flush()
    comparison = subprocess.Popen(r'java -jar ./TreeCmp_v2.0-b76/bin/treeCmp.jar' +
                                ' -r ' + tmpdir + files[0] + #reference mode makes all pairwise comparisons between files 1 and 2
                                ' -d ' + metric +
                                ' -i ' + tmpdir + files[1] +
                                ' -o ' + outdir + outfile +
                                ' -P -W',
                                shell = True)
    return (comparison, outfile)

def collect_newicks(block):
    block_name = block[0] + '-' + block[-1] + '.newick'
    with open(tmpdir + block_name,'w') as outfile:
        for gene in block:
            with open(indir + gene + '/' + gene + '.raxml.bestTree', 'r') as infile:
                outfile.write(infile.read())
    return block_name

def divide_genes(genes, cores):
    step = math.ceil(len(genes)/cores)
    for i in range(cores):
        yield genes[i*step:(i*step)+step]

sys.stdout.write('START ' + time.asctime(time.localtime(time.time())) + '\n')
sys.stdout.flush()
indir = sys.argv[1]#this should have the RAxML-ng output folders in it named as /<locus code>_faa/
outdir = sys.argv[2]
# metric = sys.argv[3]#metric for TreeCMP to use, in the abbreviated code from TreeCMP
# cutoff = sys.argv[3]
if len(sys.argv) > 3:
    cores = int(sys.argv[3])
else:
    cores = cpu_count()

metric = 'rcw gdu'
genes = os.listdir(indir)
genes.sort()
# gene_blocks = [genes[i*:(i*cores)+cores] for i in range(cores)]
gene_blocks = list(divide_genes(genes, cores))

tmpdir = outdir + 'tmp/'
os.mkdir(tmpdir)

block_names = []
name_block_dict = {}
for block in gene_blocks:
    block_names.append(collect_newicks(block))
    name_block_dict[block_names[-1]] = block

queue = [tuple([n]) for n in block_names]
queue.extend(itertools.combinations(block_names, 2))
sys.stdout.write(str(queue))
sys.stdout.flush()


running = []
q_outfile_dict = {}
while queue:
    for proc in running:
        if proc.poll() is not None:
            sys.stdout.write(str(proc.pid) + ' finished '  + time.asctime(time.localtime(time.time())) + '\n')
            sys.stdout.flush()
            running.remove(proc)
    if len(running) < cores:
        if len(queue[-1]) == 2:
            cmp = treecmp_between(queue[-1])
            running.append(cmp[0])
            q_outfile_dict[queue[-1]] = cmp[1]
            sys.stdout.write(str(cmp[0].pid) + ' ' + cmp[1] + ' started ' + time.asctime(time.localtime(time.time())) + '\n')
            sys.stdout.flush()
        elif len(queue[-1]) == 1:
            cmp = treecmp_within(queue[-1])
            running.append(cmp[0])
            q_outfile_dict[queue[-1]] = cmp[1]
            sys.stdout.write(str(cmp[0].pid) + ' ' + cmp[1] + ' started ' + time.asctime(time.localtime(time.time())) + '\n')
            sys.stdout.flush()
        else:
            sys.stdout.write('there has been an error in the queue entry')
            sys.stdout.write(str(queue))
            sys.stdout.flush()
        queue.pop(-1)
    else:
        time.sleep(10)
for proc in running:
    proc.wait()

flag = True
with open(outdir + 'SLT01_corrEvo-rawRCW.tsv', 'w') as outfile:
    for q in q_outfile_dict.keys():
        filename = outdir + q_outfile_dict[q]
        with open(filename, 'r') as cmp_file:
            if flag:
                outfile.write(cmp_file.readline())
                flag = False
            else:
                cmp_file.readline()
            if filename[-5] == 'w':
                block = name_block_dict[q_outfile_dict[q][:-6]]
                sys.stdout.write(str(block) + '\n')
                sys.stdout.flush()
                for line in cmp_file:
                    fields = line.strip().split('\t')
                    fields[1] = block[int(fields[1]) - 1]
                    fields[2] = block[int(fields[2]) - 1]
                    outfile.write('\t'.join(fields) + '\n')
            elif filename[-5] == 'b':
                block1 = name_block_dict[q_outfile_dict[q].split('__')[0]]
                block2 = name_block_dict[q_outfile_dict[q].split('__')[1][:-6]]
                sys.stdout.write(str(block1) + '\n')
                sys.stdout.write(str(block2) + '\n')
                sys.stdout.flush()
                for line in cmp_file:
                    try:
                        fields = line.strip().split('\t')
                        fields[1] = block1[int(fields[1]) - 1]
                        fields[2] = block2[int(fields[2]) - 1]
                        outfile.write('\t'.join(fields) + '\n')
                    except IndexError:
                        sys.stdout.write(str(fields))
                        sys.stdout.flush()
            else:
                sys.stdout.write('theres some kind of name error in the final write step')
                sys.stdout.flush()
sys.stdout.write('END ' + time.asctime(time.localtime(time.time())) + '\n')
sys.stdout.flush()





    










