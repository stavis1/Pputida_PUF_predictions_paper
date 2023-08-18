# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 18:45:45 2022

@author: Administrator
"""
import sys
import os
import subprocess
import time
import itertools


def job(data, mode):
    if mode == 'between':
        command = ' '.join(['TMalign',
                            '-dir1 ' + tmp_dir + str(data[0]) + '/',
                            tmp_dir + 'manifest' + str(data[0]) + '.txt',
                            '-dir2 ' + tmp_dir + str(data[1]) + '/',
                            tmp_dir + 'manifest' + str(data[1]) + '.txt',
                            '-outfmt 2',
                            '> ' + out_dir + str(data[0]) + '_' + str(data[1]) + '.out'])
    elif mode == 'within':
        command = ' '.join(['TMalign',
                            '-dir ' + tmp_dir + str(data) + '/',
                            tmp_dir + 'manifest' + str(data) + '.txt',
                            '-outfmt 2',
                            '> ' + out_dir + str(data) + '.out'])
    else:
        raise ValueError
    sys.stdout.write('command: ' + command + ' started at ' + time.asctime(time.localtime(time.time())) + '\n')
    sys.stdout.flush()
    return subprocess.Popen(command, shell = True)

    
in_dir = os.path.abspath(sys.argv[1]) + '/'
tmp_dir = os.path.abspath(os.getcwd()) + '/' + sys.argv[2]
out_dir = os.path.abspath(sys.argv[3])  + '/'
if len(sys.argv) > 4:
    cores = int(sys.argv[4])
else:
    cores = os.cpu_count()

os.mkdir(tmp_dir)
os.system('rm {d}*'.format(d = out_dir))

cifs = [f for f in os.listdir(in_dir) if f.endswith('.pdb')]
subsets = [[] for _ in range(cores)]
for i,cif in enumerate(cifs):
    subsets[i%os.cpu_count()].append(cif)

print(subsets)

for i,subset in enumerate(subsets):
    subset_dir = tmp_dir + str(i) + '/'
    os.mkdir(subset_dir)
    for cif in subset:
        os.symlink(in_dir + cif, subset_dir + cif)
    with open(tmp_dir + 'manifest' + str(i) + '.txt','w') as manifest:
        manifest.write('\n'.join(subset))


queue = [(c,'between') for c in itertools.combinations(range(len(subsets)), 2)]
queue.extend([(i,'within') for i in range(len(subsets))])

running = []
while queue:
    for proc in running:
        if proc.poll() is not None:
            sys.stdout.write(str(proc.pid) + ' finished '  + time.asctime(time.localtime(time.time())) + '\n')
            sys.stdout.flush()
            running.remove(proc)
    if len(running) < cores:
        running.append(job(queue[-1][0],queue[-1][1]))
        queue.pop(-1)
    else:
        time.sleep(10)
for proc in running:
    proc.wait()
    sys.stdout.write(str(proc.pid) + ' finished '  + time.asctime(time.localtime(time.time())) + '\n')
    sys.stdout.flush()

os.system('rm -r {d}'.format(d = tmp_dir))

