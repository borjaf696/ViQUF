#Borja :)
import sys, os
sys.path.append('')
import subprocess
import altair as alt
import pandas as pd
import progressbar
from shutil import copyfile
from Bio.Seq import Seq
from utils.utils import *
import numpy as np
import math
from mip import *
from gurobipy import *

def __process_contigs(file, seqs_file):
    unitigs_seqs = dict()
    with open(seqs_file, 'r') as f:
        for i,line in enumerate(f.readlines()):
            if line[0] == '>':
                num_unitig = line.strip().split('>')[1]
            else:
                unitigs = line.strip().split('-')
                unitigs_seqs[int(num_unitig)] = [int(t) for t in unitigs[0:len(unitigs) - 1]]
    ids, seqs = BioUtils.readFasta(file)
    contigs = dict()
    for i, id in enumerate(ids):
        seq = seqs[i]
        num, hat_freq = id.split('-')[0],id.split('-')[-1]
        contigs[int(num)] = int(hat_freq)
    return unitigs_seqs, contigs

def __process_assembly_graph(file):
    nodes = dict()
    with open(file, 'r') as f:
        for line in f.readlines():
            node_split = line.strip().split('-')
            nodes[int(node_split[1])] = [float(node_split[2]),float(node_split[3])]
    return nodes

def __build_mixed_lp(unitigs_seqs, nodes, contigs):
    # Assign nodes freqs and so on
    m = Model('qp')
    x = dict()
    unitigs_contig_map = dict()
    for key, val in unitigs_seqs.items():
        x[key] = m.addVar()
        m.addConstr(x[key] >= 0)
        for u in val:
            if u in unitigs_contig_map.keys():
                unitigs_contig_map[u].append(key)
            else:
                unitigs_contig_map[u] = [key]
    objective = QuadExpr()
    for u, belonging_contigs in unitigs_contig_map.items():
        #if len(belonging_contigs) > 1:
        #    continue
        sum_strains = 0
        for val in belonging_contigs:
            sum_strains += contigs[val]*x[val]
        objective += ((nodes[u][1] - sum_strains)*(nodes[u][1] - sum_strains)/nodes[u][1])*nodes[u][0]
        #print(objective)
    m.setObjective(objective,GRB.MINIMIZE)
    m.write('graphs/mlp.lp')
    m.optimize()
    xs = []
    for v in m.getVars():
        print(v.varName, v.x)
        xs.append(v.x)
    print('Obj:', m.objVal)
    #New freqs
    place = 0
    for key, val in contigs.items():
        print(key,' Old freq: ',val,' New freq: ', val*xs[place])
        place += 1

if __name__ == '__main__':
    fasta_files, unitigs_seq_file = 'tmp/unitigs-viaDBG-nf-std.fasta', 'tmp/unitigs-viaDBG-nf-std.fasta.unitigs'
    nodes = __process_assembly_graph('graphs/dbg_postsubsane.txt.basic')
    unitigs_seqs, contigs = __process_contigs(fasta_files, unitigs_seq_file)
    __build_mixed_lp(unitigs_seqs, nodes, contigs)
    