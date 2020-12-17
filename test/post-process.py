#Borja :)
import sys, os
sys.path.append('')
import subprocess
import altair as alt
import pandas as pd
import progressbar
from shutil import copyfile
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from utils.utils import *
from kneed import KneeLocator
from sklearn.cluster import MeanShift
import numpy as np
import math

def __process_clustering(df, algorithm, columns = None,args = None, kwds = None, ):
    assert columns is not None
    X = df[columns].values
    labels = algorithm(**kwds).fit_predict(X)
    return labels

def __split_fastas(file, labels):
    ids, seqs = BioUtils.readFasta(file)
    print(np.unique(labels))
    fastas = dict()
    for j, label in enumerate(labels):
        if label in fastas.keys():
            fastas[label][ids[j]] = str(seqs[j])+'\n'
        else:
            fastas[label] = {ids[j]:str(seqs[j])+'\n'}
    for key, val in fastas.items():
        BioUtils.export_fasta_from_dict(val, 'tmp/unitigs_'+str(key)+'.fasta')

if __name__ == '__main__':
    df = pd.read_csv('tmp/flows.csv')
    # Perform clustering (no number of clusters defined)
    labels = __process_clustering(df, MeanShift, columns = ['Freqs'], kwds = {'cluster_all':True})
    # Split fastas
    fasta_files = 'tmp/unitigs-viaDBG-nf.fasta'
    __split_fastas(fasta_files, labels)