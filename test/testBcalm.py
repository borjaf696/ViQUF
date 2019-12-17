#Borja :)
import sys, os
sys.path.append('')
import subprocess
import altair as alt
import pandas as pd
from shutil import copyfile
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from utils.utils import *

class UtilsReport:
    @staticmethod
    def exportHistogram(histogram, path):
        x = range(len(histogram))
        df = pd.DataFrame({'X':x,'Y':histogram})
        chart = alt.Chart(df).mark_bar().encode(alt.X('X',  bin=alt.Bin(maxbins=100)), y = 'Y')
        chart.save(path)


class RepresentantGraph:
    exeBcalm = '/home/bfreire/Gatb-trial/third-party/bcalm/build/bcalm'
    exeGFA = '/home/bfreire/Gatb-trial/third-party/bcalm/scripts/convertToGFA.py'
    outFile, _tail, _graphExt= 'tmp/unitigs', '.unitigs.fa', '.graph'
    seqs = []
    class GraphStruct:
        def __init__(self):
            self._graph = dict()

        def addVertex(self, u):
            self._graph[u] = []

        def addEdge(self, u, v):
            self._graph[u].append(v)

        def exportGraph(self, outputFile):
            with open(outputFile, 'w+') as fWrite:
                numKeys = len(self._graph.keys())
                fWrite.write(str(numKeys)+'\n')
                for key, val in self._graph.items():
                    fWrite.write(str(key)+' ')
                    for v in val:
                        fWrite.write(str(v)+' ')
                    fWrite.write('\n')

    def __init__(self, path = None, kmerSize = 30, abundanceMin = 1):
        self._g, self._kmerSize = self.GraphStruct(), int(kmerSize)
        self.__produceGraphFile({'in':path, 'kmerSize':kmerSize,'abundanceMin':abundanceMin,'out':self.outFile})
        self.__indexGFA(self.outFile)

        self._g.exportGraph(self.outFile+self._graphExt)
        print('End!')

    def getOutFile(self):
        return self.outFile+self._tail

    def __produceGraphFile(self, args):
        tail = '.unitigs.fa'
        cmd = [self.exeBcalm,'-in',args['in'],'-kmer-size',args['kmerSize'],'-abundance-min',args['abundanceMin'],'-out',args['out']]
        print('Bcalm cmd: ',cmd)
        Utils.executecmd(cmd)
        cmd = ['python',self.exeGFA,args['out']+tail,args['out'],args['kmerSize']]
        print('ToGFA cmd: ', cmd)
        Utils.executecmd(cmd)
        print('ToFm')
        BioUtils.fastToFm(args['out']+tail,args['out']+'.FM')

    def __indexGFA(self, file):
        numSeqs, min, max = 0, 9999999, 0
        histogram = [0]*10000
        with open(file, 'r') as f:
            for line in f.readlines():
                if line[0] == 'S':
                    dnaSeq = line.split('\t')[2]
                    unitigLength, pivote = len(dnaSeq), int(len(dnaSeq)/2)
                    self.seqs.append(Seq(dnaSeq[pivote:pivote+self._kmerSize], generic_dna))
                    numSeqs += 1
                    if unitigLength > len(histogram):
                        histogram = histogram + [0]*(unitigLength - len(histogram) + 1)
                    histogram[unitigLength] += 1
                    min = unitigLength if min > unitigLength else min
                    max = unitigLength if max < unitigLength else max
            print('Number of sequences: ', numSeqs)
            print('Exporting histogram:')
            histogram = histogram[0:max+1]
            UtilsReport.exportHistogram(histogram, 'stats/histogram.html')
            print('Histogram available!\nTotal Sequences: ',2*numSeqs)
            self._unitigs = numSeqs
            for i in range(0,2*numSeqs):
                self._g.addVertex(i)
        with open(file, 'r') as f:
            for line in f.readlines():
                if line[0] == 'L':
                    infoLine = line.split('\t')
                    ori, target = int(infoLine[1]) if infoLine[2] == '+' else int(infoLine[1])+numSeqs\
                        , int(infoLine[3]) if infoLine[4] == '+' else int(infoLine[3])+numSeqs
                    if ori > 2*numSeqs:
                        print("Ori higher: ", ori)
                    self._g.addEdge(ori, target)

if __name__=='__main__':
    print('Lets do this')
    tmpDir, resDir = 'tmp/', 'tmpresultsDir/'
    readFiles = []
    def __preprocess(path):
        Utils.mkdir(resDir)
        files = Utils.get_files(path, ['fastq'])
        files.sort()
        print('Files: ', files)
        cmd = ['karect','-correct','-matchtype=hamming','-celltype=haploid','-resultdir='+resDir]+['-inputfile='+t for t in files]
        Utils.executecmd(cmd)
        Utils.mkdir(tmpDir)
        suffix = 'Ownlatest/'
        files, outputFile = Utils.get_files(resDir,['fastq']), tmpDir+suffix+'append.fasta'
        files.sort()
        print('Files: ', files)
        newFiles = BioUtils.renameFastqSeqs(files, tmpDir)
        print('NewFiles: ', newFiles)
        Utils.mkdir(tmpDir+suffix)
        return [Utils.append_files(newFiles, outputFile), newFiles[0]]

    def __preprocess2(files):
        suffix = 'Ownlatest/'
        outputFile = tmpDir+suffix+'merge.fasta'
        return Utils.append_files(files, outputFile)

    pathIn = __preprocess(sys.argv[1])
    path, kmerSize, abundanceMin = pathIn[0], sys.argv[2], sys.argv[3]
    print('Path: ', path)
    rG = RepresentantGraph(path, kmerSize, abundanceMin)

    #BioUtils.identicalClustering(__preprocess2([pathIn[1],rG.getOutFile()]))
