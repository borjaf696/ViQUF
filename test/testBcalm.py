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
        # Primer lanzamiento para definir la abundancia
        args = {'in':path, 'kmerSize':kmerSize,'abundanceMin':abundanceMin,'out':self.outFile}
        cmd = [self.exeBcalm,'-in',args['in'],'-histo','1','-kmer-size',args['kmerSize'],'-abundance-min',args['abundanceMin'],'-out',args['out']]
        print('Bcalm cmd (first launch): ',cmd)
        Utils.executecmd(cmd)
        args['abundanceMin'] = str(self.__study_frequency_histograms(int(abundanceMin)))
        self.__produceGraphFile(args)
        self.__indexGFA(self.outFile)

        self._g.exportGraph(self.outFile+self._graphExt)
        print('End!')

    def __study_frequency_histograms(self, abundanceMin = 1):
        file = self.outFile+'.histo'
        x,y = [],[]
        with open(file, 'r') as f:
            for i, line in enumerate(f.readlines()):
                if i <= abundanceMin:
                    continue
                l_split = line.strip().split('\t')
                x.append(int(l_split[0]))
                y.append(int(l_split[1]))
        kneedle = KneeLocator(x, y, S=1.0, curve='concave', direction='increasing')
        print('Recommended abundance: ', abundanceMin+kneedle.knee)
        return abundanceMin + kneedle.knee

    def getOutFile(self):
        return self.outFile+self._tail

    def __produceGraphFile(self, args):
        tail = '.unitigs.fa'
        cmd = [self.exeBcalm,'-in',args['in'],'-histo','1','-kmer-size',args['kmerSize'],'-abundance-min',args['abundanceMin'],'-out',args['out']]
        print('Bcalm cmd: ',cmd)
        Utils.executecmd(cmd)
        cmd = ['python',self.exeGFA,args['out']+tail,args['out'],args['kmerSize']]
        print('To GFA cmd: ', cmd)
        Utils.executecmd(cmd)
        print('To FM')
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
    tmpDir, resDir = 'tmp/', 'tmpresultsDir_'
    readFiles = []
    def __preprocess_karect(path):
        global resDir

        resDir = resDir + 'Karect/'
        Utils.remove_dir(resDir)
        Utils.mkdir(resDir)
        files = Utils.get_files(path, ['fastq'])
        files.sort()
        print('Files: ', files)
        cmd = ['karect','-correct','-matchtype=hamming','-celltype=haploid','-resultdir='+resDir]+['-inputfile='+t for t in files]
        Utils.executecmd(cmd)
        Utils.remove_dir(tmpDir)
        Utils.mkdir(tmpDir)
        suffix = 'Ownlatest/'
        files, outputFile = Utils.get_files(resDir,['fastq']), tmpDir+suffix+'append.fasta'
        files.sort()
        print('New files corrected: ', files)
        newFiles = BioUtils.renameFastqSeqs(files, tmpDir)
        print('NewFiles: ', newFiles)
        Utils.mkdir(tmpDir+suffix)
        return [Utils.append_files(newFiles, outputFile), newFiles[0]]

    def __preprocess_tgs(path, path_pe = None, kmer_size = 30, technique = 'standard'):
        global resDir
        if path_pe is None:
            resDir += 'Consent/'
            suffix = 'CONSENT_corrected.fasta'
            Utils.remove_dir(resDir)
            Utils.mkdir(resDir)
            files = Utils.get_files(path)
            print('Files: ',files)
            cmd = ['third-party/CONSENT/CONSENT-correct','--in',files[0],'--out',resDir+suffix,'--type','PB']
            Utils.executecmd(cmd)
            Utils.remove_dir(tmpDir)
            Utils.mkdir(tmpDir)
        if path_pe is not None:
            RATE_KMER, KMER_USAGE, SOLID_THRESHOLD, TRIALS, ITERATIONS = kmer_size,kmer_size, 2, 1, 3
            iterative, c_iterations = False, 0
            if technique == 'standard':
                ITERATIONS = 1
            while c_iterations < ITERATIONS:
                if c_iterations < 1:
                    resDir += 'Lordec/'
                suffix = 'lordec_corrected.fasta'
                Utils.remove_dir(resDir)
                Utils.mkdir(resDir)
                files = Utils.get_files(path)
                print('Files: ',files)
                files_pe = Utils.get_files(path_pe)
                files_append = (' ').join(files_pe)
                print('Files append: ', files_append)
                Utils.remove_file(files_pe[0]+'_k'+str(KMER_USAGE)+'_s'+str(SOLID_THRESHOLD)+'.h5')
                cmd = ['lordec-correct','-t',str(TRIALS),'-2']+files_pe+['-i',files[0], '-k', str(KMER_USAGE),'-s',str(SOLID_THRESHOLD),'-o',resDir+suffix]
                Utils.executecmd(cmd)
                Utils.remove_dir(tmpDir)
                Utils.mkdir(tmpDir)
                c_iterations = c_iterations + 1
                files[0] = resDir+suffix
                # trim_split
                suffix_trim_split = 'lordec_corrected_trim_split.fasta'
                cmd = ['lordec-trim-split','-i',files[0],'-o',resDir+suffix_trim_split]
                Utils.executecmd(cmd)
                files[0] = resDir + suffix_trim_split
                #KMER_USAGE = str(int(KMER_USAGE)+int(RATE_KMER))
        suffix = 'Ownlatest/'
        files = Utils.get_files(resDir)
        print('New files corrected: ', files)
        Utils.mkdir(tmpDir+suffix)
        return files

    def __create_s_pe(tgs_file):
        global tmpDir
        insert_size = 700
        read_size = 350

        dir_reads = tmpDir+'reads/'
        Utils.mkdir(dir_reads)

        files = [dir_reads+'read1.fasta', dir_reads+'read2.fasta']
        [Utils.remove_file(f) for f in files]
        print('New paired end files in: ', files)
        with open(files[0], 'w+') as f1, open(files[1],'w+') as f2, open(tgs_file,'r') as f_read:
            id = ''
            lines = f_read.readlines()
            with progressbar.ProgressBar(max_value=len(lines)) as bar:
                for tgs_read,line in enumerate(lines):
                    bar.update(tgs_read)
                    if line[0] == '>':
                        id = line.strip()
                    if line[0] != '>':
                        read_process = line.strip()
                        for num_read,i in enumerate(range(0,len(read_process) - insert_size + 1)):
                            # Left read
                            f1.write(id+'.'+str(num_read)+'.1\n')
                            f1.write(read_process[i:i+read_size]+'\n')
                            # Right read
                            new_start = i + (insert_size - read_size)
                            f2.write(id+'.'+str(num_read)+'.2\n')
                            f2.write(BioUtils.reverse_complement(read_process[new_start:new_start+read_size])+'\n')

    def __preprocess2(files):
        suffix = 'Ownlatest/'
        outputFile = tmpDir+suffix+'merge.fasta'
        return Utils.append_files(files, outputFile)

    type = sys.argv[4]
    if type == 'ngs':
        pathIn = __preprocess_karect(sys.argv[1])
        path, kmerSize, abundanceMin = pathIn[0], sys.argv[2], sys.argv[3]
        print('Path: ', path)
        rG = RepresentantGraph(path, kmerSize, abundanceMin)
    elif type == 'tgs':
        method = 'lordec'
        if method == 'consent':
            pathIn = __preprocess_tgs(sys.argv[1])
        elif method == 'lordec':
            '''
            We demand paired-end reads to be read1.fastq, read2.fastq, here we receive the folder
            '''
            # Iterative execution
            technique = sys.argv[5]
            pathIn = __preprocess_tgs(sys.argv[1], sys.argv[2], sys.argv[3], technique=technique)
        #pathIn = ['/home/bfreire/Gatb-trial/tmpresultsDir_Consent/CONSENT_corrected.fasta']
        __create_s_pe(pathIn[0])


    #BioUtils.identicalClustering(__preprocess2([pathIn[1],rG.getOutFile()]))
