#Borja :)
import os, subprocess, sys
import altair as alt
import pandas as pd
from shutil import copyfile
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

class Utils:
    @staticmethod
    def cpfile(source, dest):
        copyfile(source, dest)

    @staticmethod
    def getcwd():
        return os.getcwd()

    @staticmethod
    def chdir(dir):
        os.chdir(dir)

    @staticmethod
    def mkdir(dir):
        if not os.path.exists(dir):
            os.mkdir(dir)
        else:
            print('Directory: ',dir,' already exists')

    @staticmethod
    def exist_dir(d):
        return os.path.isdir(d)

    @staticmethod
    def remove_file(file):
        os.remove(file)

    @staticmethod
    def get_dirs(path, bound = None):
        if bound is None:
            return [path+'/'+dI for dI in os.listdir(path) if os.path.isdir(os.path.join(path,dI))]
        else:
            dirs = []
            for dI in os.listdir(path):
                if len(dirs) == bound:
                    return dirs
                if os.path.isdir(os.path.join(path, dI)):
                    dirs.append(path+'/'+dI)
            return dirs

    @staticmethod
    def get_files(path, extension = ['fq','fastq','fasta','fa'], content = None):
        if content is None:
            return [path+'/'+t for t in os.listdir(path) if not os.path.isdir(path+'/'+t) and t.split('.')[-1] in extension]
        else:
            return [path+'/'+t for t in os.listdir(path) if not os.path.isdir(path+'/'+t) and content in t and t.split('.')[-1] in extension]

    @staticmethod
    def get_files_recursive(path, extension = ['fq','fa']):
        paths, results = [path], []
        while len(paths) > 0:
            p = paths[0]
            paths = paths[1:len(paths)]
            for t in os.listdir(p):
                n_path = p+'/'+t
                if os.path.isdir(n_path):
                    paths.append(n_path)
                elif t.split('.')[-1] in extension:
                    results.append(n_path)
        return results

    @staticmethod
    def get_files_recursive_content(path, content, threshold = None, avoid = None):
        paths, results = [path], []
        while len(paths) > 0:
            p = paths[0]
            paths = paths[1:len(paths)]
            for t in os.listdir(p):
                n_path = p + '/' + t
                if os.path.isdir(n_path):
                    paths.append(n_path)
                elif content in n_path and avoid is not None and avoid not in n_path:
                    results.append(n_path)
                elif content in n_path and avoid is None:
                    results.append(n_path)
            if threshold is not None and len(results) > threshold:
                return results
        return results

    @staticmethod
    def append_files(list_files, output_file):
        with open(output_file,'w+') as out_file:
            for f in list_files:
                with open(f,'r') as f_read:
                    for line in f_read.readlines():
                        out_file.write(line)

    @staticmethod
    def executecmd(args, out = None):
        '''
        :param path: LIST of args
        :return:
        '''
        if out is None:
            subprocess.call(args)
        else:
            with open(out, 'w+') as fpout:
                subprocess.call(args, stdout = fpout)

    @staticmethod
    def write_list(lst,output, type = 'int'):
        with open(output, 'w+') as f:
            for l in lst:
                if type != 'int':
                    line = ','.join(map(lambda x:str(x),l))
                    f.write(line+'\n')
                else:
                    f.write(str(l)+'\n')

    @staticmethod
    def export_dict(dictionary,name, format = 'csv'):
        assert format == 'csv'

        if format == 'csv':
            with open(name, 'w+') as f:
                for key, val in dictionary.items():
                    f.write(str(key))
                    for k,v in val.items():
                        f.write(','+str(v))
                    f.write('\n')

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
    outFile = 'unitigs'
    seqs = []
    class GraphStruct:
        def __init__(self):
            self.graph = dict()

        def addVertex(self, u):
            self.graph[u] = []

        def addEdge(self, u, v):
            self.graph[u].append(v)

    def __init__(self, path = None, kmerSize = 30, abundanceMin = 1):
        self._g, self._kmerSize = self.GraphStruct(), int(kmerSize)
        self.__produceGraphFile({'in':path, 'kmerSize':kmerSize,'abundanceMin':abundanceMin,'out':self.outFile})
        self.__indexGFA(self.outFile)
        print('End!')


    def __produceGraphFile(self, args):
        tail = '.unitigs.fa'
        cmd = [self.exeBcalm,'-in',args['in'],'-kmer-size',args['kmerSize'],'-abundance-min',args['abundanceMin'],'-out',args['out']]
        print('Bcalm cmd: ',cmd)
        Utils.executecmd(cmd)
        cmd = ['python',self.exeGFA,args['out']+tail,args['out'],args['kmerSize']]
        print('ToGFA cmd: ', cmd)
        Utils.executecmd(cmd)

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
                    histogram[unitigLength] += 1
                    min = unitigLength if min > unitigLength else min
                    max = unitigLength if max < unitigLength else max
            print('Exporting histogram:')
            histogram = histogram[min-1:max+1]
            UtilsReport.exportHistogram(histogram, 'stats/histogram.html')
            print('Histogram available!\nTotal Sequences: ',2*numSeqs)
            self._unitigs = numSeqs
            for i in range(0,2*numSeqs):
                self._g.addVertex(i)
            for line in f.readlines():
                if line[0] == 'L':
                    print('L: ', line)
                    infoLine = line.split('\t')
                    ori, target = int(infoLine[1]) if infoLine[2] == '+' else int(infoLine[1])+numSeqs\
                        , int(infoLine[3]) if infoLine[4] == '+' else int(infoLine[3])+numSeqs
                    self.__g.addEdge(ori, target)
        print(self.seqs[0:25])

if __name__=='__main__':
    print('Lets do this')
    path, kmerSize, abundanceMin = sys.argv[1], sys.argv[2], sys.argv[3]
    RepresentantGraph(path, kmerSize, abundanceMin)
