#Borja :)
import sys, os
sys.path.append('')
import subprocess
import pandas as pd
import altair as alt
import progressbar
from shutil import copyfile
from Bio.Seq import Seq
from utils.utils import *
import numpy as np
import math

class UtilsReport:
    @staticmethod
    def exportHistogram(histogram, path):
        x = range(len(histogram))
        df = pd.DataFrame({'X':x,'Y':histogram})
        chart = alt.Chart(df).mark_bar().encode(alt.X('X',  bin=alt.Bin(maxbins=100)), y = 'Y')
        chart.save(path)

class RepresentantGraph:
    exeBcalm = 'third-party/bcalm/build/bcalm'
    exeGFA = 'third-party/bcalm/scripts/convertToGFA.py'
    outFile, _tail, _graphExt= 'tmp/unitigs', '.unitigs.fa', '.graph'
    seqs = []
    class GraphStruct:
        def __init__(self):
            self._graph = dict()
            self._freqs = dict()
            self._active = dict()
            self._lengths = dict()

        def addVertex(self, u):
            self._graph[u] = []
            self._freqs[u] = 0

        def addEdge(self, u, v):
            self._graph[u].append(v)

        def addActive(self, u, active):
            self._active[u] = active

        def setActive(self, u, active):
            self._active[u] = active

        def addFreq(self, u, frec):
            self._freqs[u] = frec

        def addLength(self, u, length):
            self._lengths[u] = length

        def getActive(self, u):
            return self._active[u]

        def exportGraph(self, outputFile, args):
            with open(outputFile, 'w+') as fWrite:
                numKeys = len(self._graph.keys())
                fWrite.write(str(numKeys)+' '+str(args['abundanceMin'])+' \n')
                for key, val in self._graph.items():
                    if not self._active[key]:
                        continue
                    fWrite.write(str(key)+' ')
                    for v in val:
                        fWrite.write(str(v)+' ')
                    fWrite.write(str(self._freqs[key])+' ')
                    fWrite.write(str(self._lengths[key])+' \n')

    def __init__(self, path = None, kmerSize = 30, abundanceMin = 1, meta = False, filtering = True, output_prefix = 'tmp/unitigs'):
        self._g, self._kmerSize = self.GraphStruct(), int(kmerSize)
        # Primer lanzamiento para definir la abundancia
        self.outFile = output_prefix
        args = {'in':path, 'kmerSize':kmerSize,'abundanceMin':abundanceMin,'out':self.outFile}
        cmd = [self.exeBcalm,'-in',args['in'],'-histo','1','-kmer-size',args['kmerSize'],'-abundance-min',args['abundanceMin'],'-out',args['out']]
        print('Bcalm cmd (first launch): ',cmd)
        Utils.executecmd(cmd)
        args['abundanceMin'] = str(self.__study_frequency_histograms(int(abundanceMin), meta = meta, filtering = filtering))
        Utils.remove_file(self.outFile)
        self.__produceGraphFile(args)
        self.__indexGFA(self.outFile, args['abundanceMin'])
        self._g.exportGraph(self.outFile+self._graphExt, args)
        print('End!')

    def __study_frequency_histograms(self, abundanceMin = 1, meta = False, filtering = True):
        if meta:
            return 3
        if not filtering:
            return 0
        file, file_w = self.outFile+'.histo', self.outFile+'.histo.txt'
        def roundup(x):
            return int(math.ceil(x / 10.0)) * 10
        def __vector_trend(x, window_size = 15):
            assert len(x) > 2*window_size

            for i, v in enumerate(x):
                if i < window_size:
                    continue
                left_trend, right_trend = 0, 0
                for j in range(i-window_size, i):
                    left_trend += (x[j] >= v)
                for j in range(i, min(i+window_size, len(x))):
                    right_trend += (x[j] >= v)
                if left_trend <= right_trend and right_trend >= (window_size*0.5):
                    print('LT: ',left_trend,' RT: ',right_trend)
                    return i
        def __kernel_estimation(data_x, data_y):
            import matplotlib.pyplot as plt
            from datetime import datetime
            SCALE = 0.25 if meta else 1.0
            data_y_tmp = data_y[np.where(data_y < 1000)]
            data_x_tmp = data_x[0:max(data_y_tmp)]
            plt.hist(data_y, bins=125, density=True)
            plt.savefig(self.outFile+'_histo.png')
            plt.clf()
            from scipy import stats
            plt.hist(data_y_tmp, bins=100, density=True)
            gkde=stats.gaussian_kde(data_y)
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print('Kde factor: ', gkde.factor,' ',current_time)
            gkde.set_bandwidth(bw_method=gkde.factor*SCALE)
            #data_x_tmp = data_x[0:100] if meta else data_x
            kdepdf = gkde.evaluate(data_x_tmp)
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print('End evaluation. ', current_time)
            # plot estimated density
            plt.plot(data_x_tmp, kdepdf, label='kde', color="g")
            # Cross zeros
            Z = np.reshape(kdepdf.T, data_x_tmp.shape)
            diff = np.gradient(Z)
            sdiff = np.sign(diff)
            zc = np.where(sdiff[:-1] != sdiff[1:])
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print('Zero crosses: ',zc,' ',current_time)
            decision = zc[0][0]
            for i in range(len(zc[0])-1):
                first, second, third = zc[0][i], zc[0][i+1], zc[0][i+2]
                if gkde.evaluate(first) > gkde.evaluate(second) and gkde.evaluate(second) < gkde.evaluate(third):
                    decision = int((first+second)*0.5)
                    break
            plt.vlines(decision, plt.ylim()[0], plt.ylim()[1], linestyles='dashed')
            plt.title('Kernel Density Estimation')
            plt.savefig(self.outFile+'_density.png')
            return decision

        UPPER_FREQ_LIMIT = 999999 if meta else 1000
        x,y = [],[]
        X,Y = [],[]
        with open(file, 'r') as f,open(file_w, 'w') as f2:
            for i, line in enumerate(f.readlines()):
                l_split = line.strip().split('\t')
                value = UPPER_FREQ_LIMIT if int(l_split[1]) > UPPER_FREQ_LIMIT else int(l_split[1])
                if int(l_split[0]) < UPPER_FREQ_LIMIT:
                    X.append(int(l_split[0]))
                    Y += [int(l_split[0])]*value
                    if i <= abundanceMin:
                        continue
                x.append(int(l_split[0]))
                y.append(roundup(int(l_split[1])))
                f2.write(str(roundup(int(l_split[1])))+'\n')
        #kneedle = KneeLocator(x, y, S=2.0, curve='convex', direction='decreasing')
        offset = __vector_trend(y)
        print('Offset: ',offset)
        #offset = kneedle.knee
        min_freq_estimator = abundanceMin + offset
        min_freq_estimator_kde = __kernel_estimation(np.array(X), np.array(Y))
        print('Recommended abundance (kernel estimator): ',min_freq_estimator_kde)
        min_freq_estimator = min(min_freq_estimator_kde, min_freq_estimator) if meta else max(min_freq_estimator_kde, min_freq_estimator)
        #print('Handmade abundance: ', 16)
        #min_freq_estimator = 16
        print('Recommended abundance: ', min_freq_estimator)
        return min_freq_estimator

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
        #print('To FM')
        #BioUtils.fastToFm(args['out']+tail,args['out']+'.FM')

    def __indexGFA(self, file, abundanceMin):
        numSeqs, min, max = 0, 9999999, 0
        histogram = [0]*10000
        frecs, lengths, active = [], [], []
        with open(file, 'r') as f:
            lines = f.readlines()
            connect, is_connected = [False]*len(lines), [False]*len(lines)
            for line in lines:
                if line[0] == 'S':
                    line_split = line.split('\t')
                    frecs.append(line_split[5].split(':')[2].strip())
                    unitigLength = int(line_split[3].split(':')[2].strip())
                    totalFreq = int(line_split[4].split(':')[2].strip())
                    lengths.append(unitigLength)
                    if float(totalFreq / unitigLength) >= float(abundanceMin):
                        active.append(True)
                    else:
                        active.append(False)
                    numSeqs += 1
                    if unitigLength > len(histogram):
                        histogram = histogram + [0]*(unitigLength - len(histogram) + 1)
                    histogram[unitigLength] += 1
                    min = unitigLength if min > unitigLength else min
                    max = unitigLength if max < unitigLength else max
                elif line[0] == 'L':
                    infoLine = line.split('\t')
                    ori =  int(infoLine[1])
                    if infoLine[2] == '+':
                        connect[ori] = True
                    else:
                        is_connected[ori] = True
            print('Number of sequences: ', numSeqs)
            print('Exporting histogram:')
            histogram = histogram[0:max+1]
            UtilsReport.exportHistogram(histogram, 'stats/histogram.html')
            print('Histogram available!\nTotal Sequences: ',2*numSeqs)
            self._unitigs = numSeqs
            for i in range(0,2*numSeqs):
                self._g.addVertex(i)
                if i >= numSeqs:
                    self._g.addFreq(i, frecs[i - numSeqs])
                    self._g.addActive(i, active[i - numSeqs] or (connect[i - numSeqs] and is_connected[i - numSeqs]))
                    self._g.addLength(i, lengths[i - numSeqs])
                else:
                    self._g.addFreq(i, frecs[i])
                    self._g.addActive(i, active[i] or (connect[i] and is_connected[i]))
                    self._g.addLength(i, lengths[i])
        with open(file, 'r') as f:
            for line in f.readlines():
                if line[0] == 'L':
                    infoLine = line.split('\t')
                    ori, target = int(infoLine[1]) if infoLine[2] == '+' else int(infoLine[1])+numSeqs\
                        , int(infoLine[3]) if infoLine[4] == '+' else int(infoLine[3])+numSeqs
                    if self._g.getActive(ori) and self._g.getActive(target):
                        if ori > 2*numSeqs:
                            print("Ori higher: ", ori)
                        self._g.addEdge(ori, target)
MAX_AMPLICON = 4
class AmpliconsGraph:
    _nodes_map = dict()
    _nodes = []
    _active_nodes = []
    _edges = []
    _lengths = []
    _freqs = []

    def __init__(self, files, files_edges, nodes_deactivated):
        for i, file in enumerate(files):
            if i >= MAX_AMPLICON:
                break
            self.addVertices(file, i, nodes_deactivated[i])
            self.addEdges(file, i)
        for i, file in enumerate(files_edges):
            if i >= MAX_AMPLICON:
                break
            self.addExtraEdges(file, i, i+1)

    def numVertices(self):
        return len(self._nodes)

    def numEdges(self):
        edges = 0
        for set_edges in self._edges:
            edges += len(set_edges)
        return edges

    def addVertices(self, file, amplicon = 0, nodes_deactivated = None):
        num_vertices = len(self._nodes)
        with open(file, 'r+') as f:
            for line in f.readlines():
                line_split = line.strip().split(' ')
                l_line_split = len(line_split)
                if int(line_split[0]) in nodes_deactivated:
                    continue
                self._nodes_map[(amplicon,int(line_split[0]))] = num_vertices
                self._nodes.append((amplicon,int(line_split[0])))
                self._lengths.append(int(line_split[l_line_split - 1]))
                self._freqs.append(float(line_split[l_line_split - 2]))
                self._edges.append([])
                num_vertices += 1
    
    def addEdges(self, file, amplicon = 0):
        with open(file, 'r+') as f:
            for line in f.readlines():
                line_split = line.strip().split(' ')
                l_line_split = len(line_split)
                for edge in line_split[1:-2]:
                    if ((amplicon, int(line_split[0])) not in self._nodes_map.keys()) or ((amplicon, int(edge)) not in self._nodes_map.keys()):
                        continue
                    self._edges[self._nodes_map[(amplicon, int(line_split[0]))]].append(self._nodes_map[(amplicon, int(edge))])
    
    def addExtraEdges(self, file, amplicon_left = 0, amplicon_right = 1):
        with open(file, 'r+') as f:
            for line in f.readlines():
                line_split = line.strip().split(' ')
                if ((amplicon_left,int(line_split[0])) not in self._nodes_map.keys()) or ((amplicon_right,int(line_split[1])) not in self._nodes_map.keys()):
                        continue
                node_left, node_right = self._nodes_map[(amplicon_left,int(line_split[0]))], self._nodes_map[(amplicon_right,int(line_split[1]))]
                self._edges[node_left].append(node_right)

    def show_graph(self):
        for i, edges in enumerate(self._edges):
            print('Node: ', i,' ',self._nodes[i],' Length: ', self._lengths[i], ' Freq: ',self._freqs[i])
            print('Edges: ')
            for edge in edges:
                print(edge)

    def exportGraph(self, outputFile_graph, outputFile_unitigs, files_unitigs):
        print('Exporting graph in: ', outputFile_graph)
        with open(outputFile_graph, 'w+') as fWrite:
            numKeys = len(self._nodes)
            fWrite.write(str(numKeys)+' 10 \n')
            for key, edges in enumerate(self._edges):
                fWrite.write(str(key)+' ')
                for v in edges:
                    fWrite.write(str(v)+' ')
                fWrite.write(str(int(self._freqs[key]))+' ')
                fWrite.write(str(self._lengths[key])+' \n')
        print('Exporting unitigs in: ', outputFile_unitigs)
        with open(outputFile_unitigs,'w+') as fWrite:
            for amplicon,file in enumerate(files_unitigs):
                reads_dict = BioUtils.get_sequences_faf(file)
                for key, val in reads_dict.items():
                    if ((amplicon, int(key)) in self._nodes_map.keys()):
                        fWrite.write('>'+str(self._nodes_map[(amplicon, int(key))])+'\n'+str(val)+'\n')
        # Exporting map
        self._exportMap('/'.join(outputFile_graph.strip().split('/')[:-1]+['map.txt']))
    
    def _exportMap(self, output_file):
        print('Exporting Amplicons/unitigs map: ', output_file)
        with open(output_file, 'w+') as fWrite:
            for key, val in self._nodes_map.items():
                fWrite.write(str(key)+' '+str(val)+'\n')
        

if __name__=='__main__':
    tmpDir, resDir = 'tmp/', 'tmpresultsDir_'
    readFiles = []
    def __preprocess_karect(path, correct = True):
        global resDir

        resDir = resDir + 'Karect/'
        resDir_paired = resDir+'paired/'
        Utils.remove_dir(resDir)
        Utils.mkdir(resDir)
        Utils.mkdir(resDir_paired)
        files = Utils.get_files_recursive(path, ['fastq','fasta'])
        files.sort()
        print('Files: ', files)
        if len(files) > 0:
            extension = files[0].split('.')[-1]
        if correct:
            cmd = ['karect','-correct','-matchtype=hamming','-celltype=haploid','-resultdir='+resDir]+['-inputfile='+t for t in files]
            Utils.executecmd(cmd)
            files = Utils.get_files_recursive(resDir,['.fastq'])
        Utils.remove_dir(tmpDir)
        Utils.mkdir(tmpDir)
        suffix = 'Ownlatest/'
        outputFile = tmpDir+suffix+'append.fasta'
        files.sort()
        print('New files corrected: ', files)
        print('Renaiming sequences/files: ')
        # Opcion A
        # newFiles = BioUtils.renameFastqSeqs(files, tmpDir)
        # Opcion B
        if extension == 'fastq' or extension == 'fq':
            newFiles = [BioUtils.fastqtofasta(f,tmpDir+'/'+str(i)+'.fasta') for i,f in enumerate(files)]
        else:
            newFiles = files
        print('NewFiles: ', newFiles)
        Utils.mkdir(tmpDir+suffix)
        return [Utils.append_files_bash(newFiles, outputFile), newFiles[0]]

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
    print(sys.argv)
    if type == 'ngs':
        pear, correction, meta, filtering = (sys.argv[6] == '--joined'), (sys.argv[5] == '--correct'), (sys.argv[7] == '--meta'), True
        if len(sys.argv) > 8:
            filtering = not (sys.argv[8] == '--no-filter')
        if len(sys.argv) > 9:
            output_prefix = sys.argv[9]
            tmpDir = '/'.join(output_prefix.split('/')[:-1])+'/tmp'
        print('*************** Summary *********************')
        print('Pear: ', pear)
        print('Correct: ', correction)
        print('Meta: ', meta)
        print('Temporary directory: ', tmpDir)
        print('Reads directory: ', sys.argv[1])
        print('********************************************<*')
        pathIn = __preprocess_karect(sys.argv[1], correction)
        path, kmerSize, abundanceMin = pathIn[0], sys.argv[2], sys.argv[3]
        print('*********************************************')
        print('Abundance min (not in use): ', abundanceMin)
        print('Path: ', path)
        print('Kmer-size: ', kmerSize)
        print('*********************************************')
        rG = RepresentantGraph(path, kmerSize, abundanceMin, meta = meta, filtering = filtering, output_prefix = tmpDir+'/unitigs')
    elif type == 'tgs':
        print(sys.argv)
        method = 'hifi' if 'hifi' in sys.argv else 'lordec'
        if method == 'consent':
            pathIn = __preprocess_tgs(sys.argv[1])
        elif method == 'lordec':
            '''
            We demand paired-end reads to be read1.fastq, read2.fastq, here we receive the folder
            '''
            # Iterative execution
            technique = sys.argv[5]
            pathIn = __preprocess_tgs(sys.argv[1], sys.argv[2], sys.argv[3], technique=technique)
            __create_s_pe(pathIn[0])
        elif method == 'hifi':
            print('Hifi reads High Fidelity reads - correction is no need')
            path, kmerSize, abundanceMin = sys.argv[1], sys.argv[2], sys.argv[3]
            print('*********************************************')
            print('Path: ', path)
            print('Kmer-size: ', kmerSize)
            print('*********************************************')
            rG = RepresentantGraph(path, kmerSize, abundanceMin, meta = False, filtering = True, output_prefix = tmpDir+'/unitigs')
        #pathIn = ['/home/bfreire/Gatb-trial/tmpresultsDir_Consent/CONSENT_corrected.fasta']
    elif type == 'amplicons':
        KMER_SIZE = 65
        def _create_amplicons_segments(df_amplicons, reference, reads, output_dir_tmp = 'tmp_amplicons/', reference_file = None):
            # Splits, makes and align the reads to each segment the index for each amplicon
            def __build_index(ref_file, index_prefix, program = 'bowtie2'):
                BioUtils.build_align_index(ref_file, index_prefix, program)
                return index_prefix
            def __align_reads(reads, index_prefix, output_path, sam_file):
                sam_left, sam_right = sam_file+'_1.sam',sam_file+'_2.sam'
                # Left reads alignment and selection
                BioUtils.align_single_end_reads(reads[0],index_prefix, sam_left,permissive=True)
                df_left = pd.read_csv(sam_left, sep = '\t', header = None, names = [i for i in range(19)])
                df_left = df_left.loc[df_left[5] != '*',]
                reads_aligned_left = set(df_left[0])
                output_left = output_path+'_1.fasta'
                # Right reads alignmen and selection
                BioUtils.align_single_end_reads(reads[1],index_prefix, sam_right, permissive=True)
                df_right = pd.read_csv(sam_right, sep = '\t',header = None, names = [i for i in range(19)])
                df_right = df_right.loc[df_right[5] != '*',]
                reads_aligned_right = set(df_right[0])
                # TODO: Chequear si esto es necesario
                reads_set_interest = reads_aligned_left.intersection(reads_aligned_right)
                output_right = output_path+'_2.fasta'
                BioUtils.write_paired_end_reads_by_place(reads, reads_set_interest, (output_left, output_right))
                return  output_left, output_right
            def __align_reads_complete(reads, index_prefix_reference, sam_file, program = 'bowtie2'):
                sam_left, sam_right = sam_file+'_1.sam',sam_file+'_2.sam'
                # Left reads alignment against reference
                BioUtils.align_single_end_reads(reads[0],index_prefix, sam_left, permissive=True,program = program)
                df_left = pd.read_csv(sam_left, sep = '\t', header = None, names = [i for i in range(19)])
                # Right reads alignment against reference
                BioUtils.align_single_end_reads(reads[1],index_prefix, sam_right, permissive=True, program = program)
                df_right = pd.read_csv(sam_right, sep = '\t',header = None, names = [i for i in range(19)])
                return df_left, df_right
            def __distribute_per_amplicon(reads, df_sam, amplicons_limits, amplicon_file_reads, suffix = '_1.fasta', program = 'bowtie2'):
                print('Distribute per amplicon: ', suffix)
                amplicon_reads = dict()
                print(df_sam.shape)
                # TODO: Add progress bar
                for i,(index, row) in enumerate(df_sam.iterrows()):
                    if program == 'bwa':
                        if i < 2:
                            continue
                    if row[5] == '*':
                        continue
                    start, end = int(row[3]), int(row[3])+len(row[9])
                    already_in = False
                    for amplicon, pair in enumerate(amplicons_limits):    
                        next_amplicon = None
                        if amplicon < (len(amplicons_limits)-1):
                            next_amplicon = amplicons_limits[amplicon + 1]
                        if pair[0] <= start and start <= pair[1]:#and pair[1] >= end: 
                            if amplicon not in amplicon_reads.keys():
                                amplicon_reads[amplicon] = []
                            amplicon_reads[amplicon].append(int(row[0]))
                            already_in = True
                            # We avoid one read to join two different amplicons
                            # break
                # Write reads
                for key, val in amplicon_reads.items():
                    # Create the folders structure
                    tmp_dir = output_dir_tmp+str(key)
                    Utils.mkdir(tmp_dir)
                    tmp_reads_dir = tmp_dir+'/aligned_reads/'
                    Utils.mkdir(tmp_reads_dir)
                    tmp_reads_file = tmp_dir+'/aligned_reads/read'+suffix
                    BioUtils.write_single_end_by_place(reads, val, tmp_reads_file)
                    if key not in amplicon_file_reads.keys():
                        amplicon_file_reads[key] = []
                    amplicon_file_reads[key].append(tmp_reads_file)
            def __launch_first_step_viquf(aligned_reads_dir, tmp_dir_ori):
                tmp_dir = tmp_dir_ori+'/tmp'
                # Arreglar esta ruta relativa
                exe = 'python'
                cmd = [exe,'/home/bfreire/Gatb-trial/scripts/testBcalm.py',aligned_reads_dir,str(KMER_SIZE),'10','ngs','--no-correct','--no-join','--no-meta','--no-filter',tmp_dir]
                print('Launching first step of viquf')
                Utils.executecmd(cmd)
                graph_path, unitigs_path = tmp_dir+'/unitigs.scaled.graph', tmp_dir+'/unitigs.unitigs.fa'
                pair_end_reads = tmp_dir_ori+'/aligned_reads/'
                output_path, append_path = tmp_dir_ori+'/unitigs_ViQUF', tmp_dir_ori+'/tmpOwnlatest/append.fasta'
                cmd = ['./bin/output.out',tmp_dir,str(KMER_SIZE),graph_path, unitigs_path,output_path,append_path,pair_end_reads,'--debug','--virus']
                print('Executing full ViQUF')
                Utils.executecmd(cmd)
                return tmp_dir

            def __get_max_abundance(amplicon_graph_file, output_dir):
                # We want the depth per base to be able to normalize
                # TODO: In real cases we should remove low frequency cases before to avoid fake information to be part of the MEAN
                suma_abundance, suma_length, max_abundance = 0, 0, 0
                with open(amplicon_graph_file, 'r+') as f_read:
                    for i, line in enumerate(f_read.readlines()):
                        if i == 0: 
                            continue
                        line_split = line.strip().split(' ')
                        suma_abundance += float(line_split[-2])*float(line_split[-1])
                        suma_length += float(line_split[-1])
                        max_abundance = float(line_split[-2]) if float(line_split[-2]) > max_abundance else max_abundance
                if suma_abundance != 0:
                    ab_mean_base = (suma_abundance / suma_length)
                else:
                    ab_mean_base = 1
                output_file = output_dir+'/mean_amp_ab.txt'
                with open(output_file, 'w+') as f_write:
                    f_write.write(str(ab_mean_base))
                output_file = output_dir+'/max_amp_ab.txt'
                with open(output_file, 'w+') as f_write:
                    f_write.write(str(max_abundance))
                return output_file, max_abundance#ab_mean_base
            ref_index_dir = output_dir_tmp+'/ref_index/'
            index_prefix = ref_index_dir+'bowtie2_index'
            Utils.mkdir(ref_index_dir)
            # Reference index build
            __build_index(reference_file, index_prefix, program = 'bwa')
            amplicon_boundaries = [(row['start'], row['end']) for index, row in df_amplicons.iterrows()]
            # Reads alignment
            df_left, df_right = __align_reads_complete(reads,index_prefix, sam_file = ref_index_dir+'/sam_file', program = 'bwa')
            # Distribute reads by amplicons
            amplicon_file_reads = dict()
            # amplicon_file_reads = {key:(output_dir_tmp+str(key)+'/aligned_reads/read_1.fasta',output_dir_tmp+str(key)+'/aligned_reads/read_2.fasta') 
            #    for key, (index,row) in enumerate(df_amplicons.iterrows())}
            __distribute_per_amplicon(reads[0], df_left, amplicon_boundaries, amplicon_file_reads,suffix = '_1.fasta', program='bwa')
            __distribute_per_amplicon(reads[1], df_right, amplicon_boundaries, amplicon_file_reads,suffix = '_2.fasta', program = 'bwa')
            max_ab_amp, rate_ab_amplicon = dict(), dict()
            for i,(index, row) in enumerate(df_amplicons.iterrows()):
                if i >= MAX_AMPLICON:
                    break
                tmp_dir = output_dir_tmp+str(i)
                Utils.mkdir(tmp_dir)
                tmp_file = tmp_dir+'/amplicon_ref.fa'
                with open(tmp_file, 'w+') as f_write:
                    f_write.write('>amplicon_ref_'+str(i)+'\n')
                    f_write.write(reference[row['start']:row['end']])
                tmp_index_dir = tmp_dir+'/index'
                index_prefix = tmp_index_dir+'/bowtie2_index'
                Utils.mkdir(tmp_index_dir)
                __build_index(tmp_file,index_prefix)
                # Align the reads
                tmp_reads_dir = tmp_dir+'/aligned_reads/'
                aligned_reads = amplicon_file_reads[i]
                print('Reads written: ', aligned_reads)
                # Launch bcalm for general purposes
                tmp_dir_viquf = __launch_first_step_viquf(tmp_reads_dir, tmp_dir)
                tmp_dir_viquf = tmp_dir+'/tmp'
                # Get max abundance read per amplicon
                suffix, suffix_unitigs = '/unitigs.graph', '/unitigs.unitigs.fa'
                abundance_report_file, max_ab = __get_max_abundance(tmp_dir_viquf+suffix,tmp_dir_viquf)
                max_ab_amp[i] = max_ab
                print('Max abundance reported: ', abundance_report_file,' With max abundance: ', max_ab)
            max_val = 0
            for key, val in max_ab_amp.items():
                if val > max_val:
                    max_val = val
            rate_ab_amplicon = {key:max_val/val for key, val in max_ab_amp.items()}
            print('Maximal abundance: ', max_val)
            # print('Mean abundances amplicon: ', max_ab_amp)
            # print('Rate abundance amplicon: ',rate_ab_amplicon)
            # Readjust frequencies per unitigs and amplicons and reverse complement unitigs
            graphs_files, unitigs_files_complete = [], []
            for key,val in rate_ab_amplicon.items():
                tmp_dir = output_dir_tmp+str(key)+'/tmp/'
                file_read, file_scaled = tmp_dir+'/unitigs.graph', tmp_dir+'/unitigs.scaled.graph'
                graphs_files.append(file_scaled)
                with open(file_read, 'r+') as f_read, open(file_scaled, 'w+') as f_write:
                    for i, line in enumerate(f_read.readlines()):
                        if i == 0:
                            continue
                        line_split = line.strip().split(' ')
                        line_split[len(line_split)-2] = str(int(float(line_split[len(line_split)-2])*val))
                        line_write = ' '.join(line_split)+'\n'
                        f_write.write(line_write)
                # Reverse complement unitigs files
                unitig_file, unitig_file_rc = tmp_dir+'/unitigs.unitigs.fa',tmp_dir+'/unitigs.unitigs.rc.fa'
                new_seqs = BioUtils.reverse_complement_from_file(unitig_file)
                BioUtils.export_fasta_from_dict(new_seqs, unitig_file_rc)
                unitigs_files_complete.append(unitig_file_rc)
            return graphs_files, unitigs_files_complete
        def _join_amplicons(graphs_files, unitigs_files):
            ALIGNMENT_SLACK = 0.015
            print('Joining par_wise amplicons')
            print('Number of amplicons: ', len(unitigs_files))
            new_alignment_files = []
            # Nodes deactivated because of contained information
            nodes_deactivated = {i:[] for i in range(len(unitigs_files))}
            for i in range(len(unitigs_files)-1):
                unitigs_file_left, graph_file_left = unitigs_files[i], graphs_files[i]
                unitigs_file_right, graph_file_right = unitigs_files[i+1], graphs_files[i+1]
                zero_out_degree, zero_in_degree = [], []
                with open(graph_file_left, 'r+') as f:
                    for line in f.readlines():
                        line_split = line.strip().split(' ') 
                        if len(line_split) == 3:
                            zero_out_degree.append(line_split[0])
                nodes,nodes_in = [],dict()
                with open(graph_file_right,'r+') as f:
                    for line in f.readlines():
                        line_split = line.strip().split(' ')
                        nodes.append(line_split[0])
                        l_line = len(line_split)
                        for j in range(1,l_line-2):
                            nodes_in[line_split[j]] = False
                    # We get all nodes from the next amplicon
                    for node in nodes:
                        #if node not in nodes_in.keys():
                        zero_in_degree.append(node)
                u_l_dict, u_r_dict = BioUtils.get_sequences_faf(unitigs_file_left),BioUtils.get_sequences_faf(unitigs_file_right)
                '''if i == 15:
                    print('Amplicon: ',i)
                    print(zero_out_degree)
                    print(zero_in_degree)'''
                new_alignments = []
                for node_left in zero_out_degree:
                    seq_1 = u_l_dict[node_left]
                    for node_right in zero_in_degree:
                        seq_2 = u_r_dict[node_right]
                        alignment = BioUtils.align_two_seqs(seq_1, seq_2)
                        '''if i == 10 and node_right == '60':
                            print('Alignment: ', alignment)
                            if '-' in alignment.seqB[int(len(alignment.seqB)-alignment.score):]:
                                print('Contained alignment: ', alignment)'''
                        if alignment.score > (KMER_SIZE + 2*BioUtils.GAP_PENALIZATION):
                            # We remove proper alignments which end with '-' in the right seq (it means it is contained)
                            if '-' in alignment.seqB[int(len(alignment.seqB)-alignment.score):] and len(seq_2) <= len(seq_1):
                                nodes_deactivated[i+1].append(int(node_right))
                                continue
                            new_alignments.append((node_left, node_right))
                new_alignment_file = '/'.join(unitigs_file_left.split('/')[:-1]+['alignments.txt'])
                with open(new_alignment_file,'w+') as f_write:
                    for (left, right) in new_alignments:
                        f_write.write(left+' '+right+'\n')
                new_alignment_files.append(new_alignment_file)
            return graphs_files, unitigs_files, new_alignment_files, nodes_deactivated
        def _build_complete_graph(graphs_files, unitigs_files, extra_alignments_files, nodes_deactivated):
            return AmpliconsGraph(graphs_files, extra_alignments_files, nodes_deactivated)
        # reads_dir, reference, amplicons_file, type, filter
        reads_dir = sys.argv[1]
        files = Utils.get_files_recursive(reads_dir, ['fasta'])
        files.sort()
        assert len(files) == 2
        reads_left, reads_right = files[0], files[1]
        reference_file, amplicons_file = sys.argv[2], sys.argv[3]
        filtering = (sys.argv[5] == '--no-filter')
        print('*************** Summary *********************')
        print('Amplicons processing with reads: ', reads_left,' ', reads_right)
        print('Reference file: ', reference_file)
        print('Amplicons file: ', amplicons_file)
        print('Filter: ', filtering)
        print('*********************************************')
        # Get reference and amplicons
        reference = BioUtils.readReference(reference_file)
        amplicons = BioUtils.getAmplicons(amplicons_file)
        # Store amplicons tmp
        dir_tmp_amplicons = 'tmp_amplicons/'
        Utils.mkdir(dir_tmp_amplicons)
        # Join paired_end files in tmp dir
        Utils.append_files([reads_left, reads_right], dir_tmp_amplicons+'append.fasta')
        graphs_files, unitigs_files = _create_amplicons_segments(amplicons, reference, (reads_left,reads_right), reference_file=reference_file)
        # Join amplicons between them (using only fw strain)
        graphs_files, unitigs_files, extra_alignments_files, nodes_deactivated = _join_amplicons(graphs_files, unitigs_files)
        # Build graph amplicons
        graph_amplicons = _build_complete_graph(graphs_files, unitigs_files, extra_alignments_files, nodes_deactivated)
        # Get final graph file with unitigs
        graph_file, final_unitigs = dir_tmp_amplicons+'amplicons.graph', dir_tmp_amplicons+'unitigs.unitigs.fa'
        graph_amplicons.exportGraph(graph_file,final_unitigs, unitigs_files)
        sys.exit(1)

    #BioUtils.identicalClustering(__preprocess2([pathIn[1],rG.getOutFile()]))
