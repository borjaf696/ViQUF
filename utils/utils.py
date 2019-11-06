import os, subprocess
import altair as alt
import pandas as pd
import numpy as np
import shutil
from shutil import copyfile

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
    def exists(d):
        return os.path.exists(d)

    @staticmethod
    def remove_file(file):
        os.remove(file)

    @staticmethod
    def remove_dir(folder):
        shutil.rmtree(folder, ignore_errors=True)

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
    def get_files_recursive(path, extension = ['fq','fa'], threshold = None):
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
                if threshold is not None and len(results) > threshold:
                    return results
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
        outDir = '/'.join(output_file.split('/')[:-1])
        if not Utils.exists(outDir):
            Utils.mkdir(outDir)
        with open(output_file,'w+') as out_file:
            for f in list_files:
                with open(f,'r') as f_read:
                    for line in f_read.readlines():
                        out_file.write(line)
        return output_file

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
        if format == 'csv':
            with open(name, 'w+') as f:
                for key, val in dictionary.items():
                    f.write(str(key))
                    for k,v in val.items():
                        f.write(','+str(v))
                    f.write('\n')
        if format == 'fasta':
            with open(name, 'w+') as f:
                for key, val in dictionary.items():
                    f.write(str(key))
                    f.write(str(val))

    @staticmethod
    def writeMatrix(m, file = None):
        assert file is not None

        rows,cols = np.shape(m)
        with open(file, 'w+') as f:
            for i in range(rows):
                for j in range(cols):
                    f.write(str(m[i,j])+' ')
                f.write('\n')


class StatsReport:
    @staticmethod
    def show_histogram(histogram, path, y_label = 'y', x_label = 'x'):
        d = {y_label:histogram,x_label:range(len(histogram))}
        source = pd.DataFrame(d)
        chart = alt.Chart(source).mark_bar().encode(alt.X(x_label, bin = True),y=y_label)
        chart_2 = alt.Chart(source).mark_bar().encode(alt.X(x_label,  bin=alt.Bin(maxbins=100)), y = y_label)
        chart.save(path)
        chart_2.save(path+'_100_bins.html')

    @staticmethod
    def dotPlot(df, X, Y, nameFile, color = None):
        chart = alt.Chart(df)
        chart = chart.mark_point().encode(x=X, y = Y) if color is None else chart.mark_point().encode(x=X,y=Y, color=color)
        chart.save(nameFile+'.html')

    @staticmethod
    def heatmap(df, X, Y, m, nameFile):
        chart = alt.Chart(df).mark_rect().encode(x=X+':O',y=Y+':O',color=m+':Q')
        chart.save(nameFile+'.html')

    @staticmethod
    def iterative_levenshtein(s, t, costs=(1, 1, 1)):
        rows = len(s) + 1
        cols = len(t) + 1
        deletes, inserts, substitutes = costs

        dist = [[0 for x in range(cols)] for x in range(rows)]
        for row in range(1, rows):
            dist[row][0] = row * deletes
        for col in range(1, cols):
            dist[0][col] = col * inserts

        for col in range(1, cols):
            for row in range(1, rows):
                if s[row - 1] == t[col - 1]:
                    cost = 0
                else:
                    cost = substitutes
                dist[row][col] = min(dist[row - 1][col] + deletes,
                                     dist[row][col - 1] + inserts,
                                     dist[row - 1][col - 1] + cost)
        '''for r in range(rows):
            print(dist[r])'''
        return dist[rows-1][cols-1]

from Bio import SearchIO, SeqIO
from Bio.File import as_handle
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment
from Bio.SearchIO._utils import get_processor

class BioUtils:
    @staticmethod
    def identicalClustering(file, program = 'mmseqs', args = ['tmp/clusters','tmp/', '-c','1.0','--min-seq-id',
                                                              '1.0','-v','3','--seq-id-mode','0','--cov-mode','1','--remove-tmp-files'
                                                              ,'--dbtype','1']):
        mode = 'easy-cluster'
        cmd = [program,mode,file]+args
        print('Cmd: ', cmd)
        Utils.executecmd(cmd)

    @staticmethod
    def makedb(infile, name, type = 'blast', dbtype = 'prot'):
        cmd = ''
        if type == 'blast':
            exe = 'makeblastdb'
            cmd = [exe, '-in',infile, '-out', name, '-dbtype', dbtype]
        print('Cmd: ', cmd)
        Utils.executecmd(cmd)

    @staticmethod
    def mafftAlign(inFile, outmafft):
        mafft_exe = 'mafft'
        mafft_cline = MafftCommandline(mafft_exe, input=inFile, clustalout=False, thread=1)
        # Perform alingment with MAFFT
        stdout, stderr = mafft_cline()
        from io import StringIO
        from Bio import AlignIO
        align = AlignIO.read(StringIO(stdout), 'fasta')
        # Save alignment into a File
        import Bio.AlignIO
        AlignIO.write(align, outmafft, format='fasta')
        return outmafft

    @staticmethod
    def getPHMM(profilePHMM, query, outputFile = 'hmmhits2.out',e_value = 10e-30):
        #!hmmsearch - -domtblout hmmhits_domain.out mnmE_589_mafft_ref405.fastas.hmm proteins.faa > hmmhits2.out
        exe = 'hmmsearch'
        cmd = [exe, '--domtblout', 'tmp_hmmhits_domain.out', profilePHMM, query]
        print('Cmd: ',cmd, 'Output file: ', outputFile)
        Utils.executecmd(cmd, outputFile)
        return outputFile

    @staticmethod
    def parseOutput(ids = None, output_file = 'hmmhits2.out', format = 'hmmer3-text'):
        #hmmer_qresult = SearchIO.read(output_file, format)
        assert ids is not None
        results = dict()
        for id in ids:
            results[id] = {'hmm':None, 'evalue':100}
        hmmer_results_big = SearchIO.parse(output_file, format)
        for r in hmmer_results_big:
            for partial in r:
                seq, hmm = (partial.hsps[0]).hit_id, (partial.hsps[0]).query_id
                evalue = (partial.hsps[0]).evalue
                if results[seq]['evalue'] > evalue:
                    results[seq] = {'hmm':hmm, 'evalue':evalue}
        return results

    @staticmethod
    def renameFastqSeqs(files, outDir,format = 'fasta'):
        newDirs = []
        for i,f in enumerate(files):
            writeFile = outDir+str(i)+'.fasta'
            with open(f,'r') as fRead, open(writeFile,'w+') as fWrite:
                for j,line in enumerate(fRead.readlines()):
                    if line[0] == '@':
                        fWrite.write('>'+str(j/4)+'.'+str(i)+'\n')
                        write = True
                    elif write:
                        fWrite.write(line)
                        write = False
            newDirs.append(writeFile)
        return newDirs

    @staticmethod
    def fastToFm(file, fileOut):
        fileOut_2, curPos = 'placements', 0
        with open(fileOut, 'w+') as fWrite, open(fileOut_2,'w+') as fWrite2:
            for record in SeqIO.parse(file, 'fasta'):
                fWrite.write(str(record.seq)+'$')
                curPos += len(str(record.seq))+1
                fWrite2.write(str(curPos-1)+'\n')
        with  open(fileOut,'r') as fRead, open(fileOut_2,'r') as fRead2:
            line1 = fRead.readlines()[0]
            for line in fRead2.readlines():
                print(line1[int(line)])

    @staticmethod
    def readFasta(file, format = 'fasta'):
        ids, seqs = [], []
        for record in SeqIO.parse(file, format):
            seqs.append(record.seq)
            ids.append(record.id)
        return ids,seqs

    @staticmethod
    def getSeqsByIds(file, ids = None):
        assert ids is not None
        results = dict()
        with open(file, 'r') as f:
            lines = f.readlines()
            for i,r in enumerate(lines):
                if '>' in r:
                    id = r.split()[0][1:]
                    if id in ids:
                        results[id] = lines[i+1]
        return results

    @staticmethod
    def getIdsFromFa(path):
        ids = []
        with open(path, 'r') as f:
            for r in f.readlines():
                if '>' in r:
                    ids.append(r.split()[0][1:])
        return ids

    @staticmethod
    def get_sequences_faf(path):
        '''
        Retrieves the sequence from a fasta file format
        :param path:
        :return:
        '''
        sequences = dict()
        with open(path, 'r') as f:
            for r in f.readlines():
                if '>' in r:
                    id = r[1:len(r)]
                else:
                    sequences[id] = r
        return sequences

    @staticmethod
    def export_fasta_from_dict(d,output):
        with open(output,'w+') as f:
            for (key,val) in d.items():
                f.write('>'+key+'\n'+val)

    @staticmethod
    def getAllPairedEndDirs(d):
        path_files = []
        list_dirs = Utils.get_dirs(d)
        for d in list_dirs:
            list_dirs += Utils.get_dirs(d)
            if len(Utils.get_files(d)) == 2:
                path_files.append(d)
        return path_files

    @staticmethod
    def getAllPairedEndFiles(d):
        path_files = []
        list_dirs = Utils.get_dirs(d)
        for d in list_dirs:
            list_dirs += Utils.get_dirs(d)
            if len(Utils.get_files(d)) == 2:
                path_files.append(tuple(Utils.get_files(d)))
        return path_files

    @staticmethod
    def appendProteinFiles(list_files, output_file, output_file_exclude, exclude = 'hypothetical'):
        avoid_next = False
        with open(output_file, 'w+') as out_file, open(output_file_exclude, 'w+') as out_file_exclude:
            for f in list_files:
                with open(f, 'r') as f_read:
                    for line in f_read.readlines():
                        if exclude in line:
                            avoid_next = True
                            out_file_exclude.write(line)
                        elif avoid_next:
                            if '>' in line:
                                out_file.write(line)
                                avoid_next = False
                            else:
                                out_file_exclude.write(line)
                        else:
                            out_file.write(line)

    @staticmethod
    def annotateGenome(genome, output_dir = None, exe = 'prodigal'):
        outputProteins = None
        if exe == 'prodigal':
            assert output_dir is not None
            output, outputProteins = 'annotation.txt', output_dir + '/proteins.faa'
            cmd = [exe, '-i', genome, '-o', output, '-a', outputProteins]
            print('Annotation Cmd: ', cmd)
            Utils.executecmd(cmd)
        return outputProteins

    @staticmethod
    def filterFile(fileIn, fileOut, iConf):
        inf, sup = iConf[0],iConf[1]
        with open(fileOut, 'w+') as fw, open(fileIn, 'r') as fr:
            remove = False
            lines = fr.readlines()
            i = 0
            while i < len(lines):
                justOne = False
                line = lines[i]
                if line[0] == '>' and lines[i+1][0] == '>' and (len(lines[i+2]) < inf or len(lines[i+2] > sup)):
                    remove = True
                elif line[0] == '>' and lines[i+1][0] == '>':
                    remove = False
                elif line[0] == '>' and (len(lines[i+1]) < inf or len(lines[i+1] > sup)):
                    justOne = True
                elif line[0] != '>':
                    justOne = True
                if remove == False:
                    if justOne == False:
                        if lines[i+1][0] == '>':
                            fw.write(line+line[i+1]+line[i+2])
                            i += 3
                        else:
                            fw.write(line+line[i+1])
                            i += 2
                else:
                    i += 1
        Utils.remove_file(fileIn)
        Utils.cpfile(fileOut, fileIn)
        Utils.remove_file(fileOut)
