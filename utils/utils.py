import os, subprocess
import pandas as pd
import numpy as np
import altair as alt
import shutil
from shutil import copyfile

class Utils:
    @staticmethod
    def cpfile(source, dest):
        copyfile(source, dest)
        return dest

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
        print('Removing: ',file)
        if Utils.exists(file):
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
    def get_files_recursive(path, extension = ['fq','fa','fasta','fastq'], threshold = None):
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
    def append_files_bash(list_files, output_file):
        outDir = '/'.join(output_file.split('/')[:-1])
        if not Utils.exists(outDir):
            Utils.mkdir(outDir)
        cmd = ['cat']+[t for t in list_files]
        Utils.executecmd(cmd, output_file)
        return output_file

    @staticmethod
    def executecmd(args, out = None):
        '''
        :param path: LIST of args
        :return:
        '''
        print(' '.join(args))
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
from Bio.Seq import Seq
from Bio.File import as_handle
from Bio import pairwise2
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment
from Bio.SearchIO._utils import get_processor
import pandas as pd

class BioUtils:
    GAP_PENALIZATION = -10
    @staticmethod
    def reverse_complement(chain):
        return str(Seq(chain).reverse_complement())

    @staticmethod
    def reverse_complement_from_file(file):
        seqs = SeqIO.parse(file, format = 'fasta')
        fw_seqs, rc_seqs = [], []
        for record in seqs:
            fw_seqs.append(str(record.seq)+'\n')
            rc_seqs.append(str(record.seq.reverse_complement())+'\n')
        fw_seqs += rc_seqs
        seqs = dict()
        for i, seq in enumerate(fw_seqs):
            seqs[str(i)] = seq
        return seqs

    @staticmethod
    def align_two_seqs(seq_1, seq_2):
        # -1 missmatch +1 match 0 indels
        return pairwise2.align.globalms(seq_1, seq_2, 1,-10,BioUtils.GAP_PENALIZATION,0)[0]

    @staticmethod
    def identicalClustering(file, program = 'mmseqs', args = ['tmp/clusters','tmp/', '-c','1.0','--min-seq-id',
                                                              '1.0','-v','3','--seq-id-mode','0','--cov-mode','1','--remove-tmp-files'
                                                              ,'--dbtype','1']):
        mode = 'easy-cluster'
        cmd = [program,mode,file]+args
        print('Cmd: ', cmd)
        Utils.executecmd(cmd)

    @staticmethod
    def build_align_index(ref_file, db_prefix, program = 'bowtie2'):
        if program == 'bowtie2':
            exe = 'bowtie2-build'
            cmd = [exe,ref_file, db_prefix,'-q']
        if program == 'bwa':
            exe = 'bwa'
            cmd = [exe,'index', ref_file, '-p', db_prefix]
        print('Cmd: ',' '.join(cmd))
        Utils.executecmd(cmd)
    
    @staticmethod
    def align_reads(reads_left, reads_right, index_prefix, output_reads_path, sam_file, format = 'fasta', program = 'bowtie2'):
        if program == 'bowtie2':
            exe = 'bowtie2'
            if format == 'fasta':
                cmd = [exe,'-x',index_prefix,'-1',reads_left,'-2',reads_right,'--al-conc',output_reads_path,'--quiet','-f', '-S', sam_file]
            else:
                cmd = [exe,'-x',index_prefix,'-1',reads_left,'-2',reads_right,'--al-conc',output_reads_path,'--quiet', '.-S', sam_file]
        print('Cmd: ',' '.join(cmd))
        Utils.executecmd(cmd)

    @staticmethod
    def align_single_end_reads(reads, index_prefix, sam_file, format = 'fasta', program = 'bowtie2', permissive = False):
        if program == 'bowtie2':
            exe = 'bowtie2'
            # Soft parameters: --local -D 20 -R 3 -L 3 -N 1 -p 8 --gbar 1 --mp 3
            if format == 'fasta':
                cmd = [exe,'-x',index_prefix,'-U',reads,'--quiet','-f', '-S', sam_file,'--no-head']
                soft_params = []
                if permissive:
                    soft_params = ['--local','-D','20','-R','3','-L','3','-N','1','-p','8','--gbar','1','--mp','3']
                cmd += soft_params
            else:
                cmd = [exe,'-x',index_prefix,'-U',reads,'--quiet', '-S', sam_file,'--no-head']
        if program == 'bwa':
            # bwa mem testing/hcv_ref.fasta tmp/read_1.fasta tmp/read_2.fasta
            exe = 'bwa'
            if format == 'fasta':
                cmd = [exe,'mem',index_prefix,reads,'-o', sam_file,'-v','1']
        print('Cmd: ',' '.join(cmd))
        Utils.executecmd(cmd)

    # Revisar para caso de lecturas no simuladas
    @staticmethod
    def write_paired_end_reads_by_place(reads, pos, output_files):
        reads_left, reads_right = reads
        output_left, output_right = output_files
        with open(reads_left,'r+') as f_left, open(output_left,'w+') as f_left_write, open(reads_right, 'r+') as f_right, open(output_right,'w+') as f_right_write:
            lines_left, lines_right = f_left.readlines(), f_right.readlines()
            for i in pos:
                f_left_write.write('>'+str(i)+' 1\n')
                f_left_write.write(lines_left[2*i+1])
                f_right_write.write('>'+str(i)+' 2\n')
                f_right_write.write(lines_right[2*i+1]) 
    @staticmethod
    def write_single_end_by_place(reads, pos, output_files):
        reads_left = reads
        output_left = output_files
        with open(reads_left,'r+') as f_left, open(output_left,'w+') as f_left_write:
            lines_left = f_left.readlines()
            for i in pos:
                f_left_write.write('>'+str(i)+' 1\n')
                f_left_write.write(lines_left[2*i+1])

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
    def fastqtofasta(in_fastq, out_fasta):
        cmd = ['sed', '-n' ,'1~4s/^@/>/p;2~4p', in_fastq]
        Utils.executecmd(cmd, out_fasta)
        return out_fasta

    @staticmethod
    def renameFastqSeqs(files, outDir,format = 'fasta'):
        newDirs = []
        for i,f in enumerate(files):
            writeFile = outDir+str(i)+'.fasta'
            with open(f,'r') as fRead, open(writeFile,'w+') as fWrite:
                for j,line in enumerate(fRead.readlines()):
                    if j%4 == 0:
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

    @staticmethod
    def readFasta(file, format = 'fasta'):
        ids, seqs = [], []
        for record in SeqIO.parse(file, format):
            seqs.append(record.seq)
            ids.append(record.id)
        return ids,seqs
    
    @staticmethod
    def readReference(file):
        reference = ''
        for record in SeqIO.parse(file, 'fasta'):
            reference = str(record.seq)
        return reference

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
    def getAmplicons(file, format = 'csv'):
        if format == 'csv':
            df_amplicons = pd.read_csv(file, sep = ',')
        return df_amplicons

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
                r = r.strip()
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

    @staticmethod
    def splitFqFile(file, outputPattern, sep = '/'):
        left, right = outputPattern+'_1.fq',outputPattern+'_2.fq'
        l_or_r, p = False, False
        with open(file, 'r') as f, open(left,'a+') as f_left, open(right,'a+') as f_right:
            for i,line in enumerate(f.readlines()):
                if p:
                    if l_or_r:
                        f_left.write(line)
                    else:
                        f_right.write(line)
                if line[0] == '@':
                    l_split = line.rstrip('\n').split('/')
                    if len(l_split) > 1:
                        p = True
                        l_or_r = (l_split[1] == '1')
                else:
                    p = False
        return f_left,f_right
