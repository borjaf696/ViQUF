from utils import *
import sys, os
sys.path.append('')

tmp_folder = 'tmp'
append_file = 'tmp/append.fq'
output_folder = 'output'

if __name__ == '__main__':
    folder, output_file_pattern = sys.argv[1], sys.argv[2]
    files = Utils.get_files(folder)
    Utils.mkdir(tmp_folder)
    Utils.mkdir(output_folder)
    #Utils.append_files(files, append_file)
    for f in files:
        print('File: ',f)
        print('Files: ',BioUtils.splitFqFile(append_file, output_folder+'/'+output_file_pattern, sep = '/'))

