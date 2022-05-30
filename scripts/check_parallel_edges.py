# Borja :)
import sys, os
from tabnanny import check 

def check_file(graph_file):
    check_pairs = dict()
    with open(graph_file, 'r+') as f_read:
        for i, line in enumerate(f_read.readlines()):
            if i < 2:
                continue
            line_split = line.strip().split(' ')
            src, target = int(line_split[0]), int(line_split[1])
            if (src, target) in check_pairs.keys():
                print('Parallel edges found!')
                sys.exit(1)
    print('Not parallel edges')
    
if __name__ == '__main__':
    print(sys.argv)
    graph_file = sys.argv[1]
    check_file(graph_file)
