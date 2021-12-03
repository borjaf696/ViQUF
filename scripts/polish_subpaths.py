# Borja :) 
import os, sys

def _polish_subpaths(path):
    def _check_subpaths(subpaths, new_constraint):
        for key,subpath in subpaths.items():
            all_in = True
            for node in new_constraint:
                if node not in subpath:
                    all_in = False
                    break
            if all_in:
                return False
        return True
    subpaths = dict()
    new_path = str(path.split('.')[:-1][0])+'_corrected.graph'
    with open(path, 'r+') as f, open(new_path,'w+') as w:
        start_check = False
        for i,line in enumerate(f.readlines()):
            l = line.strip()
            if 'subpaths' in l:
                start_check = True
                subpath_starts = i
            elif start_check:
                new_subpath = list(map(lambda x:int(x),l.split(' ')[:-1]))
                if _check_subpaths(subpaths, new_subpath):
                    subpaths[i - subpath_starts - 1] = new_subpath
            else:
                w.write(line)
        w.write('# subpaths\n')
        for key, subpath in subpaths.items():
            new_s = list(map(lambda x: str(x), subpath))
            w.write(' '.join(new_s)+' 1.0\n')
if __name__ == '__main__':
    assert len(sys.argv) > 1

    _polish_subpaths(sys.argv[1])