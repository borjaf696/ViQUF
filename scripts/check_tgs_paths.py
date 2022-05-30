# Borja :)
from platform import node
import sys, os
sys.path.append('')
from utils.utils import BioUtils
def _process_gfa(file):
    print('Process gfa')
    seq_dir = dict()
    with open(file, 'r+') as f:
        lines = f.readlines()
        num_seqs= 0
        for line in lines:
            line_split = line.strip().split('\t')
            if line_split[0] == 'S':
                if num_seqs < int(line_split[1]):
                    num_seqs = int(line_split[1])
                    
        print('Number of seqs: ',num_seqs)
        total_seqs = 2*(num_seqs+1)
        offset = num_seqs+1
        print('Total number (rev): ',total_seqs)
        print('Reverse offset: ',offset)
        graph_struct = [[] for _ in range(0,total_seqs)]
        for line in lines:
            line_split = line.strip().split('\t')
            if line_split[0] == 'L':
                src, target = int(line_split[1]),int(line_split[3])
                src_dir, target_dir = (line_split[2] == '+'), (line_split[4] == '+')
                src += (0 if src_dir else offset)
                target += (0 if target_dir else offset)
                graph_struct[src].append(target)
            if line_split[0] == 'S':
                fw_seq, rev_seq = line_split[2], BioUtils.reverse_complement(line_split[2])
                fw_id, rev_id = int(line_split[1]), int(line_split[1]) + offset
                seq_dir[fw_id] = fw_seq
                seq_dir[rev_id] = rev_seq
        return graph_struct, seq_dir

def _process_paths(file):
    paths = []
    with open(file, 'r+') as f_read:
        lines = f_read.readlines()
        for line in lines:
            two_sides = line.strip().split('[')
            flow = int(two_sides[0])
            two_sides[1] = two_sides[1][0:len(two_sides[1])-2].replace("'",'')
            partial_path = two_sides[1].split(', ')
            partial_path = list(map(lambda x: int(x.strip()),partial_path))
            paths.append((flow, partial_path))
    return paths

def traverse_paths(graph, graph_in, seqs, paths, kmerSize = 50):
    def __get_unary_path(node):
        traversed_path = []
        in_degree, out_degree = len(graph_in[node]), len(graph[node])
        while (in_degree == 1) and (out_degree == 1):
            traversed_path.append(node)
            node = graph[node][0]
            in_degree, out_degree = len(graph_in[node]), len(graph[node])
        return traversed_path if len(traversed_path) > 0 else [node]
    contigs = []
    for p_s in paths:
        local_contig = []
        for i in range(len(p_s[1])):
            p = p_s[1][i]
            if p == -1:
                continue
            if p == len(graph):
                break
            local_contig += __get_unary_path(p)
        seq_local_contigs = ''
        for n_contig, i in enumerate(local_contig):
            if n_contig == i:
                seq_local_contigs += seqs[i]
            else:
                seq_local_contigs += seqs[i][kmerSize-1:]
        contigs.append((p_s[0],seq_local_contigs))
    return contigs

def reverse_graph(graph):
    graph_in = [[] for _ in graph]
    for node,l_neighs in enumerate(graph):
        for neigh in l_neighs:
            graph_in[neigh].append(node)
    return graph_in

if __name__ == '__main__':
    # Gfa file with unitigs information
    gfa = 'tmp/unitigs'
    path_file,kmerSize = sys.argv[1], int(sys.argv[2])
    graph, seqs = _process_gfa(gfa)
    graph_in = reverse_graph(graph)
    paths = _process_paths(path_file)
    contigs = traverse_paths(graph, graph_in, seqs, paths, kmerSize)
    with open('tmp/contigs_ilp.fa','w+') as contigs_write:
        for i,contig in enumerate(contigs):
            contigs_write.write('>'+str(i)+' freq: '+str(contig[0])+'\n')
            contigs_write.write(contig[1]+'\n')