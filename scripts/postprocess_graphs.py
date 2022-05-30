#Borja :)
import sys,os

def read_graph(graph_type, graph_file):
    print('Tipo de grafo: ', graph_type)
    print('Graph name: ',graph_file)
    corrected_graph = graph_file + '.corrected'
    endpoints = dict()
    with open(graph_file, 'r+') as f_read, open(corrected_graph,'w+') as f_write:
        lines = f_read.readlines()
        for i,line in enumerate(lines):
            if i < 2:
                f_write.write(line)
                continue
            line_split = line.strip().split(' ')
            line_split = list(map(lambda x: int(x),line_split))
            key = (line_split[0],line_split[1])
            if key not in endpoints:
                endpoints[key] = line_split[2]
            else:
                endpoints[key] += line_split[2]
        for key,val in endpoints.items():
            f_write.write(str(key[0])+' '+str(key[1])+' '+str(val)+'\n')
        

if __name__ == '__main__':
    print('Params: [std,inexact,subpaths] [file_name]')
    graph_type = sys.argv[1]
    graph_file = sys.argv[2]

    read_graph(graph_type, graph_file)