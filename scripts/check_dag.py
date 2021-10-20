## Borja :)
import sys, os

class graph():
    def __init__(self,file):
        self._edges = []
        self.__starts = []
        self.__process_file(file)
    
    def __process_file(self, file):
        with open(file, 'r+') as f_read:
            for i,line in enumerate(f_read.readlines()):
                line_split = line.split(' ')
                if i == 0:
                    num_nodes = int(line_split[0])
                    self._edges = [[] for i in range(num_nodes)]
                    self._balance =[0]*num_nodes
                    continue
                l = len(line_split)
                if l > 3:
                    node = int(line_split[0])
                    neighbors = line_split[1:l-3]
                    for j in neighbors:
                        self._balance[node] -= 1
                        self._balance[int(j)] += 1
                        self._edges[node].append(int(j))
            for i, recount in enumerate(self._balance):
                if recount > 0 and len(self._edges[i]) > 0:
                    #print(recount,' ',self._edges[i],' ', i)
                    self.__starts.append(i)
    # Naive!
    def _traverse_from(self, node, checked, paths_traversed):
        dags = []
        print(checked)
        for neigh in self._edges[node]:
            if neigh not in checked and neigh not in paths_traversed.keys():
                result = self._traverse_from( neigh,checked+[neigh], paths_traversed)
                paths_traversed[neigh] = result
                dags.append(result)
            if neigh in checked:
                return False
        acyclic = True
        for b in dags:
            acyclic &= b
        return acyclic

    def is_dag(self):
        dags = []
        paths_traversed = dict()
        for i in self.__starts:
            checked = [i]
            dags.append(self._traverse_from(i, checked, paths_traversed))
        for b in dags:
            if not b:
                return False
        return True

if __name__ == '__main__':
    g = graph(sys.argv[1])
    print(g.is_dag())