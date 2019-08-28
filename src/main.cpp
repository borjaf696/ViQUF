// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <vector>
#include <chrono>
#include <map>
#include <string>
#include <queue>
#include <stack>
#include <set>
/********************************************************************************/
/********************************************************************************/
#define INF 9999999

#define KMERS 15
#define ABUNDANCE 132
#define CLICK_LIMIT 2
#define DELTA_PATH 600

#define MIN_SIZE_REP 0.0
#define CLICK_RATIO 0.1
#define SIMILARITY_RATIO 1.0
#define IDENTITY_RATIO 1.0
#define IDENTITY_RATIO_UNION 1.0

#define PARENT_SON_LIMIT 0.1

using namespace std;

/*
 * Set operations
 */
template<typename T>
void show_set(const set<T> & set1)
{
    for (auto p: set1)
        cout << " "<<p.str();
    cout << endl;
}

template<typename T>
set<T> sustract(const set<T> & set1, const set<T> & set2)
{
    set<T> s_set;
    for (auto i = set1.begin(); i != set1.end(); ++i)
        if (set2.find(*i) == set2.end())
            s_set.emplace(*i);
    return s_set;
}

template<typename T>
set<T> getIntersection(const set<T>& set_full,const set<T>& nn_pair)
{
    set<T> intersection;
    for (auto i = nn_pair.begin(); i != nn_pair.end(); i++) {
        if (set_full.find(*i) != set_full.end())
            intersection.insert(*i);
    }
    return intersection;
}

template<typename T>
set<T> getUnion(const set<T>& set1, const set<T>& set2)
{
    set<T> result = set1;
    result.insert(set2.begin(),set2.end());
    return result;
}

template<typename T>
bool in(const set<T> & set1, T el)
{
    return (set1.find(el)!=set1.end());
}

template<typename T>
bool isSubset(const set<T> & set1, const set<T> & set2)
{
    for (auto el: set1)
        if (set2.find(el) == set2.end())
            return false;
    return true;
}

template <typename T>
bool isSubset(const set<T> & set1, const set<T> & set2, double percent)
{
    int32_t num_fails = ceil(set1.size()*(1-percent)), fails = 0;
    if (set1.size() > set2.size()+num_fails)
    {
        for (auto el:set2)
            if (set1.find(el) == set1.end())
                if (++fails >= num_fails) return false;
        return true;
    }
    for (auto el: set1)
        if (set2.find(el) == set2.end())
            if (++fails>=num_fails) return false;
    return true;
}

template<typename T>
set<T> isGetSubset(const set<T> & set1, const set<T> & set2)
{
    set<T> result;
    for (auto el: set1) {
        if (set2.find(el) == set2.end())
            return set<T>();
        result.emplace(el);
    }
    return result;
}

template<typename T>
bool isSame(const set<T> & set1, const set<T> & set2, double percent)
{
    size_t max = std::max(set1.size(), set2.size());
    int32_t num_fails = ceil(max*(1-percent)), fails = 0;
    if (set1.empty() && !set2.empty())
        return false;
    if (set2.empty() && !set1.empty())
        return false;
    for (auto el:(max ==set1.size())?set1:set2)
        if ((max==set1.size())?set2.find(el) == set2.end():set1.find(el)==set2.end())
            if (++fails >= num_fails) return false;
    return true;
}

template<typename T>
bool isSame(const set<T> & set1, const set<T> & set2)
{
    if (set1.empty() && set2.empty())
        return false;
    for (auto el: set1)
    {
        if (set2.find(el)==set2.end())
            return false;
    }
    for (auto el:set2)
    {
        if (set1.find(el) == set1.end())
            return false;
    }
    return true;
}

template<typename T>
bool same(const set<T> & set1, const set<T> & set2)
{
    for (auto el: set1)
    {
        if (set2.find(el) == set2.end())
            return false;
    }
    for (auto el: set2)
    {
        if (set2.find(el) == set2.end())
            return false;
    }
    return true;
}

template<typename T>
class SimpleUndirectedGraph{
public:
    typedef size_t vertex_graph;

    SimpleUndirectedGraph(){}
    template<typename G>
    SimpleUndirectedGraph(G iterable)
    {
        for (auto it: iterable)
            addVertex(it);
    }

    void stats()
    {
        for (auto v: _adjacency_list)
        {
            _num_edges += v.second.size();
        }
        _num_edges /= 2;
    }

    size_t addVertex(const T & node)
    {
        _adjacency_list.push_back(pair<T,set<vertex_graph >>(node, set<vertex_graph>()));
        return _num_vertices++;
    }

    void addEdge(vertex_graph src_pos, vertex_graph target_pos)
    {
        _adjacency_list[src_pos].second.emplace(target_pos);
    }

    set<vertex_graph> getNeighbors(const T & node)
    {
        return _adjacency_list[_translateT(node)];
    }

    set<T> getTNeighbors(const T & node)
    {
        set<T> output;
        for (auto n: _adjacency_list[_translateT(node)])
            output.push_back(_adjacency_list[n].first);
        return output;
    }

    set<vertex_graph> getNeighbors(const vertex_graph & node) const
    {
        return _adjacency_list[node].second;
    }

    size_t indegree(const vertex_graph & node) const
    {
        return _adjacency_list[node].second.size();
    }

    size_t outdegree(const vertex_graph & node) const
    {
        return _adjacency_list[node].second.size();
    }

    size_t degree(const vertex_graph & node) const
    {
        return _adjacency_list[node].second.size();
    }

    size_t indegree(const T & node) const
    {
        return _adjacency_list[_translateT(node)].second.size();
    }

    size_t outdegree(const T & node) const
    {
        return _adjacency_list[_translateT(node)].second.size();
    }

    size_t degree(const T & node) const
    {
        return _adjacency_list[_translateT(node)].second.size();
    }

    vector<T> getVertex() const
    {
        vector<T> vertex;
        for (auto v: _adjacency_list)
            vertex.push_back(v.first);
        return vertex;
    }

    vector<size_t> getVertex_Id() const
    {
        vector<size_t> vertex;
        size_t cont = 0;
        for(auto v: _adjacency_list)
            vertex.push_back(cont++);
        return vertex;
    }

    size_t getNumVertices() const
    {
        return _num_vertices;
    }

    size_t getNumEdges() const
    {
        return _num_edges;
    }

    T getNode(vertex_graph src) const
    {
        return _adjacency_list[src].first;
    }
private:
    size_t _translateT(const T & node)
    {
        size_t cont = 0;
        for (auto n:_adjacency_list) {
            if (n.first == node)
                return cont;
            cont++;
        }
    }
    vector<pair<T, set<size_t>>> _adjacency_list;
    size_t _num_edges = 0, _num_vertices = 0;
};

template<typename G, typename Vt>
priority_queue<pair<size_t,vector<Vt>>> findMaxClique(const G & graph, std::set<size_t> & idCliques,
                                                           map<Vt,size_t> & store_map ) {
    if (!graph.getNumEdges())
        return priority_queue<pair<size_t,vector<Vt>>>();
    priority_queue<pair<size_t,vector<Vt>>> endCliques;
    vector<Vt> maxClique;
    vector<Vt> tmpClique;
    set<Vt> already_checked;
    /*
     * Sum -> shift id to the left. Example:
     *      - 0 1 4 -> 1 0 0 1 1 -> 19
     *      - 0 3 2 -> 1 1 0 0 -> 10
     * Both cases sum = 5 but id_different. Problem with cliques larger than 64 Nodes.
     */
    int64_t sum;
    size_t min_flow;

    auto findMaxCliqueWithVertex = [&sum, &already_checked, &min_flow](const Vt vertex, const int maxCliqueSize, const G &graph,
                                                                       map<Vt,size_t> store_map, size_t * score)
    {
        vector<Vt> clique;
        clique.emplace_back(vertex);
        sum |= 1 << vertex;

        set<Vt> candidateNeighbors;

        set<Vt> visited;
        visited.emplace(vertex);

        //Testear
        set<Vt> neighbors = graph.getNeighbors(vertex);
        for (auto n: neighbors)
              candidateNeighbors.emplace(n);

        set<Vt> tmp;

        while (!candidateNeighbors.empty()) {
            const auto highestDegNeighborIt = max_element(candidateNeighbors.begin(), candidateNeighbors.end(), [graph, store_map](const Vt &lhs, const Vt &rhs) {
                if ((store_map.at(lhs)*graph.degree(lhs)) == (store_map.at(rhs) * graph.degree(rhs)))
                    return (lhs > rhs);
                return (store_map.at(lhs) * graph.degree(lhs)) < (store_map.at(rhs) * graph.degree(rhs));
            });

            const auto highestDegVert = *highestDegNeighborIt;
            min_flow = (store_map.at(highestDegVert) < min_flow)?store_map.at(highestDegVert):min_flow;
            clique.emplace_back(highestDegVert);
            (*score) += store_map.at(highestDegVert);
            //already_checked.emplace(highestDegVert);
            /*
             * Questionable
             */
            sum |= 1 << highestDegVert;
            visited.emplace(highestDegVert);

            for (auto adjVertex: graph.getNeighbors(highestDegVert))
            {
                if (candidateNeighbors.find(adjVertex) != candidateNeighbors.end() && visited.find(adjVertex) == visited.end())
                {
                    tmp.insert(adjVertex);
                }
            }
            candidateNeighbors = move(tmp);
        }
        return clique;
    };
    vector<Vt> vertices = graph.getVertex_Id();
    for (auto vertex : vertices)
    {
        //vertex = v.first;
        sum = 0;
        min_flow = 0;
        size_t score = 0;
        /*
         * Coger cada nodo como nodo de arranque
         */
        if (already_checked.find(vertex) == already_checked.end())
        {
            score += store_map.at(vertex);
            tmpClique = findMaxCliqueWithVertex(vertex, maxClique.size(), graph, store_map, &score);
            if ((tmpClique.size() >= CLICK_LIMIT) && (idCliques.find(sum) == idCliques.end())) {
                idCliques.emplace(sum);
                for (auto c: tmpClique)
                    store_map[c]-=min_flow;
                endCliques.push(pair < size_t, vector < Vt >> (score, tmpClique));
            }
        }
    }
    return endCliques;
}

class GatbGraph
{
public:
    typedef Node graphBU;

    class PairedInfoMarker
    {
    public:
        void setGraph(Graph & graph)
        {
            _markMap.clear();
            _hitMap.clear();
            // We insert all the nodes into our map.
            auto iter = [&] (Node& item)  {
                _markMap[item] = set<Node>();
                _hitMap[item] = map<Node, uint16_t>();
                GraphVector<Node> parents = graph.predecessors(item);
                for (size_t i = 0; i < parents.size(); ++i)
                    _parent_cliques[item][parents[i]] = vector<set<Node>>();
                _representativeMap[item] = item;
            };
            graph.iterator().iterate (iter);
        }

        void mark (const Node& key, Node & repr)
        {
            if (_hitMap[repr].find(key) != _hitMap[repr].end()) {
                _hitMap[repr][key]++;
            }else{
                _markMap[key].emplace(repr);
                _hitMap[repr][key] = 1;
            }
        }

        uint16_t getHitReprVal(Node & repr, Node & val)
        {
            return _hitMap[repr][val];
        }

        void setRepresentant(const Node & item, Node & representant)
        {
            _representativeMap[item] = representant;
        }

        Node getRepresentant(const Node & key)
        {
            return _representativeMap[key];
        }

        bool inRepresentants(const Node & key)
        {
            return (_representativeMap.find(key) != _representativeMap.end());
        }

        set<Node> getSet (const Node& item) const  {  return _markMap.find(item)->second;  }

        void addClique(const Node & node, const Node & parent, set<Node> clique)
        {
            _parent_cliques[node][parent].push_back(clique);
        }

        vector<set<Node>> getCliquesParent(const Node & node, const Node & parent)
        {
            return _parent_cliques[node][parent];
        }

        void removeCliquesParent( Node  node, Node parent, size_t first, size_t length)
        {
            (_parent_cliques[node][parent]).erase((_parent_cliques[node][parent]).begin() + first - length);
        }

    private:
        map<Node,set<Node>> _markMap;
        map<Node,map<Node, uint16_t>> _hitMap;
        map<Node, map<Node, vector<set<Node>>>> _parent_cliques;
        /*
         * Este mapa es temporal - mucha memoria innecesaria, fiel implementacion original
         */
        map<Node, Node> _representativeMap;
    };

    class FirstLastMarker
    {
    public:
        void mark (Node& key, Node & last, Graph g)  {
            _markMap [key] = last;
        }

        void buildTransMap(Graph g)
        {
            _currNodes = 0;
            for (auto kv:_markMap)
            {
                if (_transMap.find(kv.first) == _transMap.end() || _transMap.find(kv.second) == _transMap.end()) {
                    _transMap[kv.first] = _currNodes++;
                    if (kv.first != kv.second)
                        _transMap[kv.second] = _currNodes++;
                }
                cout << "CurrNodes: "<<_currNodes<<" "<<_transMap.size()<<endl;
            }
        }

        pair<Node, Node> firstLast (Node & key) const  {
            return {key,_markMap.find(key)->second};
        }

        size_t translateRepr(const Node & node)
        {
            map<Node, size_t>::iterator it_map = _transMap.find(node);
            return (it_map != _transMap.end())?(*it_map).second:INF;
        }

        size_t getNumRepr()
        {
            return _transMap.size();
        }

        void fillReach(vector<bool> & reach, Graph g)
        {
            vector<bool> checked(_currNodes, false);
            for (auto n: _transMap)
            {
                size_t distance = 0;
                checked[n.second] = true;
                _extent(checked, distance, n.first, n.second, n.first, g);
            }
        }

    private:
        size_t _in(Node & node)
        {
            map<Node, size_t>::iterator nodeIt = _transMap.find(node);
            if (nodeIt != _transMap.end())
                return (*nodeIt).second;
            return INF;
        }
        void _extent(vector<bool> & checked, size_t & distance, Node src, size_t src_trad, Node cur_node, Graph g)
        {
            if (distance >= (2*DELTA_PATH))
                return;
            distance++;
            GraphVector<Node> neighs = g.successors(cur_node);
            for (size_t i = 0; i < neighs.size(); ++i)
            {
                size_t trad = _in(neighs[i]);
                if (trad != INF)
                {
                    checked[(src_trad * _currNodes) + trad] = 1;
                    size_t newDistance = 0;
                    checked[trad] = true;
                    _extent(checked, newDistance, neighs[i], trad, neighs[i], g);
                }
                _extent(checked, distance, src, src_trad, neighs[i], g);
            }
        }

        map<Node,Node>  _markMap;
        map<Node, size_t> _transMap;
        size_t _currNodes;
    };

    GatbGraph(char * all_reads, char * read_left, char * read_right, bool refinement = false):_read_left(read_left),_read_right(read_right)
    {
        /*
         * Test
         */
        /*const char* sequences[] =
                {
                        "AGGCGCTAGGGAGAGGATGATGAAA",
                        "AGGCGCTCGGGTGAGGATGATGAAA",
                        "AGGCGCTTGGGTGAGGATGATGAAA"
                };*/
        cout << "Working on GatbGraph"<<endl;
        IBank * inputBank = Bank::open (all_reads);
        _g = Graph::create (inputBank,  "-kmer-size %d  -abundance-min 132  -verbose 0", KMERS);
        //_g = Graph::create(new BankStrings (sequences, ARRAY_SIZE(sequences)),  "-kmer-size %d  -abundance-min 1  -verbose 0", 7);
        if (refinement)
        {
            cout << "Simplifying the graph (same refinement as SPAdes)"<<endl;
            _g.simplify();
        }
        /*
         * Showing some information
         */
        GraphIterator<Node> it_test = _g.iterator();
        cout << "Number of vertices in graph: "<<it_test.size()<<endl;

        cout << "Setting graph"<<endl;
        _pairedInfo.setGraph(_g);
        cout << "Representative kmers"<<endl;
        _getRepresentatives();
        cout << "Adding pair end information"<<endl;
        _addPairedInfo(KMERS);
        cout << "End!"<<endl;
    }

    /*
     * Graph basic methods
     */
    size_t in_degree(Graph graph, graphBU node)
    {
        return graph.indegree(node);
    }

    size_t out_degree(Graph graph, graphBU node)
    {
        return graph.outdegree(node);
    }

    vector<graphBU> outKmerNeighbors(graphBU node)
    {
        vector<graphBU> results;
        GraphVector<graphBU> out_neighs = _g.successors(node);
        for (size_t i = 0; i < out_neighs.size(); ++i)
            results.push_back(out_neighs[i]);
        return results;
    }

    vector<graphBU> getKmerNeighbors(graphBU node)
    {
        return outKmerNeighbors(node);
    }

    vector<graphBU> getinKmerNeighbors(graphBU node)
    {
        vector<graphBU> results;
        GraphVector<graphBU> in_neighs = _g.predecessors(node);
        for (size_t i = 0; i < in_neighs.size(); ++i)
            results.push_back(in_neighs[i]);
        return results;
    }

    bool is_solid(graphBU node)

    {
        return _g.contains(node);
    }

    size_t getNumVertex()
    {
        return _g.iterator().size();
    }

    set<graphBU> getExtraInfoNode(graphBU node)
    {
        return _pairedInfo.getSet(node);
    }

    set<graphBU> getPairedInfoNode(graphBU node)
    {
        return _pairedInfo.getSet(node);
    }

    graphBU getRepresentant(graphBU node)
    {
        return _pairedInfo.getRepresentant(node);
    }

    pair<graphBU, graphBU> getFirstLast(graphBU node)
    {
        return _firstlastInfo.firstLast(node);
    }

    uint16_t getRepresentantHitsNode(graphBU key, graphBU val)
    {
        return _pairedInfo.getHitReprVal(key, val);
    }

    size_t getRepresentantTranslation(graphBU key)
    {
        return _firstlastInfo.translateRepr(key);
    }

    string printNode(graphBU node)
    {
        return _g.toString(node);
    }

    void addHaplotype(graphBU son, graphBU parent, set<graphBU> clique)
    {
        _pairedInfo.addClique(son, parent, clique);
    }

    vector<set<graphBU>> getCliquesWithParent(graphBU son, graphBU parent)
    {
        return _pairedInfo.getCliquesParent(son, parent);
    }

    void removeCliquesParent( Node  node, Node parent, size_t first, size_t length)
    {
        _pairedInfo.removeCliquesParent(node, parent, first, length);
    }

private:

    void _modify_info()
    {
        size_t num_vertex = getNumVertex(), num_representants = _firstlastInfo.getNumRepr();
        cout << "Num vertex: "<<num_vertex<<" Num representants: "<<_firstlastInfo.getNumRepr()<<endl;
        //Progression
        cout << "STAGE: Processing cliques\n";
        //Incremental information:
        size_t cont = 0;

        GraphIterator<Node> it = _g.iterator();
        for (it.first(); !it.isDone(); it.next())
        {
            Node v = it.item();
            auto node_extra = getExtraInfoNode(v);
            if (node_extra.empty())
                continue;
            auto endpoints = getKmerNeighbors(v);
            for (auto endpoint: endpoints)
            {
                size_t curr_node = 0;
                auto neigh_extra = getExtraInfoNode(endpoint);
                if (neigh_extra.empty())
                    continue;
                set <graphBU> local_vect = getUnion(node_extra,neigh_extra),last_set, removed_representant;
                //Con la nueva version esto no vale de nada
                map<graphBU, vector<graphBU>> store_map;
                map<graphBU, graphBU> translate_map;
                map<graphBU, size_t> representant_hits;
                set<graphBU> local_set;
                size_t max_size = 0;

                set<graphBU> to_rep;
                set<graphBU> endpoints;
                map<graphBU, size_t> representant_siege;
                for (auto s:local_vect) {
                    graphBU representative = getRepresentant(s);
                    if (store_map.find(representative) == store_map.end()) {
                        to_rep.emplace(representative);
                        (in(node_extra, s)?
                                representant_siege[representative] = 0:representant_siege[representative] = 1);
                        removed_representant.emplace(representative);
                        store_map[representative] = vector<graphBU>(1, s);
                        representant_hits[representative] = getRepresentantHitsNode(representative, v) + getRepresentantHitsNode(representative, endpoint);
                        translate_map[representative] = representative;
                        local_set.emplace(representative);
                    }else {
                        store_map[representative].push_back(s);
                        if (removed_representant.find(representative) != removed_representant.end())
                        {
                            if ((in(node_extra, s) && representant_siege[representative] == 1) ||
                                (in(neigh_extra, s) && representant_siege[representative] == 0))
                                removed_representant.erase(representative);
                        }
                    }
                }
                for (auto k:getKmerNeighbors(endpoint)) {
                    for (auto s: getExtraInfoNode(k)){
                        graphBU representative = getRepresentant(s);
                        if (representant_hits.find(representative) != representant_hits.end())
                            representant_hits[representative] += getRepresentantHitsNode(representative, k);
                        else{
                            to_rep.emplace(representative);
                            store_map[representative] = vector<graphBU>(1,s);
                            representant_hits[representative] = getRepresentantHitsNode(representative, k);
                            translate_map[representative] = representative;
                            local_set.emplace(representative);
                        }
                    }
                }
                for (auto k:getinKmerNeighbors(v)) {
                    for (auto s: getExtraInfoNode(k)){;
                        graphBU representative = getRepresentant(s);
                        if (representant_hits.find(representative) != representant_hits.end())
                            representant_hits[representative] += getRepresentantHitsNode(representative, k);
                        else{
                            store_map[representative] = vector<graphBU>(1,s);
                            translate_map[representative] = representative;
                            representant_hits[representative] = getRepresentantHitsNode(representative, k);
                            local_set.emplace(representative);
                        }
                    }
                }
                removed_representant.clear();
                for (auto p:store_map)
                    (p.second.size() > max_size)?max_size=p.second.size():max_size=max_size;
                for (auto p:store_map) {
                    if (p.second.size() < floor(max_size*MIN_SIZE_REP))
                        removed_representant.emplace(p.first);
                }
                local_vect = to_rep;
                for (auto s: removed_representant) {
                    local_vect.erase(s);
                    local_set.erase(s);
                }
                last_set = local_vect;
                typedef size_t vertex_graph;
                SimpleUndirectedGraph<graphBU> local_graph;
                //Cambiar por vector!
                map<size_t, vertex_graph> local_node_map;
                map<vertex_graph, size_t> representative_size;
                for (auto s:local_set)
                {
                    vertex_graph source = local_graph.addVertex(s);
                    representative_size[source] = representant_hits[s];
                    local_node_map[getRepresentantTranslation(s)] = source;
                }
                for (auto s:last_set) {
                    vertex_graph source = local_node_map[getRepresentantTranslation(s)];
                    for (auto t:local_vect) {
                        if (s == t)
                            continue;
                        bool reached = false;
                        pair<graphBU, graphBU> s_pair, t_pair;

                        s_pair = getFirstLast(s), t_pair = getFirstLast(t);
                        reached = (_reach[getRepresentantTranslation(s_pair.second)*num_representants+getRepresentantTranslation(t_pair.first)]
                                   || _reach[getRepresentantTranslation(t_pair.second)*num_representants+getRepresentantTranslation(s_pair.first)]);
                        //if (node_info.id < 20)
                        {
                            cout << "S_pair (s): " << printNode(s_pair.first) << " " << printNode(s_pair.second) << " Translation: "
                                 << getRepresentantTranslation(s_pair.second) << " "
                                 << getRepresentantTranslation(s_pair.first)
                                 << " Tpair (t):" << printNode(t_pair.second) << " " << printNode(t_pair.first) << " Translation: "
                                 << getRepresentantTranslation(t_pair.first)
                                 << " " << getRepresentantTranslation(t_pair.second) << " " << reached << endl;
                        }
                        if (reached) {
                            vertex_graph target = local_node_map[getRepresentantTranslation(t)];
                            local_graph.addEdge(source, target);
                        }
                    }
                    local_vect.erase(s);
                }
                map <vertex_graph, vector<Node>> clique_representation;
                set <vertex_graph> nodes_erase;
                set<size_t> idCliques;
                priority_queue < pair < size_t, vector < vertex_graph >> > output;
                size_t num_vertex_local = local_graph.getNumVertices(),
                        num_edges_local = local_graph.getNumEdges();
                if (((num_vertex_local * (num_vertex_local - 1)) / 2) == (num_edges_local)) {
                    if ((node_extra.size()*PARENT_SON_LIMIT) > neigh_extra.size())
                        continue;
                    set<graphBU> uniqueHaplotype;
                    for (auto s:last_set) {
                        for (auto p:store_map[s])
                            uniqueHaplotype.emplace(p);
                    }
                    #pragma omp critical(update)
                    {
                        addHaplotype(endpoint, v, uniqueHaplotype);
                    }
                } else {
                    output = findMaxClique < SimpleUndirectedGraph<graphBU>, vertex_graph> (local_graph, idCliques, representative_size);
                    set <vertex_graph> edge_transversed;
                    float total_size = 0.0, average_size = 0.0;
                    size_t num_clicks = 0;
                    if (output.empty()) {
                        vector<graphBU> vertices = local_graph.getVertex();
                        for (auto v2: vertices){
                            set<graphBU> local_haplotype;
                            for (auto n: store_map[translate_map[v2]]) {
                                local_haplotype.emplace(n);
                            }
                            #pragma omp critical(update)
                            {
                                addHaplotype(endpoint, v, local_haplotype);
                            }
                            total_size += local_haplotype.size();
                            num_clicks++;
                        }
                    } else {
                        while (!output.empty() & (edge_transversed.size() != local_graph.getNumVertices())) {
                            pair <size_t, vector<vertex_graph>> top_click = output.top();
                            output.pop();
                            vector <vertex_graph> clique = top_click.second;
                            set<graphBU> local_haplotype;
                            for (size_t i = 0; i < clique.size(); ++i) {
                                Node node_clique = local_graph.getNode(clique[i]);
                                for (auto n: store_map[translate_map[node_clique]]){
                                    local_haplotype.emplace(n);
                                }
                                edge_transversed.emplace(clique[i]);
                            }
                            #pragma omp critical(update)
                            {
                                addHaplotype(endpoint, v, local_haplotype);
                            }
                            total_size += local_haplotype.size();
                            num_clicks++;
                        }
                        if (output.empty() & (edge_transversed.size() != local_graph.getNumVertices())) {
                            vector<size_t> vertices = local_graph.getVertex_Id();
                            for (auto v2: vertices){
                                if (edge_transversed.find(v2) == edge_transversed.end()) {
                                    set<graphBU> local_haplotype;
                                    graphBU node_clique = local_graph.getNode(v2);
                                    for (auto n: store_map[translate_map[node_clique]]) {
                                        local_haplotype.emplace(n);
                                    }
                                    #pragma omp critical(update)
                                    {
                                        addHaplotype(endpoint, v, local_haplotype);
                                    }
                                    total_size += local_haplotype.size();
                                    num_clicks++;
                                }
                            }
                        }
                    }
                    average_size = total_size / num_clicks;
                    vector<set<graphBU>> cliquesAvail = getCliquesWithParent(endpoint,v);
                    set<graphBU> pairedInfoNeigh = getPairedInfoNode(endpoint),
                            pairedInfoNode = getPairedInfoNode(v);
                    size_t removes = 0, numberCliques = cliquesAvail.size();
                    vector<bool> fake_cliques(numberCliques, false);
                    vector<size_t> caso(numberCliques,6);
                    vector<size_t> join_cliques;
                    vector<set<Node>> intersections;
                    for (auto clique: cliquesAvail) {
                        intersections.push_back(getIntersection(clique, pairedInfoNeigh));
                    }
                    for (uint16_t i = 0; i < numberCliques; ++i)
                    {
                        if (fake_cliques[i])
                            continue;
                        size_t num_neighs = getKmerNeighbors(v).size(),
                                clique_size = cliquesAvail[i].size();
                        set<Node> parent_none = getIntersection(cliquesAvail[i],pairedInfoNode), neigh_none = intersections[i];
                        if (neigh_none.empty() || parent_none.empty()) {
                            fake_cliques[i] = true;
                            if (parent_none.empty())
                                caso[i] = 7;
                            else {
                                caso[i] = 1;
                                total_size -= clique_size;
                            }
                            removes++;
                        }else {
                            if (clique_size < (average_size * CLICK_RATIO)) {
                                fake_cliques[i] = true;
                                caso[i] = 4;
                            } else {
                                bool neigh_full = isSubset(cliquesAvail[i],pairedInfoNeigh, SIMILARITY_RATIO),
                                        parent_full = isSubset(cliquesAvail[i],pairedInfoNode, SIMILARITY_RATIO),
                                        parent_same = isSame(parent_none,pairedInfoNode,IDENTITY_RATIO);
                                if ((!neigh_full && parent_full) && !parent_same && (num_neighs > 1)) {
                                    fake_cliques[i] = true;
                                    caso[i] = 2;
                                    total_size -= clique_size;
                                    removes++;
                                } else {
                                    for (uint16_t j = 0;
                                         j < numberCliques; ++j) {
                                        if (i == j || fake_cliques[j])
                                            continue;
                                        if (isSubset(intersections[i], intersections[j]) &&
                                            (intersections[i].size() != intersections[j].size())) {
                                            fake_cliques[i] = true;
                                            caso[i] = 3;
                                            total_size -= clique_size;
                                            removes++;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    average_size = total_size / (num_clicks-removes);
                    for (size_t i = 0; i < fake_cliques.size(); i++)
                    {
                        if (fake_cliques[i] && caso[i] == 7){
                            if (_g.indegree(endpoint) == 1)
                                caso[i] = 4;
                        }
                        if (fake_cliques[i] && caso[i] == 4){
                            if ( cliquesAvail[i].size() > (average_size*CLICK_RATIO)) {
                                fake_cliques[i] = false;
                                caso[i] = 6;
                            }
                        }
                    }
                    size_t count = 0, count2 = 0;
                    for (auto i: fake_cliques) {
                        if (i) {
                            #pragma omp critical(update)
                            {
                                _pairedInfo.removeCliquesParent(endpoint, v, count, count2++);
                            }
                        }
                        count++;
                    }
                    fake_cliques.clear();
                }
            }
        }
        _reach.clear();
    }

    void _addPairedInfo(size_t kmer_size)
    {
        /*
         * Sequence iterator
         */
        IBank *bank1 = Bank::open(_read_left);
        IBank *bank2 = Bank::open(_read_right);
        PairedIterator <Sequence> *itPair = new PairedIterator<Sequence>(bank1->iterator(), bank2->iterator());
        ProgressIterator <std::pair<Sequence, Sequence>> progress_iter(itPair, "paired-end", bank1->estimateNbItems());
        /*
         * KmerSize model -> for kmers iteration
         */
        Kmer<128>::ModelDirect kmerModel (kmer_size);
        Kmer<128>::ModelDirect::Iterator kmerItLeft (kmerModel), kmerItRight (kmerModel);
        /*
         * Adds counter
         */
        size_t addsCounter = 0;
        for (progress_iter.first(); !progress_iter.isDone(); progress_iter.next())
        {
            Sequence& s1 = itPair->item().first;
            Sequence& s2 = itPair->item().second;
            kmerItLeft.setData(s1.getData());
            kmerItRight.setData(s2.getData());
            //We will replicate with kmerItRight
            kmerItRight.first();
            for (kmerItLeft.first(); !kmerItLeft.isDone();kmerItLeft.next())
            {
                string k1 = kmerModel.toString(kmerItLeft->value()), k2 = kmerModel.toString(kmerItRight->value());
                Node n_left = _g.buildNode(k1.c_str()), n_right = _g.buildNode(k2.c_str());
                if (_g.contains(n_left) && _g.contains(n_right))
                {
                    addsCounter++;
                    Node right_repr = _pairedInfo.getRepresentant(n_right);
                    _pairedInfo.mark(n_left, right_repr);
                }
                kmerItRight.next();
            }
        }
        cout << "Number of adds: "<<addsCounter<<endl;
    }

    void __extent(Node & n, Node & r, Node & p, set<Node> & visited, bool mark = true)
    {
        stack<Node> cur_stack, rep_stack, prev_stack;
        cur_stack.push(n);rep_stack.push(r);prev_stack.push(p);
        set<Node> current_path;
        while (!cur_stack.empty()) {
            //Stack managing
            Node node = cur_stack.top(), representant = rep_stack.top(), prev_node = prev_stack.top();
            cur_stack.pop();rep_stack.pop();prev_stack.pop();
            //Node stats
            size_t indegree = _g.indegree(node), outdegree = _g.outdegree(node);
            /*for (auto n: current_path)
                cout << " "<<_g.toString(n)<<" ";
            cout <<endl;
            cout << "Current node: "<<_g.toString(node)<<endl;
            cout << "Representant: "<<_g.toString(representant)<<endl;
            if (current_path.find(node) != current_path.end())
            {
                cout << "Ciclo found in: "<<_g.toString(node)<<endl;
                cout << "PrevNode: "<<_g.toString(prev_node)<<endl;
                cout << "Representant: "<<_g.toString(representant)<<endl;
                cout << "Degree: "<<indegree<<" "<<outdegree<<endl;
                GraphVector<Node> pred = _g.predecessors(node);
                for (size_t i = 0; i < pred.size(); ++i) {
                    cout << "Prev: "<<_g.toString(pred[i]) << endl;
                    cout << "Degree: "<<_g.indegree(pred[i])<<" "<<_g.outdegree(pred[i])<<endl;
                }
                pred = _g.successors(node);
                for (size_t i = 0; i < pred.size(); ++i) {
                    cout << "Next: "<<_g.toString(pred[i]) << endl;
                    cout << "Degree: "<<_g.indegree(pred[i])<<" "<<_g.outdegree(pred[i])<<endl;
                }
                exit(1);
            }*/
            current_path.emplace(node);
            if (indegree == 1 && outdegree <= 1) {
                GraphVector<Node> neighs = _g.successors(node);
                _pairedInfo.setRepresentant(node, representant);
                for (size_t i = 0; i < neighs.size(); ++i) {
                    cur_stack.push(neighs[i]);rep_stack.push(representant);prev_stack.push(node);
                };
            } else if (outdegree > 1 && indegree <= 1) {
                if (mark) {
                    _pairedInfo.setRepresentant(node, representant);
                    _firstlastInfo.mark(representant, node, _g);
                } else
                    mark = true;
                current_path.clear();
                if (visited.find(node) == visited.end()) {
                    visited.emplace(node);
                    GraphVector<Node> neighs = _g.successors(node);
                    for (size_t i = 0; i < neighs.size(); ++i){
                        cur_stack.push(neighs[i]);rep_stack.push(neighs[i]);prev_stack.push(neighs[i]);}
                }
            } else if (indegree > 1 && outdegree <= 1) {
                _firstlastInfo.mark(representant, prev_node, _g);
                _pairedInfo.setRepresentant(node, node);
                current_path.clear();
                if (visited.find(node) == visited.end()) {
                    visited.emplace(node);
                    GraphVector<Node> neighs = _g.successors(node);
                    for (size_t i = 0; i < neighs.size(); ++i) {
                        cur_stack.push(neighs[i]);rep_stack.push(node);prev_stack.push(node);
                    }
                }
            } else if (indegree > 1 && outdegree > 1) {
                _firstlastInfo.mark(representant, prev_node, _g);
                current_path.clear();
                if (visited.find(node) == visited.end()) {
                    visited.emplace(node);
                    _pairedInfo.setRepresentant(node, node);
                    _firstlastInfo.mark(node, node, _g);
                    GraphVector<Node> neighs = _g.successors(node);
                    for (size_t i = 0; i < neighs.size(); ++i){
                        cur_stack.push(neighs[i]);rep_stack.push(neighs[i]);prev_stack.push(neighs[i]);}
                }
            }
        }
    }

    void _getRepresentatives()
    {
        vector<Node> start_points;
        GraphIterator<Node> it = _g.iterator();
        set<Node> visited;
        for (it.first(); !it.isDone();it.next()) {
            Node current = it.item();
            if (_g.indegree(current) == 0) {
                start_points.push_back(current);
            } else if (_g.outdegree(current) > 1 && _g.indegree(current) == 1){
                start_points.push_back(current);
            }
        }
        cout << "Number of start points: "<<start_points.size()<<endl;
        size_t cont = 0;
        for (auto node : start_points) {
            cout << "Launch: "<<(cont++)<<endl;
            __extent(node, node, node, visited, false);
        }
        cout << "Pending cycles"<<endl;
        it = _g.iterator();
        for (it.first(); !it.isDone(); it.next())
        {
            if (!_pairedInfo.inRepresentants(it.item()))
            {
                cout << _g.toString(it.item())<<" Not transversed!"<<endl;
            }
        }
        _firstlastInfo.buildTransMap(_g);
        cout << "Number of representants: "<<_firstlastInfo.getNumRepr()<<endl;
    }

    void _fillReachMatrix()
    {
        set<Node> transversed;
        _firstlastInfo.fillReach(_reach, _g);
    }

    Graph _g;
    PairedInfoMarker _pairedInfo;
    FirstLastMarker _firstlastInfo;
    vector<bool> _reach;
    char * _read_left, * _read_right;
};
typedef Node graphBU;
size_t in_degree(Graph graph, graphBU node)
{
    return graph.indegree(node);
}

size_t out_degree(Graph graph, graphBU node)
{
    return graph.outdegree(node);
}

vector<graphBU> outKmerNeighbors(Graph graph, graphBU node)
{
    vector<graphBU> results;
    GraphVector<graphBU> out_neighs = graph.successors(node);
    for (size_t i = 0; i < out_neighs.size(); ++i)
        results.push_back(out_neighs[i]);
    return results;
}

vector<graphBU> inKmerNeighbors(Graph graph, graphBU node)
{
    vector<graphBU> results;
    GraphVector<graphBU> in_neighs = graph.predecessors(node);
    for (size_t i = 0; i < in_neighs.size(); ++i)
        results.push_back(in_neighs[i]);
    return results;
}

bool is_solid(Graph graph, graphBU node)
{
    return graph.contains(node);
}

void addPairedInfo(Graph graph, char * read_left, char * read_right, size_t kmer_size) {
    class GraphMarker
    {
    public:
        GraphMarker (const Graph& graph) : graph(graph)
        {
            // We insert all the nodes into our map.
            auto iter = [&] (const Node& item)  {
                this->markMap[item] = set<Node>();
                this->hitMap[item] = map<Node, uint16_t>();
            };
            graph.iterator().iterate (iter);
        }
        void mark (const Node& item, Node & val)
        {
            if (hitMap[item].find(val) != hitMap[item].end()) {
                hitMap[item][val]++;
            }else{
                markMap[item].emplace(val);
                hitMap[item][val] = 1;
            }
        }
        set<Node> getSet (const Node& item) const  {  return markMap.find(item)->second;  }
    private:
        const Graph& graph;
        map<Node,set<Node>> markMap;
        map<Node,map<Node, uint16_t>> hitMap;
    };
    /*
     * Sequence iterator
     */
    IBank *bank1 = Bank::open(read_left);
    IBank *bank2 = Bank::open(read_right);
    PairedIterator <Sequence> *itPair = new PairedIterator<Sequence>(bank1->iterator(), bank2->iterator());
    ProgressIterator <std::pair<Sequence, Sequence>> progress_iter(itPair, "paired-end", bank1->estimateNbItems());
    /*
     * KmerSize model -> for kmers iteration
     */
    Kmer<128>::ModelDirect kmerModel (kmer_size);
    Kmer<128>::ModelDirect::Iterator kmerItLeft (kmerModel), kmerItRight (kmerModel);
    /*
     * Marker
     */
    GraphMarker marker(graph);
    /*
     * Adds counter
     */
    size_t addsCounter = 0;
    for (progress_iter.first(); !progress_iter.isDone(); progress_iter.next())
    {
        Sequence& s1 = itPair->item().first;
        Sequence& s2 = itPair->item().second;
        kmerItLeft.setData(s1.getData());
        kmerItRight.setData(s2.getData());
        //We will replicate with kmerItRight
        kmerItRight.first();
        for (kmerItLeft.first(); !kmerItLeft.isDone();kmerItLeft.next())
        {
            string k1 = kmerModel.toString(kmerItLeft->value()), k2 = kmerModel.toString(kmerItRight->value());
            Node n_left = graph.buildNode(k1.c_str()), n_right = graph.buildNode(k2.c_str());
            if (graph.contains(n_left) && graph.contains(n_right))
            {
                addsCounter++;
                marker.mark(n_left, n_right);
            }
            kmerItRight.next();
        }
    }
    cout << "Number of adds: "<<addsCounter<<endl;
}

void getRepresentatives(Graph graph)
{
    class GraphMarker
    {
    public:
        GraphMarker (const Graph& graph) : graph(graph)
        {
        }
        void mark (Node& key, Node & last)  {  markMap [key] = last;  }
        pair<Node, Node> lastFirst (BranchingNode& key) const  {  return {key,markMap.find(key)->second};  }
    private:
        const Graph& graph;
        map<Node,Node>  markMap;
    };
    GraphIterator<BranchingNode> it = graph.iteratorBranching ();
    set<BranchingNode> representatives_seted;
    /*
     * Marker
     */
    GraphMarker marker(graph);
    for (it.first(); !it.isDone(); it.next())
    {
        /*
         * Lets do this
         */
        GraphIterator<Node> path = graph.simplePath(it.item(), (graph.outdegree(it.item()) > 1)?DIR_INCOMING:DIR_OUTCOMING);
        path.first();
        Node first_node = path.item(), last_node;
        path.next();
        if (path.isDone()) {
            graph.setRepresentant(it.item(), it.item());
            marker.mark(first_node, first_node);
            /*
             * Test
             */
            Node fake_node = graph.buildNode(graph.toString(it.item()).c_str());
            if (fake_node == it.item())
            {
                cout << "Representant: "<<graph.toString(it.item())<<" "<<graph.toString(graph.getRepresentant(it.item()))<<endl;
                cout << "Representant: "<<graph.toString(fake_node)<<" "<<graph.toString(graph.getRepresentant(fake_node))<<endl;
            }
        } else {
            vector<Node> local_unitig;
            set<Node> cycle;
            local_unitig.push_back(first_node);
            while (!path.isDone()) {
                graph.setRepresentant(path.item(), first_node);
                if (cycle.find((*path)) != cycle.end()){
                    cout << "Cycle found!" <<it.rank()<< endl;
                    break;
                }else
                    cycle.emplace((*path));
                local_unitig.push_back(path.item());
                path.next();
            }
            marker.mark(first_node, local_unitig[local_unitig.size()-2]);
        }
        /*
         * Check 1 - representative is only assigned once
         */
        if (representatives_seted.find((*it)) != representatives_seted.end())
            cout << "Representant assigned twice: "<<graph.toString(it.item())<<endl;
        else
            representatives_seted.emplace((*it));
    }
    cout << "Number of representants: "<<representatives_seted.size()<<endl;
}

void iterateGraph(Graph graph)
{
    GraphIterator<Node> it = graph.iterator ();
    for (it.first(); !it.isDone(); it.next())
    {
        Node & current = it.item();
        if (graph.outdegree(current) > 1 || graph.indegree(current) > 1) {
            std::cout << "Branching node: " << std::endl;
            std::cout << graph.toString(current) << std::endl;
        }
    }
}

void parse_args_gatb(int argc, char **argv, std::string * path_to_file, std::string * dir_pairs, std::string * path_to_write
        ,std::string * path_unitigs, std::string * program, bool * pair_end, bool * thirdPartyCount, bool * do_correction
        ,bool * do_polish, bool * meta)
{
    /*OptionsParser parser("GatbParser");
    parser.push_back(new OptionOneParam("-s","Path single end file",false));
    parser.push_back(new OptionOneParam("-p","Paired_end folder",false));
    parser.push_back(new OptionOneParam("-o","Output folder",false));
    parser.push_back(new OptionOneParam("-u","Unitigs path to write",false));
    parser.push_back(new OptionOneParam("-b","Do correction",false));
    parser.push_back(new OptionOneParam("-c","K-mer counting third party program",false));*/

    OptionsParser parser ("KmerTest");
    parser.push_back (new OptionOneParam (STR_URI_INPUT,      "bank input",     true));
    parser.push_back (new OptionOneParam (STR_KMER_SIZE,      "kmer size",      true));
    parser.push_back (new OptionOneParam (STR_MINIMIZER_SIZE, "minimizer size", true));
}
int main (int argc, char* argv[])
{
    size_t kmerSize = 120;
    /*std::string path_to_file(""), dir_pairs(""), output_path, path_unitigs, program="dsk";
    bool pair_end = false, thirdPartyCount = true, do_correction = false, do_polish = false, meta = false;
    parse_args_gatb(argc,argv, &path_to_file, &dir_pairs, &output_path,&path_unitigs,&program,
                      &pair_end, &thirdPartyCount,&do_correction,&do_polish,&meta);*/
    const char * non_sequences[] =
            {
                "AGCCCGATTATACCGA"
            };
    const char* sequences[] =
            {
                    //      x <- difference here
                    "AGGCGCTAGGGAGAGGATGATGAAA",
                    "AGGCGCTCGGGAGAGGATGATGAAA",
                    "AGGCGCTTGGGAGAGGATGATGAAA"
            };
    try
    {
        char * all_reads = argv[1], * read_left = argv[2], * read_right = argv[3];
        cout << "Start building: GraphGatb"<<endl;
        auto start_time = std::chrono::high_resolution_clock::now();
        GatbGraph g(all_reads, read_left, read_right);
        auto end_time = std::chrono::high_resolution_clock::now();
        cout << "End building"<<endl;
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        cout << "Total time: "<<ms<<endl;
        exit(1);
        // We create the graph.
        cout << "Start building"<<endl;
        start_time = std::chrono::high_resolution_clock::now();
        IBank * inputBank = Bank::open (all_reads);
        Graph graph = Graph::create (inputBank,  "-kmer-size %d  -abundance-min 132  -verbose 0", kmerSize);
        //Graph graph = Graph::create(new BankStrings (sequences, ARRAY_SIZE(sequences)),  "-kmer-size %d  -abundance-min 1  -verbose 0", kmerSize);
        end_time = std::chrono::high_resolution_clock::now();
        cout << "End building"<<endl;
        ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        cout << "Time build: "<<ms<<endl;
        /*
         * Based on SPAdes approach
         *  1. Tip removal
         *  2. Bulges removal
         *  3. EC removal
         */
        cout << "Start polishing"<<endl;
        start_time = std::chrono::high_resolution_clock::now();
        graph.simplify();
        end_time = std::chrono::high_resolution_clock::now();
        cout << "End polishing"<<endl;
        ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        cout << "Time Simplify: "<<ms<<endl;
        /*
        * Iterate graph
        */
        /*cout << "Start iterating"<<endl;
        start_time = std::chrono::high_resolution_clock::now();
        iterateGraph(graph);
        end_time = std::chrono::high_resolution_clock::now();
        cout << "End iterating"<<endl;
        ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        cout << "Time Iterate: "<<ms<<endl;*/
        /*
         * Check some stats
         */
        GraphIterator<Node> it_test = graph.iterator();
        cout << "Number of vertices in graph: "<<it_test.size()<<endl;
        cout << "Start representatives"<<endl;
        start_time = std::chrono::high_resolution_clock::now();
        getRepresentatives(graph);
        end_time = std::chrono::high_resolution_clock::now();
        cout << "End representatives"<<endl;
        ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        cout << "Time Representatives: "<<ms<<endl;
        /*
         * Adding paired_end
         */
        cout << "Start adding"<<endl;
        start_time = std::chrono::high_resolution_clock::now();
        addPairedInfo(graph, read_left, read_right, kmerSize);
        end_time = std::chrono::high_resolution_clock::now();
        cout << "End adding"<<endl;
        ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        cout << "Time Add: "<<ms<<endl;
        exit(1);
        /*
         * Testing abundances:
         */
        Node node = graph.buildNode (sequences[0]), fake_node = graph.buildNode (non_sequences[0]);
        if (graph.contains(node))
        {
            std::cout << "Node exists!"<<graph.toString(node) << std::endl;
            int abundance = graph.queryAbundance(node);
            std::cout << "Abundancia: "<<abundance<<std::endl;
            std::cout << "Indegree: "<<in_degree(graph, node)<<std::endl;
            std::cout << "Degree: "<<out_degree(graph, node)<<std::endl;
            graph.setRepresentant(node, node);
        }
        if (!graph.contains(fake_node))
            std::cout << "Node does not exist!"<<std::endl;
        // We retrieve the branching neighbors for the node.
        std::cout << "Starting branching nodes"<<std::endl;
        GraphVector<BranchingNode> branchingNeighbors = graph.successorsBranching (node);
        GraphVector<BranchingNode> branchingInNeighbors = graph.predecessorsBranching (node);
        std::cout << "We found " << branchingNeighbors.size() << " branching neighbors from node " << graph.toString(node) << std::endl;
        // We loop over the branching neighbors. Here, we should have 3 branching neighbors, being the same GGGAGAG
        for (size_t i=0; i<branchingNeighbors.size(); i++)
        {
            std::cout << graph.toString (branchingNeighbors[i])  << std::endl;
        }
        std::cout << "End branching nodes test"<<std::endl;
        //Getting representatives:
        std:cout<<"Starting representatives nodes test"<<std::endl;
        GraphIterator<Node> it = graph.iterator ();
        for (it.first(); !it.isDone(); it.next())
        {
            Node& current = it.item();
            Node & last = (graph.outdegree(current) != 1 || graph.indegree(current) != 1)
                    ?current:graph.successorsBranching(current)[0];
            Node & first = (graph.outdegree(current) != 1 || graph.indegree(current) != 1)
                           ?current:graph.predecessorsBranching(current)[0];
            std::cout << "Node: "<<graph.toString(current)
                <<" First: " <<graph.toString(first)<<" Last: "<<graph.toString(last)<<std::endl;
        }
        std::cout<<"End representatives nodes test"<<std::endl;
        //Trying unitig graphs
        std::cout << "Starting unitig graph trial"<<std::endl;

        std::cout << "End unitig graph trial"<<std::endl;
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }
    return EXIT_SUCCESS;
}