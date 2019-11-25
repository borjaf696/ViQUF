// Created by borja on 8/11/19.
//
#ifndef GATB_TRIAL_GRAPH_H
#include "graph.h"
#endif

size_t split(const std::string &txt, std::vector<std::string> &strs, char ch)
{
    size_t pos = txt.find( ch );
    size_t initialPos = 0;
    strs.clear();

    // Decompose statement
    while( pos != std::string::npos ) {
        strs.push_back( txt.substr( initialPos, pos - initialPos ) );
        initialPos = pos + 1;

        pos = txt.find( ch, initialPos );
    }

    // Add the last one
    //strs.push_back( txt.substr( initialPos, std::min( pos, txt.size() ) - initialPos + 1 ) );

    return strs.size();
}

DBG::DBG(const DBG & dbg):_g(dbg._g),_numedges(dbg._numedges),_pairedEndInformation(dbg._pairedEndInformation)
{}

DBG::DBG(char * file):_numedges(0)
{
    auto start = chrono::steady_clock::now();
    ifstream inFile;
    inFile.open(file);

    if (!inFile)
    {
        cout << "Unable to open file"<<endl;
        exit(1);
    }

    string line;
    size_t cont = 0;
    while (std::getline(inFile, line))
    {
        if (cont++ == 0)
        {
            istringstream iss(line);
            iss >> _num_nodes;
            cout << "Number of nodes: "<<_num_nodes<<endl;
            _g = vector<vector<size_t>>(_num_nodes, vector<size_t>());
            _pairedEndInformation = vector<unordered_set<size_t>>(_num_nodes, unordered_set<size_t>());
            _abundance_map_node = vector<unordered_map<OwnNode_t, size_t>>(_num_nodes, unordered_map<OwnNode_t, size_t>());
            _reachibility = vector<bool>(_num_nodes*_num_nodes, false);
        }else {
            vector <string> results;
            split(line, results, ' ');
            istringstream iss(results[0]);
            size_t node;
            iss >> node;
            if (results.size() > 1) {
                for (auto s: results) {
                    size_t target;
                    istringstream iss(s);
                    iss >> target;
                    if (target == node)
                        continue;
                    size_t neighbor;
                    iss = istringstream(s);
                    iss >> neighbor;
                    _g[node].push_back(neighbor);
                }
            }
        }
    }
    inFile.close();
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds (graph construction from file) : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
    /*
     * Getting full reachability matrix
     */
    _complete_reach_matrix();
}

size_t DBG::edges()
{
    return _numedges;
}


size_t DBG::vertices()
{
    return _g.size();
}

vector<size_t> DBG::getNeighbor(size_t v , bool dir)
{
    if (dir)
        return _g[v];
    return _g[v];
}

unordered_set<size_t> DBG::getPairedEndInformation(size_t v)
{
    return _pairedEndInformation[v];
}

void DBG::addPair(size_t v, size_t p)
{
    _pairedEndInformation[v].emplace(p);
}

size_t DBG::addPair(size_t v, unordered_set <size_t> pairInformation)
{
    _pairedEndInformation[v] = pairInformation;
    return pairInformation.size();
}

void DBG::build_process_cliques()
{
    for (size_t i = 0; i < _num_nodes; ++i)
    {
        if (_pairedEndInformation[i].empty())
            continue;
        for (auto n: getNeighbor(i))
        {
            if (_pairedEndInformation[n].empty())
                continue;
            Pairedendinformation_t node_pei = getPairedEndInformation(i), neigh_pei = getPairedEndInformation(n);
            Pairedendinformation_t union_pei;
            set_union(node_pei.begin(), node_pei.end(),
                    neigh_pei.begin(), neigh_pei.end(),
                    inserter(union_pei,union_pei.begin()));
            UG local_graph;
            for (auto pn: union_pei)
                local_graph.addVertex(pn);
            local_graph.print();
        }
    }
}
/*
 * Private
 */
void DBG::_complete_reach_matrix()
{
    std::function<void(size_t, size_t,size_t, size_t&)> __reach = [this, &__reach](size_t  src, size_t target,
            size_t steps, size_t &branches)->void{
        if (steps++ == MAX_PATH)
            return;
        if (branches == MAX_BRANCHES)
            return;
        //cout << src<<" - "<<target<<" - "<<steps<<" - "<<branches<<endl;
        this->_reachibility[src*_num_nodes+target] = true;
        auto neighbors = getNeighbor(target);
        branches += (neighbors.size() > 1)?1:0;
        for (auto neigh: neighbors)
            __reach(src, neigh, steps, branches);
    };
    auto start = chrono::steady_clock::now();
    cout << "Building reachiability matrix..."<<endl;
    for (size_t i = 0; i < _num_nodes; ++i)
    {
        size_t distance = 0, branches = 0;
        __reach(i,i,distance, branches);
    }
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds (Completing reach matrix) : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
}
/*
 * Show information
 */
string DBG::toString()
{
    return string("Not supported yet");
}

void DBG::print()
{
    size_t vertex = 0;
    for (auto v:_g)
    {
        cout << "Node: "<<vertex<<" - ";
        for (auto n:v)
            cout << " "<<n;
        cout << endl;
        cout << "Paired-end information: ";
        for (auto n: _pairedEndInformation[vertex])
            cout << " "<<n;
        cout <<endl;
        vertex++;
    }
}

/*
 * Undirected Graphs
 */
void UG::addVertex(OwnNode_t node)
{
    _g_nodes.push_back(Node(_num_vertex++, node));
    _g_edges.push_back(vector<OwnNode_t>());
}

void UG::addVertexWithAbundances(OwnNode_t node, size_t ab)
{
    _g_nodes.push_back(Node(_num_vertex++,node,ab));
    _g_edges.push_back(vector<OwnNode_t>());
}

bool UG::setAbundance(OwnNode_t node, size_t ab)
{
    for (size_t i = 0; i < _g_nodes.size(); ++i)
    {
        if (_g_nodes[i]._val == node)
        {
            _g_nodes[i]._abundance = ab;
            return true;
        }
    }
    return false;
}

void UG::addEdge(OwnNode_t i, OwnNode_t j)
{
    _num_edges++;
    _g_edges[i].push_back(j);
    _g_edges[j].push_back(i);
}

vector<OwnNode_t > UG::getNeighbors(OwnNode_t i)
{
    return _g_edges[i];
}

vector<UG::Node> UG::getNodeNeighbors(OwnNode_t node)
{
    vector<Node> results;
    for (auto i: getNeighbors(node))
        results.push_back(_g_nodes[i]);
    return results;
}

priority_queue<pair<size_t,vector<OwnNode_t>>> UG::report_maximal_cliques()
{
    if (!_num_edges)
        return std::priority_queue<pair<size_t,vector<OwnNode_t>>>();
    priority_queue<pair<size_t,vector<OwnNode_t>>> endCliques;
    vector<OwnNode_t> maxClique, tmpClique;
    unordered_set<OwnNode_t> already_checked;
    set<size_t> idCliques;
    /*
     * Original version!
     */
    int64_t sum;
    auto __get_maximal_clique = [this, &already_checked, &sum](OwnNode_t vertex,const int maxCliqueSize, size_t * score) {
        vector<OwnNode_t> clique;
        clique.emplace_back(vertex);
        sum |= 1 << _g_nodes[vertex]._id;

        set<OwnNode_t> candidateNeighbors;

        unordered_set<OwnNode_t> visited;
        visited.emplace(vertex);

        for (auto n:getNeighbors(vertex))
            candidateNeighbors.emplace(n);

        set<OwnNode_t> tmp;

        while (!candidateNeighbors.empty()) {
            const auto highestDegNeighborIt = std::max_element(candidateNeighbors.begin(), candidateNeighbors.end(), [this](const OwnNode_t &lhs, const OwnNode_t &rhs) {
                if ((degree(lhs)) == (degree(rhs)))
                    return _g_nodes[lhs]._id > _g_nodes[rhs]._id;
                return (degree(lhs)) < (degree(rhs));
            });

            const auto highestDegVert = *highestDegNeighborIt;
            clique.emplace_back(highestDegVert);
            (*score) += degree(highestDegVert);
            sum |= 1 << highestDegVert;
            visited.emplace(highestDegVert);

            for (auto n: getNeighbors(highestDegVert))
            {
                if (candidateNeighbors.find(n)!=candidateNeighbors.end() && visited.find(n) == visited.end()) {
                    tmp.insert(n);
                }
            }
            candidateNeighbors = std::move(tmp);
        }
        return clique;

    };
    for (size_t vertex = 0; vertex < _g_edges.size(); vertex++)
    {
        sum = 0;
        size_t score = 0;
        if (already_checked.find(vertex) == already_checked.end())
        {
            tmpClique = __get_maximal_clique(vertex, maxClique.size(), &score);
            if ((tmpClique.size() >= CLIQUE_LIMIT) && (idCliques.find(sum) == idCliques.end())) {
                idCliques.emplace(sum);
                endCliques.push(pair < size_t, vector < OwnNode_t >> (score, tmpClique));
            }
        }
    }
    return endCliques;
}