// Created by borja on 8/11/19.
//
#include "graph.h"
/*
 * Esto es una castaña
 */
void write_histogram(vector<size_t> histogram, string outputfile = "histogram.txt")
{
    ofstream outputFile(outputfile.c_str());
    for (size_t i = 0; i < histogram.size(); ++i)
        outputFile << i<<" "<<histogram[i]<<"\n";
    outputFile.close();
}

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

DBG::DBG(const DBG & dbg):_g_nodes(dbg._g_nodes),_g_edges(dbg._g_edges),_numedges(dbg._numedges),_reachability(dbg._reachability)
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
            _g_edges = vector<vector<OwnNode_t>>(_num_nodes, vector<size_t>());
            _g_in_edges = vector<vector<OwnNode_t>>(_num_nodes, vector<size_t>());
            _reachability = vector<unordered_set<size_t>>(_num_nodes, unordered_set<size_t>());
            _num_nodes = 0;
        }else {
            vector <string> results;
            split(line, results, ' ');
            istringstream iss(results[0]);
            size_t node;
            iss >> node;
            _g_nodes.push_back(UG_Node(_num_nodes, _num_nodes));
            _num_nodes++;
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
                    _g_edges[node].push_back(neighbor);
                    _g_in_edges[neighbor].push_back(node);
                }
            }
        }
    }
    cout << "Sanity (number of nodes): "<<_num_nodes<<endl;
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
    return _g_nodes.size();
}

vector<size_t> DBG::getNeighbor(size_t v , bool dir)
{
    if (v >= _g_nodes.size())
        throw std::invalid_argument( "Out of bounds");
    if (dir)
        return _g_edges[v];
    return _g_in_edges[v];
}

unordered_set<size_t> DBG::getPairedEndInformation(size_t v)
{
    if (v >= _g_nodes.size())
        throw std::invalid_argument( "Out of bounds");
    return _g_nodes[v].get_paired_information();
}

size_t DBG::getAbundance(size_t node, size_t pair)
{
    if (node >= _g_nodes.size())
        throw std::invalid_argument( "Out of bounds");
    return _g_nodes[node].get_abundance(pair);
}

void DBG::addPair(size_t v, size_t p)
{
    _g_nodes[v].add_paired_information(p);
}

size_t DBG::addPair(size_t v, Pairedendinformation_t pairInformation)
{
    _g_nodes[v].add_paired_information(pairInformation);
    return pairInformation.size();
}

void DBG::build_process_cliques(DBG & apdbg)
{
    cout << "Number of nodes (in the AdBG): "<< _num_nodes - 1<<endl;
    bool show = false;
    size_t empty_pairs = 0;
    OwnNode_t node_test = 13341;
    string progress_bar;
    for (size_t i = 0; i < _num_nodes-1; ++i)
    {
        vector<size_t> local_histogram = vector<size_t>(10000);
        show = (_g_nodes[i]._id == node_test);
        //show = (_g_nodes[i]._val == node_test);
        Pairedendinformation_t node_pei = getPairedEndInformation(i);
        UG_Node parent(_g_nodes[i]._val, node_pei);
        //cout << " Node: "<<i<<" Show: "<<show<<" Pe: "<<node_pei.size()<<endl;
        /*if (show) {
            cout << "Size Node (paired-end): " << node_pei.size() << " Neighbors: " << getNeighbor(i).size() << endl;
        }*/
        vector<OwnNode_t> neighbors = getNeighbor(i);
        if (neighbors.size() == 0) {
            apdbg.addNode(parent, show);
        }
        /*if (node_pei.empty()) {
            for (auto n: neighbors)
            {
                UG_Node son(_g_nodes[n]._val, getPairedEndInformation(n));
                apdbg.addNode(parent, son);
            }
            continue;
        }*/
        if (node_pei.empty())
            empty_pairs++;
        for (auto n: neighbors)
        {
            //show = (_g_nodes[n]._val == node_test);
            Pairedendinformation_t neigh_pei = getPairedEndInformation(n);
            /*if (show)
                cout << "Size (paired-end): "<<neigh_pei.size()<<endl;*/
            /*if (_g_nodes[n].get_paired_information().empty()) {
                UG_Node son(_g_nodes[n]._val, getPairedEndInformation(n));
                apdbg.addNode(parent, son);
                continue;
            }*/
            Pairedendinformation_t union_pei = union_set(node_pei, neigh_pei);
            if (COMPLETE)
            {
                vector<OwnNode_t> parents = getNeighbor(i, false), grand_son = getNeighbor(n);
                if (parents.size() != 0) {
                    for (auto p:parents)
                        union_pei = union_set(union_pei, getPairedEndInformation(p));
                }
                if (grand_son.size() != 0){
                    for (auto s:grand_son)
                        union_pei = union_set(union_pei, getPairedEndInformation(s));
                }
            }
            UG local_graph;
            unordered_map<size_t, size_t> traslation_map, reverse_map;
            for (auto pn: union_pei){
                local_histogram[std::min((size_t) 999,getAbundance(n,pn) + getAbundance(i,pn))]++;
                size_t id = local_graph.addVertex(pn);
                traslation_map[pn] = id;
                reverse_map[id] = pn;
            }
            for (auto pn: union_pei)
            {
                for (auto tn: union_pei)
                {
                    if (_reach(pn,tn)) {
                        local_graph.addEdge(traslation_map[pn], traslation_map[tn]);
                    }
                }
            }
            if (show) {
                cout << "Graph: "<<endl;
                local_graph.print();
            }
            priority_queue<pair<size_t,vector<OwnNode_t>>> cliques = local_graph.report_maximal_cliques(reverse_map, *this);
            unordered_set<OwnNode_t> nodes_checked;
            if (show)
                cout << " Cliques size: "<<cliques.size()<<endl;
            vector<pair<UG_Node,UG_Node>> potential_nodes;
            while (!cliques.empty() && nodes_checked.size() < union_pei.size())
            {
                pair <size_t, vector<OwnNode_t>> top_click = cliques.top();
                cliques.pop();
                Pairedendinformation_t clique, p_info_a, p_info_b;
                bool in_a = true, in_b = true;
                //cout << endl;
                for (auto node: top_click.second)
                {
                    nodes_checked.emplace(node);
                    UG::UG_Node real_node = local_graph[node];
                    clique.emplace(real_node._val);
                    //cout << node << " ";
                    if (node_pei.find(real_node._val) != node_pei.end()) {
                        in_a = true;
                        p_info_a.emplace(real_node._val);
                    }
                    if (neigh_pei.find(real_node._val) != neigh_pei.end()) {
                        in_b = true;
                        p_info_b.emplace(real_node._val);
                    }
                }
                //cout << endl;
                bool check_both = true;
                if (!node_pei.empty() and !neigh_pei.empty())
                    check_both = false;
                if (check_both or (in_a && in_b))
                {
                    UG_Node parent_l(_g_nodes[i]._val, p_info_a), son_l(_g_nodes[n]._val,p_info_b);
                    potential_nodes.push_back(pair<UG_Node,UG_Node>(parent_l,son_l));
                    //apdbg.addNode(parent_l, son_l);
                }
            }
            /*
             * Procesar las paired-end y eliminar los subconjuntos
             */
            vector<bool> offset_erase = vector<bool>(potential_nodes.size(), false);
            for (size_t i = 0; i < potential_nodes.size(); ++i)
            {
                for (size_t j = i + 1; j < potential_nodes.size(); ++j)
                {
                    Pairedendinformation_t p1 = potential_nodes[i].first._paired_info, p2 = potential_nodes[j].first._paired_info
                            , s1 = potential_nodes[i].second._paired_info, s2 = potential_nodes[j].second._paired_info;
                    if (p1 != p2)
                    {
                        if (is_subset(p1,p2))
                            offset_erase[i] = true;
                        if (is_subset(p2,p1))
                            offset_erase[j] = true;
                    }
                    if (s1 != s2)
                    {
                        if (is_subset(s1,s2))
                            offset_erase[i] = true;
                        if (is_subset(s2,s1))
                            offset_erase[j] = true;
                    }
                }
            }
            size_t num_removes = 0;
            for (size_t pos = 0; pos < offset_erase.size(); pos++) {
                if (offset_erase[pos]) {
                    potential_nodes.erase(potential_nodes.begin() + pos - num_removes++);
                }
            }
            for (auto s: potential_nodes)
                apdbg.addNode(s.first,s.second);
        }
        /*
         * Crear directorios antes!
         */
        //string histogram_output = "stats/histograms_nodes/histogram_node"+to_string(i)+".txt";
        //write_histogram(local_histogram, histogram_output);
    }
    /*cout << "Displaying ADPBG"<<endl;
    apdbg.print();*/
    cout << "Number of empty pairs: "<<empty_pairs<<" out of "<<_g_nodes.size()<<endl;
}

size_t DBG::out_degree(OwnNode_t node_id)
{
    return _g_edges[node_id].size();
}

size_t DBG::in_degree(OwnNode_t node_id)
{
    return _g_in_edges[node_id].size();
}

vector<vector<OwnNode_t>> DBG::export_unitigs()
{
    vector<vector<OwnNode_t>> unitigs;
    vector<vector<OwnNode_t>> vector_unitigs(_g_nodes.size());
    vector<UG_Node> starting_points = _get_starting_nodes();
    queue<UG_Node> extensible_points;
    vector<bool> checked(_g_nodes.size(), false);
    size_t higher_out = 0, higher_in = 0, zero = 0, zero_nopair = 0;
    for (auto i: starting_points) {
        extensible_points.push(i);
    }
    cout << "Retrieving unitigs"<<endl;
    while(!extensible_points.empty())
    {
        UG_Node node = extensible_points.front();
        extensible_points.pop();
        cout << "Lanzó!"<<endl;
        if (checked[node._id])
            continue;
        vector<OwnNode_t> current_unitig;
        vector<bool> local_checked(_g_nodes.size(), false);
        size_t repetitive_call = 0;
        _extension(extensible_points, checked, current_unitig, node,
                higher_out, zero, higher_in, zero_nopair, repetitive_call, unitigs, local_checked);
        //unitigs.push_back(current_unitig);
    }
    cout << "Higher out degree: "<<higher_out<<endl;
    cout << "Higher in degree: "<<higher_in<<endl;
    cout << "Zero case: "<<zero<<endl;
    cout << "Zero no pair case: "<<zero_nopair<<endl;
    return unitigs;
}

/*
 * Private
 */
OwnNode_t DBG::_get_possible_extension(OwnNode_t node_id)
{
    Pairedendinformation_t node_pi = getPairedEndInformation(node_id);
    for (auto p: node_pi) {
        vector<OwnNode_t> possible_index = _get_posible_pos(p);
        if (possible_index.size() == 1)
        {
            OwnNode_t candidate = possible_index[0];
            vector<OwnNode_t> outNeighs = getNeighbor(candidate), inNeighs = getNeighbor(candidate, false);
            bool out = true;
            while (out && inNeighs.size() == 1)
            {
                candidate = inNeighs[0];
                outNeighs = getNeighbor(candidate);
                out = (outNeighs.size() == 1);
                inNeighs = getNeighbor(candidate, false);
            }
            return candidate;
        }
    }
    return NO_NEIGH;
}
void DBG::_extension(queue <UG_Node> & queue, vector<bool> & checked, vector <OwnNode_t> & unitig,
        UG_Node node, size_t & higher_out, size_t & zero, size_t & higher_in, size_t & zero_nopair,
        size_t & repetitive_call, vector<vector<OwnNode_t>> & unitigs, vector<bool> & local_checked)
{
    checked[node._id] = true;
    local_checked[node._id] = true;
    unitig.push_back(_g_nodes[node._id]._val);
    OwnNode_t parent = NO_NEIGH;
    vector<OwnNode_t> neighs = getNeighbor(node._id);
    bool extension = true;
    //cout << "Node: "<<node._id<<" "<<node._val<<endl;
    if (neighs.size() > 1 && in_degree(node._id) == 0)
    {
        for (auto i:neighs)
        {
            if (!checked[i])
            {
                _extension(queue, checked, unitig, _g_nodes[i], higher_out, zero, higher_in, zero_nopair, repetitive_call, unitigs, local_checked);
            }
        }
        return;
    }
    while (neighs.size() == 1)
    {
        parent = node._id;
        node = _g_nodes[neighs[0]];
        //cout << " "<<node._id<<" "<<node._val<<endl;
        if (in_degree(node._id) == 1 && !local_checked[node._id]) {
            unitig.push_back(node._val);
            local_checked[node._id] = true;
            neighs = getNeighbor(node._id);
        /*
         * REVISAR
         */
        /*if (!local_checked[neighs[0]])
        {
            unitig.push_back(_g_nodes[neighs[0]]._val);
            local_checked[neighs[0]] = true;
            neighs = getNeighbor(neighs[0]);*/
        } else {
            break;
        }
    }
    if (neighs.size() == 0) {
        //cout <<" Node with no neighbors: "<<node._id<<" "<<node._val<<endl;
        if (getPairedEndInformation(node._id).empty())
            zero_nopair++;
        OwnNode_t candidate = _get_possible_extension(node._id);
        if (candidate != NO_NEIGH ) {
            if (!local_checked[_g_nodes[candidate]._id])
            {
                //cout << "Candidato seleccionado: " << candidate <<" Id: "<<_g_nodes[candidate]._id<< endl;
                _extension(queue, checked, unitig, _g_nodes[candidate],
                           higher_out, zero, higher_in, zero_nopair, repetitive_call, unitigs, local_checked);
                return;
            }
        } else
            cout << "No hay candidato possible!"<<endl;
        zero++;
    }
    if (neighs.size() > 1 && in_degree(node._id) > 1)
    {
        (neighs.size()>1)?higher_out++:higher_in++;
        /*
         * Cuidado que esto es proclive a errores.
         */
        bool control = (parent != NO_NEIGH);
        if (control)
        {
            OwnNode_t node_pivot = _g_nodes[neighs[0]]._val;
            for (auto i: neighs)
            {
                Pairedendinformation_t nei_pi = getPairedEndInformation(i), node_pi = getPairedEndInformation(parent);
                Pairedendinformation_t union_pi = intersect(nei_pi, node_pi);
                if (union_pi.empty())
                    control = false;
                if (_g_nodes[i]._val != node_pivot)
                    control = false;
            }
        }
        if (control)
        {
            for (auto i:neighs)
            {
                if (!checked[i])
                {
                    _extension(queue, checked, unitig, _g_nodes[i],
                            higher_out, zero, higher_in, zero_nopair, repetitive_call, unitigs, local_checked);
                }
            }
            return;
        } else {
            /*cout << "Caso: "<< neighs.size()<<" Paired-end: "<<getPairedEndInformation(node._id).size()<<endl;
            vector<OwnNode_t> out_neighs = getNeighbor(node._id), in_neighs = getNeighbor(node._id, false);
            for (auto n:out_neighs)
                cout << "OUT: "<<n<<" "<<_g_nodes[n]._val<<" ";
            cout << endl;
            for (auto n: in_neighs)
                cout << "IN: "<<n<<" "<<_g_nodes[n]._val<<" ";
            cout << endl;*/
            for (auto i:neighs)
            {
                if (!checked[i])
                    queue.push(_g_nodes[i]);
            }
        }
    }
    unitigs.push_back(unitig);
    cout << "LLego!"<<endl;
}
vector<DBG::UG_Node> DBG::_get_starting_nodes()
{
    vector<UG_Node> suspicious_nodes;
    for (auto i: _g_nodes)
    {
        if (!in_degree(i._id))
            suspicious_nodes.push_back(i);
    }
    cout << "Suspicious starts: "<<suspicious_nodes.size()<<endl;
    return suspicious_nodes;
}
bool DBG::_reach(OwnNode_t i, OwnNode_t j)
{
    return (_reachability[i].find(j) != _reachability[i].end()) || (_reachability[j].find(i) != _reachability[j].end());
}
void DBG::_print_reachability()
{
    cout << "Exporting reachability information..."<<endl;
    for (size_t i = 0; i < _reachability.size(); ++i)
    {
        cout << "Node: "<<i<<": ";
        for (auto n: _reachability[i])
            cout << " "<<n<<", ";
        cout <<endl;
    }
}
void DBG::_complete_reach_matrix()
{
    std::function<void(size_t, size_t,size_t, size_t&)> __reach = [this, &__reach](size_t  src, size_t target,
            size_t steps, size_t &branches)->void{
        if (steps++ == MAX_PATH)
            return;
        if (branches == MAX_BRANCHES)
            return;
        //cout << src<<" - "<<target<<" - "<<steps<<" - "<<branches<<endl;
        if (src != target)
            this->_reachability[src].emplace(target);
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
    //_print_reachability();
}
/*
 * Show information
 */
string DBG::toString()
{
    return string("Not supported yet");
}

void DBG::stats()
{
    float neighs = 0.0, pair_info_avg = 0.0;
    size_t num_nodes_no_pair = 0, num_nodes_pair = 0, no_pair_no_neigh = 0;
    for (auto v:_g_nodes) {
        neighs += _g_edges[v._id].size();
        if (_g_nodes[v._id].get_paired_information().empty()) {
            if (v._id < (_g_nodes.size() / 2)) {
                num_nodes_no_pair++;
                auto neighs = getNeighbor(v._id, false);
                if (!neighs.size())
                    no_pair_no_neigh++;
            }
        }else {
            pair_info_avg += _g_nodes[v._id].get_paired_information().size();
            num_nodes_pair++;
        }
    }
    cout << "Avg neighs: "<<(neighs/_g_nodes.size())<<endl;
    cout << "Avg paired end information: "<<(pair_info_avg/num_nodes_pair)<<endl;
    cout << "Number of nodes with no paired-end information "<<num_nodes_no_pair<<" out of "<<_g_nodes.size()<<endl;
    cout << "No pair and no neighbors: "<<no_pair_no_neigh<<endl;
}

void DBG::_subsane()
{
    for (auto v:_g_nodes)
    {

    }
}

void DBG::print()
{
    for (auto v:_g_nodes)
    {
        cout << "Node: "<<v._id<<"/"<<v._val<<" Number of neighbors: "<<_g_edges[v._id].size()<<" - ";
        for (auto n:_g_edges[v._id])
            cout << " "<<n;
        cout << endl;
        cout << "Paired-end information: ";
        for (auto n: _g_nodes[v._id].get_paired_information())
            cout << " "<<n<<"/"<<_g_nodes[v._id].get_abundance(n);
        cout <<endl;
    }
    cout << "End print"<<endl;
}

/*
 * Undirected Graphs
 */
OwnNode_t UG::addVertex(OwnNode_t node)
{
    _g_nodes.push_back(Node(_num_vertex, node));
    _g_edges.push_back(unordered_set<OwnNode_t>());
    return _num_vertex++;
}

OwnNode_t UG::addVertexWithAbundances(OwnNode_t node, size_t ab)
{
    _g_nodes.push_back(Node(_num_vertex,node,ab));
    _g_edges.push_back(unordered_set<OwnNode_t>());
    return _num_vertex++;
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
    _g_edges[i].emplace(j);
    _g_edges[j].emplace(i);
}

unordered_set<OwnNode_t > UG::getNeighbors(OwnNode_t i)
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

priority_queue<pair<size_t,vector<OwnNode_t>>> UG::report_maximal_cliques(unordered_map<size_t, size_t> reverse_map,
        DBG & dbg)
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
    auto __get_maximal_clique = [this, &already_checked, &sum, &reverse_map, &dbg](OwnNode_t vertex,const int maxCliqueSize, size_t * score) {
        vector<OwnNode_t> clique;
        clique.emplace_back(vertex);
        //Avoid repetitions
        already_checked.emplace(vertex);
        sum |= 1 << _g_nodes[vertex]._id;

        set<OwnNode_t> candidateNeighbors;

        unordered_set<OwnNode_t> visited;
        visited.emplace(vertex);

        for (auto n:getNeighbors(vertex))
            candidateNeighbors.emplace(n);
        set<OwnNode_t> tmp;

        auto highestDegVert = vertex;

        while (!candidateNeighbors.empty()) {
            Pairedendinformation_t pivot = dbg.getPairedEndInformation(reverse_map[highestDegVert]);
            /*auto highestDegNeighborIt = std::max_element(candidateNeighbors.begin(), candidateNeighbors.end(), [this](const OwnNode_t &lhs, const OwnNode_t &rhs) {
                if ((degree(lhs)) == (degree(rhs)))
                    return _g_nodes[lhs]._id > _g_nodes[rhs]._id;
                return (degree(lhs)) < (degree(rhs));
            });*/
            auto highestDegNeighborIt = std::max_element(candidateNeighbors.begin(), candidateNeighbors.end(),
                    [this, & pivot, &dbg, &reverse_map](const OwnNode_t &lhs, const OwnNode_t &rhs) {
                auto inter1 = intersect(pivot, dbg.getPairedEndInformation(reverse_map[lhs]))
                        , inter2 = intersect(pivot, dbg.getPairedEndInformation(reverse_map[rhs]));
                if ((inter1.size()) == (inter2.size()))
                    return _g_nodes[lhs]._id > _g_nodes[rhs]._id;
                return (inter1.size()) < (inter2.size());
            });
            highestDegVert = *highestDegNeighborIt;
            clique.emplace_back(highestDegVert);
            already_checked.emplace(highestDegVert);
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
        /*for (auto c:clique)
            cout << c << ", ";
        cout << endl;*/
        return clique;
    };
    for (size_t vertex = 0; vertex < _g_nodes.size(); vertex++)
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