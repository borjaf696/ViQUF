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
            _g_nodes_frequency = vector<float>(_num_nodes, 0);
            _g_edges = vector<vector<OwnNode_t>>(_num_nodes, vector<size_t>());
            _g_edges_reads = vector<vector<OwnNode_t>>(_num_nodes, vector<size_t>());
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
                    _g_edges_reads[node].push_back(0);
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
}

size_t DBG::edges()
{
    return _numedges;
}


size_t DBG::vertices()
{
    return _g_nodes.size();
}

void DBG::set_length(OwnNode_t n, size_t length)
{
    _g_nodes[n].setLength(length);
}

size_t DBG::getLength(OwnNode_t n)
{
    return _g_nodes[n]._length;
}

vector<size_t> DBG::getNeighbor(size_t v , bool dir, size_t filter)
{
    if (v >= _g_nodes.size())
        throw std::invalid_argument( "Out of bounds");
    if (filter != NO_NEIGH)
    {
        vector<size_t> output, tmp;
        if (dir)
            tmp = _g_edges[v];
        else
            tmp = _g_in_edges[v];
        for (auto n: tmp) {
            if (n != filter)
                output.push_back(n);
        }
        return output;
    }
    if (dir)
        return _g_edges[v];
    return _g_in_edges[v];
}

float DBG::getFreqNode(OwnNode_t node)
{
    return _g_nodes_frequency[node];
}

unordered_set<size_t> DBG::getPairedEndInformation(size_t v)
{
    if (v >= _g_nodes.size())
        throw std::invalid_argument( "Out of bounds");
    return _g_nodes[v].get_paired_information();
}

float DBG::getAbundance(OwnNode_t u1)
{
    return _g_nodes[u1]._abundance;
}

size_t DBG::getAbundance(size_t node, size_t pair)
{
    if (node >= _g_nodes.size())
        throw std::invalid_argument( "Out of bounds");
    return _g_nodes[node].get_abundance(pair);
}

void DBG::addPair(size_t v, size_t p, bool forward, bool forward_right)
{
    /*if (v == 1190)
        cout <<"V: "<<v<<" P: "<<p<<" Forward: "<<forward<<" Forward_Right: "<<forward_right<<endl;*/
    _g_nodes[v].add_paired_information(p);
}

size_t DBG::addPair(size_t v, Pairedendinformation_t pairInformation)
{
    _g_nodes[v].add_paired_information(pairInformation);
    return pairInformation.size();
}

void DBG::build_process_cliques(DBG & apdbg,
        const vector<string> & sequence_map, size_t num_unitigs_fw)
{
    cout << "Number of nodes (in the AdBG): "<< _num_nodes - 1<<endl;
    bool show = false, turn = true;
    size_t empty_pairs = 0;
    OwnNode_t node_test = 742;
    string progress_bar;
    for (size_t i = 0; i < _num_nodes; ++i)
    {
        vector<size_t> local_histogram = vector<size_t>(10000);
        //show = (_g_nodes[i]._id == node_test);
        show = (_g_nodes[i]._val == node_test);
        Pairedendinformation_t node_pei = getPairedEndInformation(i);
        UG_Node parent(_g_nodes[i]._val, node_pei);
        //cout << " Node: "<<i<<" "<<_g_nodes[i]._val<<endl;
        vector<OwnNode_t> neighbors = getNeighbor(i);
        if (neighbors.size() == 0)
        {
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
            Pairedendinformation_t neigh_pei = getPairedEndInformation(n), new_node_pei = node_pei;
            Pairedendinformation_t neigh_pei_copy = neigh_pei;
            bool content = (neigh_pei.size() > node_pei.size())?is_subset(node_pei,neigh_pei):is_subset(neigh_pei,node_pei);
            show = (show || (_g_nodes[n]._val == node_test));
            if (show) {
                turn = _g_nodes[i]._val >= num_unitigs_fw/2;
                cout << "Node: "<<_g_nodes[i]._val<<" "<<
                sequence_map.at((turn) ? _g_nodes[i]._val - num_unitigs_fw / 2 : _g_nodes[i]._val)<<" "<<out_degree(_g_nodes[i]._id)<<endl;
                cout << "Size Node (paired-end): " << node_pei.size() << " Neighbors: " << getNeighbor(i).size() << " Indegree: "<<in_degree(_g_nodes[n]._id)<<endl;
                cout << "Paired-end information:"<<endl;
                for (auto n: node_pei)
                    cout <<" "<<n;
                cout <<endl;
                turn = _g_nodes[n]._val >= num_unitigs_fw/2;
                cout << " Neigh: " << n << " " << _g_nodes[n]._val << " "
                     << sequence_map.at((turn) ? _g_nodes[n]._val - num_unitigs_fw / 2 : _g_nodes[n]._val)<<" " <<neigh_pei.size()<< endl;
                cout << "Paired-end information: "<<endl;
                for (auto n: neigh_pei)
                    cout << " "<<n;
                cout << endl;
            }
            Pairedendinformation_t union_pei = union_set(new_node_pei, neigh_pei_copy);
            if (show)
            {
                cout << "Union original: "<<endl;
                display_set(union_pei);
            }
            if (COMPLETE != 0)
            {
                vector<OwnNode_t> in_neighs = getNeighbor(_g_nodes[n]._id, false, _g_nodes[i]._id),
                    out_neighs = getNeighbor(_g_nodes[i]._id, true, _g_nodes[n]._id);
                Pairedendinformation_t extra_neigh, extra_node;
                for (auto j:in_neighs)
                    extra_neigh = union_set(extra_neigh,getPairedEndInformation(_g_nodes[j]._id));
                for (auto j: out_neighs)
                    extra_node = union_set(extra_node,getPairedEndInformation(_g_nodes[j]._id));
                if (show) {
                    cout << " extra neigh (informacion from extra parents): " << endl;
                    display_set(extra_neigh);
                    cout << " extra node (informacion from extra sons): " << endl;
                    display_set(extra_node);
                }
                new_node_pei = minus_set(new_node_pei, minus_set(extra_neigh,neigh_pei));
                neigh_pei_copy = minus_set(neigh_pei_copy, minus_set(extra_node, node_pei));
            }
            union_pei = union_set(new_node_pei, neigh_pei_copy);
            if (show)
            {
                cout << "Union post-procesado: "<<endl;
                display_set(union_pei);
            }
            /*if (show)
                cout << "Size (paired-end): "<<neigh_pei.size()<<endl;*/
            /*if (_g_nodes[n].get_paired_information().empty()) {
                UG_Node son(_g_nodes[n]._val, getPairedEndInformation(n));
                apdbg.addNode(parent, son);
                continue;
            }*/
            UG local_graph, local_graph_2 = UG(true);
            unordered_map<size_t, size_t> traslation_map, reverse_map;
            unordered_map<size_t, float> frequencies_map;
            for (auto pn: union_pei){
                local_histogram[std::min((size_t) 999,getAbundance(n,pn) + getAbundance(i,pn))]++;
                size_t id = local_graph.addVertex(pn);
                if (show)
                    cout <<"Pn: "<<pn<<" Abundance: "<<getAbundance(pn)<<endl;
                local_graph_2.addVertexWithAbundances(pn, getAbundance(pn));
                traslation_map[pn] = id;
                frequencies_map[id] = getAbundance(pn);
                reverse_map[id] = pn;
            }
            /*vector<vector<bool>> reachs = vector<vector<bool>>(union_pei.size(),vector<bool>(union_pei.size(),false))
                    , reached_by = vector<vector<bool>>(union_pei.size(),vector<bool>(union_pei.size(),false));*/
            vector<size_t> reached_by = vector<size_t>(union_pei.size(), 0);
            std::function<void(size_t, size_t,size_t, size_t&, unordered_set<OwnNode_t> &)> __reach = [this, &__reach,
                                                                           &traslation_map,
                                                                           &local_graph_2,
                                                                           &union_pei,&show](size_t  src, size_t target,
                                                                                           size_t steps, size_t &branches
                                                                                           , unordered_set<OwnNode_t> & checked)->void{
                checked.emplace(target);
                bool swap = false;
                unordered_set<OwnNode_t> n_checked;
                size_t new_steps = 0, new_branches = 0;
                if (src != target && (in(union_pei, target))) {
                    if (show)
                        cout <<"¡Change! Src: "<<src<<" Target: "<<target<<endl;
                    local_graph_2.addEdge(traslation_map[src], traslation_map[target]);
                    src = target;
                    n_checked = checked;
                    swap = true;
                }
                steps += getLength(target);
                if (steps >= MAX_PATH)
                    return;
                if (branches == MAX_BRANCHES)
                    return;
                auto neighbors = getNeighbor(target);
                branches += (neighbors.size() > 1)?1:0;
                for (auto neigh: neighbors) {
                    if (checked.find(neigh) == checked.end())
                    {
                        if (!swap)
                            __reach(src, neigh, steps, branches, checked);
                        else
                            __reach(src, neigh, new_steps, new_branches, n_checked);
                    }
                }};
            for (auto pn: union_pei)
            {
                for (auto tn: union_pei)
                {
                    if (_reach_e(pn,tn))
                        reached_by[traslation_map[tn]]++;
                    if (_reach(pn,tn)) {
                        local_graph.addEdge(traslation_map[pn], traslation_map[tn]);
                    } else {
                        /*if (show) {
                            if (pn != tn) {
                                cout << "pn: " << pn << " tn: " << tn << " reach: " << _reach(pn, tn) << endl;
                                if (pn < num_unitigs_fw / 2)
                                    cout << "Unitig: " << sequence_map.at(pn) << endl;
                                if (tn < num_unitigs_fw / 2)
                                    cout << "Unitig2: " << sequence_map.at(tn) << endl;
                            }
                        }*/
                    }
                    /*if (show) {
                        cout << "pn: " << pn << " tn: " << tn << " reach: " << _reach(pn, tn) << endl;
                        cout << "Node (g): "<<traslation_map[pn]<<" - "<<traslation_map[tn]<<endl;
                    }*/
                }
            }
            for (auto pn:union_pei)
            {
                size_t steps = 0, branches = 0;
                if (!reached_by[traslation_map[pn]])
                {
                    unordered_set<OwnNode_t> checked;
                    __reach(pn, pn, steps, branches, checked);
                }
            }
            if (show) {
                cout << "Graph: "<<endl;
                local_graph.print();
                cout << "DAG: "<<endl;
                local_graph_2.print();
                cout << "Min cost flow: "<<endl;
                local_graph_2.report_min_cost_flow(reached_by, show);
            }
            priority_queue<pair<size_t,vector<OwnNode_t>>>
                    cliques = local_graph.report_maximal_cliques(reverse_map, *this, frequencies_map, show);
            vector<vector<OwnNode_t>> cliques_store;
            unordered_set<OwnNode_t> nodes_checked;
            float max_size = 0.0;
            if (show)
                cout << " Cliques size: "<<cliques.size()<<endl;
            vector<pair<UG_Node,UG_Node>> potential_nodes;
            while (!cliques.empty() && nodes_checked.size() < union_pei.size())
            {
                pair <size_t, vector<OwnNode_t>> top_click = cliques.top();
                cliques.pop();
                max_size = (max_size < top_click.second.size()) ? top_click.second.size() : max_size;
                cliques_store.push_back(top_click.second);
            }
            vector<bool> offset_erase, offset_node, offset_neigh;
            for (auto top_click: cliques_store)
            {
                if (top_click.size() <  SIZE_RELATION*max_size)
                    continue;
                if (show)
                    cout << "New Clique:" << endl;
                Pairedendinformation_t clique, p_info_a, p_info_b;
                bool in_a = false, in_b = false, isa = true, isb = true;
                for (auto node: top_click)
                {
                    nodes_checked.emplace(node);
                    UG::UG_Node real_node = local_graph[node];
                    clique.emplace(real_node._val);
                    if (show) {
                        bool reverse = false;
                        size_t u = real_node._val;
                        if (u >= (num_unitigs_fw/2))
                        {
                            u -= (num_unitigs_fw/2);
                            reverse = true;
                        }
                        string seq = sequence_map.at(u);
                        //seq = ((reverse)?Sequence(seq.c_str()).getRevcomp():seq);
                        cout << node << ":" << real_node._val<<":"<<_g_nodes[real_node._val]._abundance<<":"<<real_node._id<< " " << seq << " "<<reverse<<" ";
                    }
                    if (node_pei.find(real_node._val) != node_pei.end()) {
                        in_a = true;
                        p_info_a.emplace(real_node._val);
                    } else {
                        isa = false;
                        if (show) {
                            size_t u = real_node._val;
                            if (u >= (num_unitigs_fw/2))
                                u -= (num_unitigs_fw/2);
                            string seq = sequence_map.at(u);
                        }
                    }
                    if (neigh_pei.find(real_node._val) != neigh_pei.end()) {
                        in_b = true;
                        p_info_b.emplace(real_node._val);
                    } else
                        isb = false;
                }
                if (show) {
                    cout << endl;
                    cout << "ISA: "<<isa<<" ISB: "<<isb<<endl;
                }
                bool check_both = true;
                if (!node_pei.empty() and !neigh_pei.empty())
                    check_both = false;
                if (check_both or (in_a && in_b))
                {
                    if ((in_a && in_b) && ((isa || isb) & !(isa && isb)) && !content) {
                        offset_erase.push_back(true);
                        offset_node.push_back(isa);
                        offset_neigh.push_back(isb);
                    } else
                        offset_erase.push_back(false);
                    UG_Node parent_l(_g_nodes[i]._val, p_info_a), son_l(_g_nodes[n]._val,p_info_b);
                    potential_nodes.push_back(pair<UG_Node,UG_Node>(parent_l,son_l));
                }
            }
            /*
             * Procesar las paired-end y eliminar los subconjuntos
             */
            if (offset_erase.size() > 1) { //Esto es una nhapa de cuidado
                size_t num_removes = 0;
                for (size_t it = 0; it < offset_erase.size(); ++it) {
                    if (offset_erase[it]) {
                        Pairedendinformation_t p1 = (offset_node[it]) ? potential_nodes[it].second._paired_info :
                                                    potential_nodes[it].first._paired_info;
                        for (size_t it_2 = 1; it_2 < potential_nodes.size(); ++it_2) {
                            Pairedendinformation_t p2 = (offset_node[it]) ? potential_nodes[it_2].second._paired_info :
                                                        potential_nodes[it_2].first._paired_info;
                            if (p1 != p2) {
                                offset_erase[it] = is_subset(p1, p2);
                            }
                        }
                    }
                }
                num_removes = 0;
                for (size_t pos = 0; pos < offset_erase.size(); pos++) {
                    if (offset_erase[pos]) {
                        potential_nodes.erase(potential_nodes.begin() + pos - num_removes++);
                    }
                }
            }
            for (auto s: potential_nodes) {
                if (show)
                {
                    s.first.print_info();
                    s.second.print_info();
                }
                apdbg.addNode(s.first, s.second);
            }
            if (show) {
                cout << "Finally cliques to add: " << potential_nodes.size() << endl;
            }
            if (potential_nodes.size() == 0)
            {
                Pairedendinformation_t empty;
                UG_Node parent_l(_g_nodes[i]._val, empty), son_l(_g_nodes[n]._val,empty);
                apdbg.addNode(parent_l, son_l);
            }
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

vector<vector<OwnNode_t>> DBG::export_unitigs(bool show)
{
    vector<vector<OwnNode_t>> unitigs;
    vector<vector<OwnNode_t>> vector_unitigs(_g_nodes.size());
    vector<UG_Node> starting_points = _get_starting_nodes();
    queue<UG_Node> extensible_points;
    vector<bool> checked(_g_nodes.size(), false);
    size_t higher_out = 0, higher_in = 0, zero = 0, zero_nopair = 0, solved = 0;
    for (auto i: starting_points) {
        extensible_points.push(i);
    }
    cout << "Retrieving unitigs"<<endl;
    while(!extensible_points.empty())
    {
        UG_Node node = extensible_points.front();
        extensible_points.pop();
        if (checked[node._id])
            continue;
        vector<OwnNode_t> current_unitig;
        vector<bool> local_checked(_g_nodes.size(), false);
        size_t repetitive_call = 0;
        current_unitig.push_back(_g_nodes[node._id]._val);
        _extension(extensible_points, checked, current_unitig, node,
                higher_out, zero, higher_in, zero_nopair, solved,
                repetitive_call, unitigs, local_checked, show);
        //unitigs.push_back(current_unitig);
    }
    size_t not_checked = 0;
    for (size_t i = 0; i < checked.size(); ++i)
    {
        if (!checked[i]) {
            not_checked++;
        }
    }
    cout << "Higher out degree: "<<higher_out<<endl;
    cout << "Higher in degree: "<<higher_in<<endl;
    cout << "Zero case: "<<zero<<endl;
    cout << "Zero no pair case: "<<zero_nopair<<endl;
    cout << "Solved: "<<solved<<endl;
    cout << "Not checked nodes: "<<not_checked<<endl;
    return unitigs;
}

/*
 * Private
 */
void DBG::post_process_pairs(const vector<string> & sequence_map)
{
    for (size_t i = 0; i < _g_nodes.size(); ++i)
        _g_nodes[i].post_process_pairs(sequence_map);
}
size_t DBG::_find_edge(OwnNode_t u1, OwnNode_t u2, bool direction)
{
    if (direction)
        for (size_t i = 0; i < _g_edges[u1].size();++i)
        {
            if (_g_edges[u1][i] == u2)
                return i;
        }
    else
        for (size_t i = 0; i < _g_in_edges[u2].size();++i)
        {
            if (_g_in_edges[u2][i] == u1)
                return i;
        }
    cout << "U1: " << u1 << " U2: " << u2 << endl;
    exit(1);
    return NO_NEIGH;
}

void DBG::add_read(OwnNode_t u1, size_t quantity)
{
    _g_nodes[u1].increase_abundance(quantity);
}

void DBG::add_read(OwnNode_t u1, OwnNode_t u2, size_t quantity, bool direction)
{
    size_t pos = _find_edge(u1,u2);
    OwnNode_t node_increment = (direction)?u1:u2;
    _g_nodes[node_increment].increase_abundance(quantity);
    if (pos != NO_NEIGH)
    {
        _g_edges_reads[u1][pos]++;
        /*if (_g_edges_reads[u1][pos] > _g_nodes_frequency[u1])
            _g_nodes_frequency[u1] = _g_edges_reads[u1][pos];*/
        if (_g_edges_reads[u1][pos] > _g_nodes_frequency[u2])
            _g_nodes_frequency[u2] = _g_edges_reads[u1][pos];
    }
}

OwnNode_t DBG::_get_possible_extension(OwnNode_t node_id)
{
    Pairedendinformation_t node_pi = getPairedEndInformation(node_id);
    vector<size_t> curr_length, nodes;
    for (auto p: node_pi) {
        size_t curr_length_var = 0;
        vector<OwnNode_t> possible_index = _get_posible_pos(p);
        if (possible_index.size() == 1)
        {
            OwnNode_t candidate = possible_index[0], best_choice = possible_index[0];
            vector<OwnNode_t> outNeighs = getNeighbor(candidate), inNeighs = getNeighbor(candidate, false);
            bool out = true;
            while (out && inNeighs.size() == 1)
            {
                curr_length_var++;
                candidate = inNeighs[0];
                outNeighs = getNeighbor(candidate);
                out = (outNeighs.size() == 1);
                inNeighs = getNeighbor(candidate, false);
                if (out && inNeighs.size() == 1)
                    best_choice = candidate;
            }
            nodes.push_back(best_choice);
            curr_length.push_back(curr_length_var);
        }
    }
    if (nodes.size() > 0)
    {
        OwnNode_t candidate = nodes[0];
        size_t cur_max = curr_length[0];
        for (size_t i = 0; i < curr_length.size(); ++i)
            candidate = (cur_max < curr_length[i]) ? nodes[i] : candidate;
        if (cur_max > 0)
            return candidate;
    }
    return NO_NEIGH;
}
void DBG::_extension(queue <UG_Node> & queue, vector<bool> & checked, vector <OwnNode_t> & unitig,
        UG_Node node, size_t & higher_out, size_t & zero, size_t & higher_in, size_t & zero_nopair, size_t &solved,
        size_t & repetitive_call, vector<vector<OwnNode_t>> & unitigs, vector<bool> & local_checked, bool show)
{
    /*
     * Revisar punto por punto la insercion de los unitigs
     */
    checked[node._id] = true;
    local_checked[node._id] = true;
    OwnNode_t parent = NO_NEIGH;
    vector<OwnNode_t> neighs = getNeighbor(node._id);
    bool extension = true;
    //cout << "Node: "<<node._id<<" "<<node._val<<" OutDegree: "<<neighs.size()<<" Indegree: "<<in_degree(node._id)<<endl;
    if (neighs.size() > 1)
    {
        for (auto i: neighs)
        {
            if (!checked[_g_nodes[i]._id])
            {
                if (out_degree(i) == 1) {
                    vector <OwnNode_t> uni_save = unitig;
                    uni_save.push_back(_g_nodes[i]._val);
                    _extension(queue, checked, uni_save, _g_nodes[i],
                               higher_out, zero, higher_in, zero_nopair, solved,
                               repetitive_call, unitigs, local_checked);
                }else {
                    vector <OwnNode_t> uni_save = unitig;
                    uni_save.push_back(_g_nodes[i]._val);
                    unitigs.push_back(uni_save);
                    queue.push(_g_nodes[i]);
                }
            }
        }
    }
    while (neighs.size() == 1)
    {
        parent = node._id;
        node = _g_nodes[neighs[0]];
        neighs = getNeighbor(node._id);
        //cout << " "<<node._id<<" "<<node._val<<" "<<neighs.size()<<" "<<in_degree(node._id)<<endl;
        unitig.push_back(node._val);
        if (in_degree(node._id) == 1 && neighs.size() == 1)
            checked[node._id] = true;
        if (in_degree(node._id) != 1)
            break;
    }
    neighs = getNeighbor(node._id);
    if (neighs.size() == 0) {
        if (show && unitig.size() > 3) {
            cout << "Unitig: "<<endl;
            for (auto u:unitig)
                cout << " "<<u;
            cout << endl;
            cout << "Output = 0" << endl;
            cout << "Node: " << node._id << " " << _g_nodes[node._id]._val << endl;
            cout << "Paired-end: "<<endl;
            for (auto n: _g_nodes[node._id]._paired_info)
                cout << " "<<n;
            cout << endl;
            cin.get();
        }
        checked[node._id] = true;
        //cout <<" Node with no neighbors: "<<node._id<<" "<<node._val<<endl;
        if (getPairedEndInformation(node._id).empty())
        {
            /*for (size_t i = unitig.size() - 1; i > 0; --i)*/
            zero_nopair++;
        }
        OwnNode_t candidate = NO_NEIGH;//_get_possible_extension(node._id);
        if (candidate != NO_NEIGH )
        {
            if (!local_checked[_g_nodes[candidate]._id])
            {
                //cout << "Candidato seleccionado: " << candidate <<" Id: "<<_g_nodes[candidate]._id<< endl;
                /*_extension(queue, checked, unitig, _g_nodes[candidate],
                           higher_out, zero, higher_in, zero_nopair,solved,
                           repetitive_call, unitigs, local_checked);*/
                unitigs.push_back(unitig);
                return;
            }
        }
        zero++;
    }
    if (neighs.size() > 1 && in_degree(node._id) > 1)
    {
        (neighs.size()>1)?higher_out++:higher_in++;
        bool control = (parent != NO_NEIGH);
        //cout << "Control: "<<control<<endl;
        if (control)
        {
            OwnNode_t node_pivot = _g_nodes[neighs[0]]._val;
            for (auto i: neighs)
            {
                Pairedendinformation_t nei_pi = getPairedEndInformation(i), node_pi = getPairedEndInformation(parent);
                Pairedendinformation_t intersection_pi = intersect(nei_pi, node_pi);
                if (intersection_pi.empty())
                    control = false;
                if (_g_nodes[i]._val != node_pivot)
                    control = false;
            }
        }
        if (control)
        {
            if (!checked[node._id])
            {
                cout << "Caso control!"<<endl;
                //for ()
                checked[node._id] = true;
                for (auto i:neighs) {
                    /*vector<OwnNode_t> unitig_saved = unitig;
                    _extension(queue, checked, unitig, _g_nodes[i],
                            higher_out, zero, higher_in, zero_nopair,solved,
                            repetitive_call, unitigs, local_checked);
                    unitig = unitig_saved;*/
                    queue.push(_g_nodes[i]);
                }
            }
        } else {
            /*cout << "Caso: "<< neighs.size()<<" Paired-end: "<<getPairedEndInformation(node._id).size()<<endl;
            vector<OwnNode_t> out_neighs = getNeighbor(node._id), in_neighs = getNeighbor(node._id, false);
            for (auto n:out_neighs)
                cout << "OUT: "<<n<<" "<<_g_nodes[n]._val<<" ";
            cout << endl;
            for (auto n: in_neighs)
                cout << "IN: "<<n<<" "<<_g_nodes[n]._val<<" ";
            cout << endl;*/
            if (!checked[node._id]) {
                queue.push(_g_nodes[node._id]);
                /*for (auto u: unitig)
                    cout << " "<<u;
                cout<<endl;
                for (auto n:getNeighbor(node._id))
                    cout << " "<<_g_nodes[n]._id<<" "<<_g_nodes[n]._val;
                cout << endl;
                cin.get();*/
            }
            /*for (auto i:neighs)
            {
                if (!checked[i])
                    queue.push(_g_nodes[i]);
            }*/
        }
    } else {
        if (neighs.size() > 1)
        {
            if (show && unitig.size() > 3)
            {
                cout << "Unitig: "<<endl;
                for (auto u:unitig)
                    cout << " "<<u;
                cout << endl;
                cout << "Caso out > 1" << endl;
                cout << "Node: " << node._id << " " << _g_nodes[node._id]._val << endl;
                cout << "Paired-end: "<<endl;
                for (auto n: _g_nodes[node._id]._paired_info)
                    cout << " "<<n;
                cout << endl;
                cout << "Neighbors (out):"<<endl;
                for (auto n:neighs)
                    cout << "Neigh (val): "<<_g_nodes[n]._val<<" "<<_g_nodes[n]._id<<endl;
                cout << endl;
                cin.get();
            }
            OwnNode_t node_pivot = _g_nodes[neighs[0]]._val;
            bool control = true;
            for (auto i: neighs)
            {
                if (_g_nodes[i]._val != node_pivot)
                    control = false;
            }
            if (control) {
                cout << "Caso same neigh!!!" << endl;
                for (auto i : neighs)
                {
                    vector<OwnNode_t> unitig_saved = unitig;
                    _extension(queue, checked, unitig, _g_nodes[i],
                               higher_out, zero, higher_in, zero_nopair,solved,
                               repetitive_call, unitigs, local_checked, show);
                    unitigs.push_back(unitig);
                    unitig = unitig_saved;
                }
                return;
            }
            if (!checked[node._id]) {
                checked[node._id] = true;
                for (auto i: neighs)
                    queue.push(_g_nodes[i]);
            }
            /*for (auto i:neighs) {
                if (!checked[i])
                    queue.push(_g_nodes[i]);
            }*/
        }
        if (neighs.size() == 1)
        {
            if (!checked[node._id])
                queue.push(_g_nodes[node._id]);
        }
        if (in_degree(node._id) > 1)
        {
            if (!checked[node._id]) {
                queue.push(_g_nodes[node._id]);
            }
        }
    }
    unitigs.push_back(unitig);
}
vector<DBG::UG_Node> DBG::_get_starting_nodes()
{
    vector<UG_Node> suspicious_nodes;
    for (auto i: _g_nodes)
    {
        if (in_degree(i._id) == 0)
            suspicious_nodes.push_back(i);
        else if (out_degree(i._id) > 1)
            suspicious_nodes.push_back(i);
        else if (in_degree(i._id) > 1)
            suspicious_nodes.push_back(i);
    }
    cout << "Suspicious starts: "<<suspicious_nodes.size()<<endl;
    return suspicious_nodes;
}
bool DBG::_reach_e(OwnNode_t i, OwnNode_t j)
{
    return (_reachability[i].find(j) != _reachability[i].end());
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
    std::function<void(size_t, size_t,size_t, size_t&,unordered_set<OwnNode_t> &,bool)> __reach = [this, &__reach](size_t  src, size_t target,
            size_t steps, size_t &branches, unordered_set<OwnNode_t> & checked, bool show)->void{
        if (show)
            cout <<"From: "<<src<<" target: "<<target<<" Steps: "<<steps<<" Branches: "<<branches<<endl;
        steps += getLength(target)-Parameters::get().kmerSize + 1;
        if (steps > MAX_PATH)
            return;
        if (branches == MAX_BRANCHES)
            return;
        if (src != target)
            this->_reachability[src].emplace(target);
        auto neighbors = getNeighbor(target);
        checked.emplace(target);
        branches += (neighbors.size() > 1)?1:0;
        for (auto neigh: neighbors)
        {
            if (checked.find(neigh) == checked.end())
                __reach(src, neigh, steps, branches,checked, show);
        }
    };
    auto start = chrono::steady_clock::now();
    cout << "Building reachiability matrix..."<<endl;
    for (size_t i = 0; i < _num_nodes; ++i)
    {
        bool show = false;
        unordered_set<OwnNode_t> checked;
        size_t distance = 0, branches = 0;
        if (i == INF)
            show = true;
        __reach(i,i,distance, branches,checked, show);
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

void DBG::subsane(const vector<string> & sequence_map)
{
    /*
     * Set frequencies for nodes based on frequencies of the edges
     * TODO: remove prints
     */
    for (size_t i = 0; i < _g_edges.size(); ++i)
    {
        bool show = false;
        float node_freq = _g_nodes_frequency[i];
        for (size_t j = 0; j < _g_edges[i].size(); ++j)
        {
            float neigh_freq = _g_nodes_frequency[_g_edges[i][j]];
            // ¿Max o min?
            //node_freq = max(neigh_freq, node_freq);
            if (show)
            {
                cout <<"Node: "<< _g_nodes[i]._val<<" "<<node_freq<<endl;
                cout << "Neighbor: "<<_g_nodes[_g_edges[i][j]]._val<<" Frequency: "<<neigh_freq<<endl;
                cout << "Edge frequency: "<<_g_edges_reads[i][j]<<endl;
            }
            node_freq = max(neigh_freq, node_freq);
            if (_g_edges_reads[i][j]*RATIO < node_freq)
            {
                OwnNode_t connected_node = _g_edges[i][j];
                _g_edges_reads[i].erase(_g_edges_reads[i].begin() + j);
                _g_edges[i].erase(_g_edges[i].begin()+j);
                size_t pos = _find_edge(_g_nodes[i]._id, connected_node, false);
                _g_in_edges[connected_node].erase(_g_in_edges[connected_node].begin() + pos);
            }
        }
    }
    /*
     * Change from freq to score
     */
    for (size_t i = 0; i < _g_nodes_frequency.size()/2; ++i)
    {
        _g_nodes[i]._abundance /= sequence_map[i].size();
        _g_nodes[i+sequence_map.size()]._abundance /= sequence_map[i].size();
    }
    /*
    * Getting full reachability matrix - after polishing the matrix
    */
    _complete_reach_matrix();
}

void DBG::print()
{
    for (auto v:_g_nodes)
    {
        cout << "Node: "<<v._id<<"/"<<v._val;
        cout << " Frequency: "<<_g_nodes_frequency[v._id]<<" - ";
        cout <<" Number of neighbors: "<<_g_edges[v._id].size()<<" - ";
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

OwnNode_t UG::addVertexWithAbundances(OwnNode_t node, float ab)
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
    if (!_directed) {
        _g_edges[i].emplace(j);
        _g_edges[j].emplace(i);
    }else
        _g_edges[i].emplace(j);
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
        DBG & dbg, unordered_map<size_t, float> & frequency_map, bool show)
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
    auto __get_maximal_clique = [this, &frequency_map, &already_checked, &sum, &reverse_map, &dbg, &show]
            (OwnNode_t vertex,const int maxCliqueSize, size_t * score) {
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
        float max_freq = INF;
        if (show)
            cout << "First node: "<<highestDegVert<<endl;
        while (!candidateNeighbors.empty())
        {
            Pairedendinformation_t pivot = dbg.getPairedEndInformation(reverse_map[highestDegVert]);
            float freq_pivot = frequency_map[highestDegVert];
            if (freq_pivot != 0)
                max_freq = (max_freq > freq_pivot)?freq_pivot:max_freq;
            /*auto highestDegNeighborIt = std::max_element(candidateNeighbors.begin(), candidateNeighbors.end(),
                 [this, &frequency_map, &dbg, &max_freq, &show](const OwnNode_t &lhs, const OwnNode_t &rhs) {
                     size_t freq_lhs = frequency_map[lhs], freq_rhs = frequency_map[rhs];
                     if ( freq_lhs == freq_rhs)
                         return _g_nodes[lhs]._id < _g_nodes[rhs]._id;
                     return freq_lhs < freq_rhs;
                 });*/
            /*
             * Differencia minima
             */
            auto highestDegNeighborIt = std::max_element(candidateNeighbors.begin(), candidateNeighbors.end(),
                    [this, &frequency_map, &dbg, &freq_pivot, &show](const OwnNode_t &lhs, const OwnNode_t &rhs) {
                float freq_lhs = frequency_map[lhs], freq_rhs = frequency_map[rhs];
                float diff_lhs = std::abs(int(freq_lhs - freq_pivot)), diff_rhs = std::abs(int(freq_rhs - freq_pivot));
                /*if ((degree(lhs)) == (degree(rhs)))
                    return _g_nodes[lhs]._id > _g_nodes[rhs]._id;
                return (degree(lhs)) < (degree(rhs));*/
                if ( diff_lhs == diff_rhs)
                    return _g_nodes[lhs]._id > _g_nodes[rhs]._id;
                return (diff_lhs) > (diff_rhs);
            });
            /*auto highestDegNeighborIt = std::max_element(candidateNeighbors.begin(), candidateNeighbors.end(),
                    [this, & pivot, &dbg, &reverse_map](const OwnNode_t &lhs, const OwnNode_t &rhs) {
                auto inter1 = intersect(pivot, dbg.getPairedEndInformation(reverse_map[lhs]))
                        , inter2 = intersect(pivot, dbg.getPairedEndInformation(reverse_map[rhs]));
                if ((inter1.size()) == (inter2.size()))
                    return _g_nodes[lhs]._id > _g_nodes[rhs]._id;
                return (inter1.size()) < (inter2.size());
            });*/
            highestDegVert = *highestDegNeighborIt;
            if (show)
                cout << "Seleccion: "<<highestDegVert<<" "<<max_freq<<" "<<frequency_map[highestDegVert]<<endl;
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
        /*
         * ¿Es necesario actualizar?
         * */
        /*for (auto c:clique)
            frequency_map[c] = (frequency_map[c] < max_freq)?0:frequency_map[c] - max_freq;*/
        if (show)
            cout << endl;
        return clique;
    };
    auto own_order = [this, &frequency_map] (Node i,Node j) { return (frequency_map[i._id] > frequency_map[j._id]);};
    vector<Node> aux = _g_nodes;
    std::sort(aux.begin(), aux.end(), own_order);
    for (size_t vertex = 0; vertex < _g_nodes.size(); vertex++)
    {
        sum = 0;
        size_t score = 0;
        size_t r_vertex = aux[vertex]._id;
        if (already_checked.find(r_vertex) == already_checked.end())
        {
            tmpClique = __get_maximal_clique(r_vertex, maxClique.size(), &score);
            //cout << endl;
            if ((tmpClique.size() >= CLIQUE_LIMIT) && (idCliques.find(sum) == idCliques.end())) {
                idCliques.emplace(sum);
                endCliques.push(pair < size_t, vector < OwnNode_t >> (score, tmpClique));
            }
        }
    }
    return endCliques;
}

void UG::report_min_cost_flow(const vector<size_t> & reached_by, bool show)
{
    Graph l_g;
    vector<Graph::Node> node_vect;
    Graph::NodeMap<OwnNode_t> translation_map(l_g);
    Graph::Node source = l_g.addNode(), target = l_g.addNode(), s_source = l_g.addNode(), s_target = l_g.addNode();
    size_t id_source = reached_by.size(), id_target = reached_by.size() + 1;
    // Build the source/sink and synthetic sources and sinks - to avoid infinite and insufficient
    translation_map[source] = id_source; translation_map[target] = id_target;
    translation_map[s_source] = INF; translation_map[s_target] = INF;
    /*
     * Out flow and in flow
     */
    unordered_map<OwnNode_t, float> in_flow, out_flow;
    in_flow[id_source] = 0.0;out_flow[id_source] = 0.0;
    in_flow[id_target] = 0.0;out_flow[id_target] = 0.0;
    for (auto n:_g_nodes)
    {
        Graph::Node l_node = l_g.addNode();
        translation_map[l_node] = n._id;
        node_vect.push_back(l_node);
        in_flow[n._id] = 0.0;out_flow[n._id] = 0.0;
    }
    // Add source and target
    node_vect.push_back(source);node_vect.push_back(target);
    vector<Graph::Arc> arc_vector;
    ArcMap<Capacity> capacities(l_g);
    ArcMap<Weight> weights(l_g);

    /*
     * Edge from sink to source
     */
    size_t max_capacity = 0;
    for (size_t i = 0; i < _g_edges.size(); ++i)
    {
        cout << reached_by[i]<<endl;
        cout << _g_edges[i].size()<<endl;
        Graph::Node n_1 = node_vect[i];
        if (!reached_by[i])
        {
            Graph::Arc a = l_g.addArc(source, n_1);
            capacities[a] = INF;
            weights[a] = 1;
            arc_vector.push_back(a);
            a = l_g.addArc(n_1, source);
            capacities[a] = _g_nodes[i]._abundance;
            weights[a] = 1;
            max_capacity = (max_capacity < _g_nodes[i]._abundance)?_g_nodes[i]._abundance:max_capacity;
            // Update in and out flow
            in_flow[id_source] += _g_nodes[i]._abundance;
            out_flow[i] += _g_nodes[i]._abundance;
        }
        if (_g_edges[i].size() == 0)
        {
            Graph::Arc a = l_g.addArc(n_1,target);
            capacities[a] = INF;
            weights[a] = 1;
            // Back flow
            a = l_g.addArc(target,n_1);
            capacities[a] = _g_nodes[i]._abundance;
            weights[a] = 1;
            arc_vector.push_back(a);
            // Update in and out flow
            out_flow[id_target] += _g_nodes[i]._abundance;
            in_flow[i] += _g_nodes[i]._abundance;
        }
        for (auto n: _g_edges[i])
        {
            Graph::Node n_2 = node_vect[n];
            Graph::Arc a = l_g.addArc(n_1, n_2);
            capacities[a] = INF;
            weights[a] = 1;
            arc_vector.push_back(a);
            // Back flow
            a = l_g.addArc(n_2, n_1);
            capacities[a] = _g_nodes[n]._abundance;
            weights[a] = 1;
            // Update in and out flow
            out_flow[n] += _g_nodes[n]._abundance;
            in_flow[i] += _g_nodes[n]._abundance;
        }
    }
    /*
     * Readjust flow
     */
    for (auto k_v: in_flow)
    {
        float exogenous_flow = out_flow[k_v.first] - k_v.second;
        if (exogenous_flow > 0)
        {
            Graph::Arc a = l_g.addArc(s_source, node_vect[k_v.first]);
            capacities[a] = exogenous_flow;
            weights[a] = 1;
        } else if (exogenous_flow < 0)
        {
            Graph::Arc a = l_g.addArc(s_source, node_vect[k_v.first]);
            capacities[a] = exogenous_flow;
            weights[a] = 1;
        }
    }
    /*
     * Edge from global sink to source
     */
    Graph::Node g_source = l_g.addNode(), g_target = l_g.addNode();
    translation_map[g_source] = INF; translation_map[g_target] = INF;
    Graph::Arc a = l_g.addArc(g_target, g_source);
    capacities[a] = INF;
    weights[a] = 0;
    arc_vector.push_back(a);
    ////////
    NS ns(l_g);
    ns.costMap(weights).upperMap(capacities).stSupply(source, target, max_capacity);
    ArcMap<Capacity> flows(l_g);
    NS::ProblemType status = ns.run();
    switch (status) {
        case NS::INFEASIBLE:
            cerr << "insufficient flow" << endl;
            break;
        case NS::OPTIMAL:
            ns.flowMap(flows);
            for (auto a: arc_vector)
            {
                cout << "Flow: " << ns.flow(a) << endl;
            }
            cerr << "cost=" << ns.totalCost() << endl;
            break;
        case NS::UNBOUNDED:
            cerr << "infinite flow" << endl;
            break;
        default:
            break;
    }
}