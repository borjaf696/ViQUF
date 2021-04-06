// Created by borja on 8/11/19.
//
#ifndef GATB_TRIAL_GRAPH_H
    #include "graph.h"
#endif
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

DBG::DBG(char * file)
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
            vector<string> results;
            split(line, results, ' ');
            istringstream iss(results[0]),iss2(results[1]);
            size_t _num_nodes;
            iss >> _num_nodes;
            iss2 >> _min_abundance;
            cout << "Number of nodes: "<<_num_nodes<<" Minimal depth: "<<_min_abundance<<endl;
            _g_nodes = vector<UG_Node>(_num_nodes, UG_Node());
            _g_nodes_frequency = vector<float>(_num_nodes, 0);
            _g_edges = vector<vector<OwnNode_t>>(_num_nodes, vector<size_t>());
            _g_edges_reads = vector<vector<OwnNode_t>>(_num_nodes, vector<size_t>());
            _g_in_edges = vector<vector<OwnNode_t>>(_num_nodes, vector<size_t>());
            _reachability = vector<unordered_set<size_t>>(_num_nodes, unordered_set<size_t>());
            _num_nodes = 0;_numedges = 0;
        } else {
            vector <string> results;
            split(line, results, ' ');
            istringstream iss(results[0]);
            size_t node;
            iss >> node;
            //_g_nodes.push_back(UG_Node(_num_nodes, node));
            _g_nodes[node].setId(node);
            _g_nodes[node].setVal(node);
            _num_nodes++;
            if (results.size() > 1) {
                for (size_t i = 1; i < results.size(); ++i) {
                    if (i == results.size() - 2) {
                        float abundance;
                        istringstream iss(results[i]);
                        iss >> abundance;
                        _g_nodes[node].set_abundance(abundance);
                    } else if (i == results.size() - 1){
                        size_t length;
                        istringstream iss(results[i]);
                        iss >> length;
                        _g_nodes[node].setLength(length, _min_abundance);
                    } else {
                        size_t target;
                        istringstream iss(results[i]);
                        iss >> target;
                        if (target == node)
                            continue;
                        size_t neighbor;
                        iss = istringstream(results[i]);
                        iss >> neighbor;
                        _g_edges[node].push_back(neighbor);
                        _g_edges_reads[node].push_back(0);
                        _g_in_edges[neighbor].push_back(node);
                        _numedges++;
                    }
                }
            }
        }
    }
    for (size_t i = 0; i < _g_nodes.size(); ++i)
        _g_nodes[i].assayStatus();
    cout << "Sanity (number of nodes): "<<_num_nodes<<endl;
    inFile.close();
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds (graph construction from file) : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
}

void DBG::polish(bool apdbg)
{
    std::ofstream nodes_deactivated;
    if (Parameters::get().debug)
    {
        if(!apdbg)
            nodes_deactivated.open("graphs/nodes_deactivated.txt", std::ofstream::out);
        else
            nodes_deactivated.open("graphs/apdbg_nodes_deactivated.txt", std::ofstream::out);
    }
    bool changes = true;
    while (changes) {
        changes = false;
        /*
         * Remove isolated nodes
         */
        for (size_t i = 0; i < _g_nodes.size(); ++i) {
            if (getNodeState(i)) {
                if (_g_edges[i].size() == 0 && _g_in_edges[i].size() == 0 && _g_nodes[i]._active && _g_nodes[i]._length < MIN_LENGTH_PATH) {
                    if (Parameters::get().debug)
                        nodes_deactivated << "Node: " << (_g_nodes[i]._val) << " " << (_g_nodes[i]._id) << " "
                                          << _g_nodes[i]._length
                                          << " No out/in neighs" << endl;
                    _g_nodes[i]._active = false;
                    changes = true;
                }
            }
        }
        /*
         * MOD -> Reactivating tips
         */
        //if (Parameters::get().metagenomic) {
            /*
             * Remove tips
             */
            auto check_tip = [this, & nodes_deactivated](OwnNode_t node) {
                auto neighs = getNeighbor(node);
                vector <OwnNode_t> stored_nodes(1, node);
                size_t cur_length = _g_nodes[node]._length;
                while (true) {
                    if (neighs.size() > 1)
                        return false;
                    if (cur_length > MIN_TIP_LENGTH)
                        return false;
                    if (neighs.size() == 0) {
                        for (auto n:stored_nodes) {
                            if (Parameters::get().debug)
                                nodes_deactivated << "Node: " << (_g_nodes[n]._val) << " " << (_g_nodes[n]._id)
                                                  << " Tip" << endl;
                            _g_nodes[n]._active = false;
                        }
                        return true;
                    }
                    /*cout << "Neighs size: "<<neighs.size()<<endl;
                    cout << "Neigh: "<<neighs[0]<<" "<<_g_nodes[neighs[0]]._val<<" "<<_g_nodes[neighs[0]]._length<<" "<<cur_length<<endl;*/
                    cur_length += (_g_nodes[neighs[0]]._length - Parameters::get().kmerSize + 1);
                    stored_nodes.push_back(neighs[0]);
                    neighs = getNeighbor(neighs[0]);
                }
            };
            for (size_t i = 0; i < _g_nodes.size(); ++i) {
                if (!_g_nodes[i]._active)
                    continue;
                auto neighs = getNeighbor(i);
                if (neighs.size() > 1)
                    for (auto n: neighs)
                        changes |= check_tip(n);
            }
        }
    //}
    /*
     * Ultra greedy approach by joining splits with out-neighbors = 0 to same nodes with in-neighbors = 0 by checking their pe information (subsets?)
     */
    if (apdbg)
    {
        if (Parameters::get().debug)
            nodes_deactivated << "Joining repeated unsupported nodes!"<<endl;
        for (auto node: _g_nodes)
        {
            if (out_degree(node._id) == 0)
            {
                if (Parameters::get().debug){
                    nodes_deactivated << "Node: "<<node._id<<" With val: "<<node._val<<endl;
                    nodes_deactivated << "Number of splits: "<<_split_map[node._val].size()<<" "<<in_degree(node._id)<<endl;
                    node.export_info(nodes_deactivated);
                    nodes_deactivated << "Rest: "<<endl;
                }
                for (auto split: _split_map[node._val]){
                    if (in_degree(_g_nodes[split]._id) == 0){
                        UG_Node split_node = _g_nodes[split];
                        split_node.export_info(nodes_deactivated);
                        if (is_subset(split_node._paired_info,node._paired_info) || is_subset(node._paired_info,split_node._paired_info))
                        {
                            if (Parameters::get().debug)
                                nodes_deactivated<<"Is Subset! Joining"<<endl;
                            // Pointing from the other
                            _g_edges[node._id] = _g_edges[split_node._id];
                            _g_edges_reads[node._id] = _g_edges_reads[split_node._id];
                            // Change in_edges information
                            for (auto neighs: _g_edges[split_node._id])
                                _g_in_edges[neighs].push_back(node._id);
                            _g_nodes[split]._active = false;
                        }
                    }
                }
                //cin.get();
            }
        }
    }
    // Lets ensure capacities are on point - REVISAR (Ahora mismo usamos mínimo y corregimos todo)
    if (apdbg)
    {
        if (Parameters::get().debug)
            nodes_deactivated << "Readjusting capacities!"<<endl;
        auto starting_points = _get_starting_nodes_basic();
        size_t max_capacity = INF;
        vector<OwnNode_t> cur_unitig;
        vector<float> freqs;
        stack<OwnNode_t> nodes_to_process;
        vector<bool> checked = vector<bool>(_g_nodes.size(), false);
        for (auto sp:starting_points)
            nodes_to_process.push(sp._id);
        if (Parameters::get().debug)
            nodes_deactivated << "Readjusting capacities from "<<starting_points.size()<<endl;
        while(!nodes_to_process.empty()){
            OwnNode_t node = nodes_to_process.top();
            nodes_to_process.pop();
            cur_unitig.push_back(node);
            if (out_degree(node) > 1 || in_degree(node) > 1 || out_degree(node) == 0)
            {
                if (freqs.size() > 0){
                    sort(freqs.begin(), freqs.end());
                    // Lets use median
                    max_capacity = Maths::median(freqs);
                }
                if (Parameters::get().debug){
                    nodes_deactivated << "Unitig size: "<<cur_unitig.size()<<" "<<in_degree(cur_unitig[0])<<endl;
                    nodes_deactivated << "Median capacity: "<<max_capacity<<endl;
                }
                // Changing first node connection
                if (in_degree(cur_unitig[0]) == 1 && max_capacity != INF){
                    OwnNode_t node_local = cur_unitig[0];
                    cout << "Caso indegree "<<cur_unitig.size()<<" "<<node_local<<" in_degree "<<in_degree(node_local)<<" "<<_g_in_edges[node_local].size()<<endl;
                    cout<<endl;
                    if (Parameters::get().debug)
                        nodes_deactivated << "(initial case) changing capacity for: "<<_g_in_edges[node_local][0]<<" with neigh: "<<node_local<<
                            " Original freq: "<<_g_edges_reads[_g_in_edges[node_local][0]][0]<<endl;
                    for (size_t neighs_parents = 0; neighs_parents < _g_edges[_g_in_edges[node_local][0]].size();neighs_parents++){
                        if (_g_edges[_g_in_edges[node_local][0]][neighs_parents] == node_local)
                            _g_edges_reads[_g_in_edges[node_local][0]][neighs_parents] = max_capacity;
                    }
                }
                for (size_t i = 0; i < cur_unitig.size()-1; ++i)
                {
                    OwnNode_t node_local = cur_unitig[i], neigh = cur_unitig[i+1];
                    if (_g_edges[node_local][0] == neigh && max_capacity != INF)
                    {
                        if (Parameters::get().debug)
                            nodes_deactivated << "changing capacity for: "<<node_local<<" with neigh: "<<_g_edges[node_local][0]<<
                                " Original freq: "<<_g_edges_reads[node_local][0]<<endl;
                        _g_edges_reads[node_local][0] = max_capacity;
                    }
                }
                cur_unitig.clear();
                freqs.clear();
                /*if (in_degree(node) > 1)
                    cur_unitig.push_back(node);*/
                max_capacity = INF;
            } else {
                freqs.push_back(_g_edges_reads[node][0]);
                max_capacity = (max_capacity > _g_edges_reads[node][0])?_g_edges_reads[node][0]:max_capacity;
            }
            if (!checked[node]){
                checked[node] = true;
                for (auto n: getNeighbor(node))
                    nodes_to_process.push(n);
            }
        }
    }
    if (Parameters::get().debug)
        nodes_deactivated.close();
}

size_t DBG::edges()
{
    size_t num_edges = 0;
    for (auto n:_g_nodes)
        num_edges += (n._active)?getNeighbor(n._id).size():0;
    return num_edges;
}


size_t DBG::vertices(bool actives)
{
    if (actives) {
        size_t num_nodes = 0;
        for (auto n: _g_nodes)
            num_nodes += (n._active) ? 1 : 0;
        return num_nodes;
    } else
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
            if (n != filter & _g_nodes[n]._active)
                output.push_back(n);
        }
        return output;
    }
    if (dir)
    {
        vector<size_t> output;
        for (auto n: _g_edges[v])
        {
            if (_g_nodes[n]._active)
                output.push_back(n);
        }
        return output;
    }
    return _g_in_edges[v];
}

float DBG::getFreqNode(OwnNode_t node)
{
    return _g_nodes[node]._abundance;//_g_nodes_frequency[node];
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
    /*if (v == 0 && p == 3)
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
    size_t empty_pairs = 0, nodes_process_sf = 0;
    OwnNode_t node_test = INF;
    string progress_bar;
    /*
     * Output to File
     */
    std::ofstream rejected_nodes;
    if (Parameters::get().debug) {
        rejected_nodes.open("debug/rejected_nodes.txt", std::ofstream::out );
        rejected_nodes<<"Rejected kmers:"<<endl;
    }
    /*
     * Initialize split mapping
     */
    apdbg.initialize_split_map();
    auto starting_points = _get_starting_nodes_basic();
    queue<OwnNode_t> nodes_process;
    vector<bool> checked = vector<bool>(_g_nodes.size(), false);
    for (auto starting_point: starting_points) {
        checked[starting_point._id] = true;
        nodes_process.push(starting_point._id);
    }
    while (true)
    {
        while (!nodes_process.empty())
    {
    /*for (size_t i = 0; i < _num_nodes; ++i)
    {*/
        OwnNode_t i = nodes_process.front();
        nodes_process.pop();
        nodes_process_sf++;
        if (!_g_nodes[i]._active)
            continue;
        size_t val_node = _g_nodes[i]._val, node_id = _g_nodes[i]._id;
        if (((size_t) _g_nodes[i]._abundance) < _quantile_L)
        {
            if (Parameters::get().debug)
                rejected_nodes << "Kmer: "<<_g_nodes[i]._val<<" "<<((size_t)_g_nodes[i]._abundance)<<" "
                    <<Common::return_unitig(sequence_map, _g_nodes[i]._val)<<endl;
            continue;
        }
        //show = (_g_nodes[i]._id == node_test);
        show = (val_node == node_test);
        Pairedendinformation_t node_pei;
        /*
         * Polish node pairs based on abundance
         */
        for (auto pn: getPairedEndInformation(i))
        {
            /*if (getAbundance(pn) < _quantile_L || getAbundance(pn) > _quantile_H)
                continue;*/
            node_pei.emplace(pn);
        }
        UG_Node parent(_g_nodes[i]._val, node_pei);
        vector<OwnNode_t> neighbors = getNeighbor(i);
        if (neighbors.size() == 0)
        {
            parent.setLength(_g_nodes[i]._length);
            parent.set_abundance(_g_nodes[i]._abundance);
            apdbg.addNode(parent, show);
        }
        //if (show)
            cout << " Node: "<<i<<" "<<_g_nodes[i]._val
                <<" Neighbors: "<<neighbors.size()<<" InNeighbors: "<<in_degree(i)
                <<" Nodes processed: "<<nodes_process_sf<<endl;
        if (node_pei.empty())
            empty_pairs++;
        for (size_t n_i = 0; n_i < neighbors.size(); ++n_i)
        {
            auto n = neighbors[n_i];
            if (!checked[n]) {
                checked[n] = true;
                nodes_process.push(n);
            }
            if (_g_nodes[n]._active == false)
                continue;
            auto in_parents = getNeighbor(_g_nodes[n]._id, false);
            size_t val_neigh = _g_nodes[n]._val;
            //size_t strain_freq = Maths::my_round(_g_edges_reads[i][n_i], true);
            size_t strain_freq = round(_g_edges_reads[i][n_i]);
            float strain_freq_check = strain_freq;
            /*
             * Output to File
             */
            std::ofstream basic_information, paths, rejected_pairs;
            if (Parameters::get().debug)
            {
                basic_information.open("debug/"+to_string(_g_nodes[i]._val)+"_"+to_string(_g_nodes[n]._val)+"_basic.txt", std::ofstream::out);
                paths.open("debug/"+to_string(_g_nodes[i]._val)+"_"+to_string(_g_nodes[n]._val)+"_paths.txt", std::ofstream::out);
                rejected_pairs.open("debug/"+to_string(_g_nodes[i]._val)+"_"+to_string(_g_nodes[n]._val)+"_rejected_pairs.txt", std::ofstream::out);
            }
            if (_g_nodes[n]._abundance < _quantile_L)
            {
                if (Parameters::get().debug)
                    rejected_nodes << "Kmer: "<<_g_nodes[n]._val<<" "<<((size_t)_g_nodes[n]._abundance)<<" "
                                   <<Common::return_unitig(sequence_map, _g_nodes[n]._val)<<endl;
                continue;
            }
            /*
             * Polish node pairs based on abundance
             */
            Pairedendinformation_t neigh_pei;
            for (auto pn: getPairedEndInformation(n))
            {
                /*if (getAbundance(pn) < _quantile_L || getAbundance(pn) > _quantile_H)
                {
                    if (Parameters::get().debug)
                        rejected_pairs << "Kmer (suspicious): "<<_g_nodes[pn]._val<<" "<<((size_t)_g_nodes[pn]._abundance)<<" "
                                   <<Common::return_unitig(sequence_map, _g_nodes[pn]._val)<<" NEIGH"<<endl;
                    continue;
                }*/
                neigh_pei.emplace(pn);
            }
            Pairedendinformation_t neigh_pei_copy = neigh_pei, new_node_pei = node_pei;
            
            //Maybe negligible
            bool content = (neigh_pei.size() > node_pei.size())?is_subset(node_pei,neigh_pei):is_subset(neigh_pei,node_pei);
            if (Parameters::get().debug) {
                basic_information << "Node: "<<_g_nodes[i]._val<<" ";
                basic_information << Common::return_unitig(sequence_map, _g_nodes[i]._val);
                basic_information<<" "<<out_degree(_g_nodes[i]._id)<<endl;
                basic_information << "Size Node (paired-end): " << node_pei.size()
                    << " Neighbors: " << getNeighbor(i).size() << " Indegree: "<<in_degree(_g_nodes[i]._id)<<endl;
                basic_information << "Paired-end information:"<<endl;
                for (auto n: node_pei)
                    basic_information <<" "<<n;
                basic_information <<endl;
                basic_information << " Neigh: " << n << " " << _g_nodes[n]._val << " Frequency of the join: "<<strain_freq<<" ("<<_g_edges_reads[i][n_i]<<") ";
                basic_information<<Common::return_unitig(sequence_map,_g_nodes[n]._val);
                basic_information<<" " <<neigh_pei.size()<< endl;
                basic_information << "Paired-end information: "<<endl;
                for (auto n: neigh_pei)
                    basic_information << " "<<n;
                basic_information << endl;
                basic_information << "Frequency of the strain (estimated): "<<strain_freq<<endl;
            }
            Pairedendinformation_t union_pei = union_set(new_node_pei, neigh_pei_copy);
            if (Parameters::get().debug)
            {
                basic_information << "Union original: "<<endl;
                basic_information << set_to_string(union_pei);
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
                if (Parameters::get().debug) {
                    basic_information << " extra neigh (informacion from extra parents): " << endl;
                    basic_information << set_to_string(extra_neigh);
                    basic_information << " extra node (informacion from extra sons): " << endl;
                    basic_information << set_to_string(extra_node);
                }
                new_node_pei = minus_set(new_node_pei, minus_set(extra_neigh,neigh_pei));
                neigh_pei_copy = minus_set(neigh_pei_copy, minus_set(extra_node, node_pei));
            }
            union_pei = union_set(new_node_pei, neigh_pei_copy);
            union_pei.emplace(val_node);
            union_pei.emplace(val_neigh);
            /*
            * Corregir - siempre estan insertados a dia de hoy
            */
            bool parent_in_union = in(union_pei, _g_nodes[i]._val), son_in_union = in(union_pei, _g_nodes[n]._val);
            bool cautious_execution = (parent_in_union || son_in_union);
            if (Parameters::get().debug)
            {
                basic_information << "Union post-procesado: "<<endl;
                basic_information << set_to_string(union_pei);
                basic_information << "Parent in union: "<<parent_in_union<<endl<<"Son in union: "<<son_in_union<<endl;
                basic_information << ((cautious_execution)?"GUIDED EXECUTION":"REGULAR EXECUTION")<<endl;
            }
            UG local_graph, local_graph_2 = UG(true);
            unordered_map<size_t, size_t> traslation_map, reverse_map;
            unordered_map<size_t, float> frequencies_map;
            unordered_map<size_t, bool> true_nodes;
            float max_flow = 0.0, min_flow = INF;
            /*
             * RESCORING version
             */
            vector<bool> rescoring = vector<bool>(union_pei.size(), false);
            for (auto pn: union_pei)
            {
                size_t id = local_graph.addVertex(pn);
                if (RESCORING) {
                    if (/*(pn == _g_nodes[i]._val)*/ (node_pei.find(pn) != node_pei.end() && neigh_pei.find(pn) != neigh_pei.end()))
                        rescoring[id] = true;
                }
                //float node_abundance = (float) Maths::my_round(getAbundance(pn), true);
                float node_abundance = round((float) getAbundance(pn));
                if (Parameters::get().debug)
                    basic_information << "Pair: "<<pn<<" "<<getAbundance(pn)<<" "
                        <<abs(node_abundance - strain_freq)<<" "<<node_abundance<<" Id: "<<id<<" Belongs to both: "<<rescoring[id]<<endl;
                local_graph_2.addVertexWithAbundances(pn, ((pn == val_node) || (pn == val_neigh)?strain_freq:node_abundance));
                traslation_map[pn] = id;
                true_nodes[id] = true;
                frequencies_map[id] = getAbundance(pn);
                reverse_map[id] = pn;
            }
            /*vector<vector<bool>> reachs = vector<vector<bool>>(union_pei.size(),vector<bool>(union_pei.size(),false))
                    , reached_by = vector<vector<bool>>(union_pei.size(),vector<bool>(union_pei.size(),false));*/
            vector<size_t> reached_by = vector<size_t>(union_pei.size(), 0);
            /*
             * Processing reachability information.
             * Minimum reached by - to sources (only in no cautious mode) otherwise parent/son as sources
             */
            size_t minimum_reached_by = INF;
            if (!cautious_execution)
            {
                if (Parameters::get().debug)
                    basic_information << "Getting reachability information"<<endl;
                for (auto pn: union_pei) {
                    for (auto tn: union_pei) {
                        if (tn == pn)
                            continue;
                        if (_reach_e(pn, tn)) {
                            reached_by[traslation_map[tn]]++;
                            if (Parameters::get().debug) {
                                basic_information << pn << " reachs " << tn << " Nodos: " << traslation_map[pn] << "->"
                                                  << traslation_map[tn] << endl;
                                basic_information << Common::return_unitig(sequence_map, _g_nodes[pn]._val)
                                                  << " -> " << Common::return_unitig(sequence_map, _g_nodes[tn]._val)
                                                  << endl;
                            }
                        }
                        if (_reach(pn, tn)) {
                            local_graph.addEdge(traslation_map[pn], traslation_map[tn]);
                        }
                    }
                }
                for (auto pn: union_pei)
                    minimum_reached_by = (reached_by[traslation_map[pn]] < minimum_reached_by)?reached_by[traslation_map[pn]]:minimum_reached_by;
                if (Parameters::get().debug)
                    basic_information << "Minimum reached by: "<< minimum_reached_by<<endl;
            } else {
                if (Parameters::get().debug)
                    basic_information << "Parent/Son as starting points "<<((parent_in_union)?"Parent":"Son")<<endl;
                reached_by = vector<size_t>(union_pei.size(), INF);
                minimum_reached_by = 0;
                if (parent_in_union)
                    reached_by[traslation_map[val_node]] = minimum_reached_by;
                else
                    reached_by[traslation_map[val_neigh]] = minimum_reached_by;
            }
            /*
             * Frequency distributor
             */
            auto __distribute = [this](float f0, size_t node, const vector<OwnNode_t> & neighbors, const vector<float> & in_flows, const vector<float> & offsets)
            {
                vector<OwnNode_t> neighs;
                vector<float> freqs, in_flows_local = in_flows;
                // Add the expected frequency as in_flow with no parent, we will remove it afterwards
                if (f0 != 0)
                    in_flows_local.push_back(f0);
                unordered_map<int, OwnNode_t> __traslation_map;
                __traslation_map[-1] = INF;
                size_t num_neigh = 0;
                for (size_t i = 0; i < neighbors.size(); ++i)
                {
                    float remain_freq = (float) _g_edges_reads[node][i] - offsets[i];
                    if (remain_freq > 0)
                    {
                        freqs.push_back(remain_freq);
                        neighs.push_back(neighbors[i]);
                        __traslation_map[num_neigh++] = i;
                    }
                }
                // The best way to express this is to define were the 1's are by index (n_vars^n)
                size_t n = in_flows_local.size(), n_vars = neighs.size();
                float num_its = pow(n_vars,n);
                float min_diff = INF, cur_diff = 0, control = INF;
                if (node == control)
                    cout <<"Node: "<<node <<" Sizes: "<<n<<" "<<n_vars<<" "<<num_its<<endl;
                vector<float> best_comb = vector<float>(n,-1), cur_comb = vector<float>(n,0);
                //cout << "Working over possibilities: "<<num_its<<endl;
                for (size_t i = 0; i < num_its; ++i)
                {
                    cur_comb[n-1] = i % n_vars;
                    
                    if (node == control){
                    cout << "Cur comb:"<<endl;
                    for (auto c:cur_comb)
                        cout << " "<<c;
                    cout << endl;}
                    vector<float> diff_per_freq = freqs;
                    for (size_t j = 0; j < n; ++j)
                        diff_per_freq[cur_comb[j]] -= in_flows_local[j];
                    // Interval confidence
                    for (size_t j = 0; j < n_vars; ++j)
                        diff_per_freq[j] *= ((diff_per_freq[j] < MIN_FLOW_PATH)&&(diff_per_freq[j] > (-MIN_FLOW_PATH)))?0:1;
                    cur_diff = 0;
                    for (auto freq: diff_per_freq){
                        if (node == control)
                            cout  << "Differences: "<<freq<<", ";
                        cur_diff += abs(freq);
                    }
                    // Quadratic difference
                    cur_diff = pow((cur_diff),2);
                    if (node == control)
                        cout <<endl<<"Cur diff: "<<cur_diff<<endl;
                    if (cur_diff < min_diff)
                    {
                        min_diff = cur_diff;
                        best_comb = cur_comb;
                    }
                    int j = (n-1), carry_over = 1;
                    while (carry_over & j >= 0)
                    {
                        if (cur_comb[j] + 1 == n_vars)
                        {
                            carry_over = 1;
                            cur_comb[j] = 0;
                            j -= 1;
                        } else {
                            carry_over = 0;
                            cur_comb[j] += 1;
                        }
                    }
                }
                for (size_t i = 0; i < best_comb.size(); ++i){
                    best_comb[i] = __traslation_map[best_comb[i]];
                }
                if (node == control){
                    cout << "Best comb:"<<endl;
                    for (auto c:best_comb)
                        cout << " "<<c;
                    cout << endl;
                    cout << "Min diff: "<<min_diff<<endl;
                    cout << "Freqs:"<<endl;
                    for (auto f: freqs)
                        cout << "Freq: "<<f<<endl;
                    cout << "Inflows: "<<endl;
                    for (auto i_f:in_flows_local)
                        cout << "Inflow: "<<i_f<<endl;
                    cout << "Strain freq: "<<f0<<endl;
                    //cin.get();
                }
                // Remove last assingment
                if (f0 != 0)
                    best_comb.pop_back();
                return best_comb;
            };
            /*
             * Pre Distributor - check which neighs is pointed but which parent.
             */
            auto __check_pe = [this](OwnNode_t node, vector<OwnNode_t> in_parents, vector<OwnNode_t> neighbors)
            {
                bool show = (false), distributed = false;
                vector<vector<OwnNode_t>> distributed_parents = vector<vector<OwnNode_t>>(neighbors.size(), vector<OwnNode_t>());
                std::pair<bool, vector<vector<OwnNode_t>>> empty_sol = {false,vector<vector<OwnNode_t>>(neighbors.size(), vector<OwnNode_t>())};
                for (size_t i = 0; i < neighbors.size(); ++i){
                    if (show)
                        cout << "Neighbor: "<<neighbors[i] <<endl;
                    vector<OwnNode_t> nodes_to_check = vector<OwnNode_t>(1, neighbors[i]);
                    for (auto n: getNeighbor(neighbors[i]))
                        if (in_degree(n) == 1)
                            nodes_to_check.push_back(n);
                    vector<OwnNode_t> parent_assigned;
                    for (size_t j = 0; j < in_parents.size(); ++j)
                    {
                        if (show)
                            cout << "Parent: "<<in_parents[j] << endl;
                        auto in_parent = in_parents[j];
                        auto p_e = _g_nodes[in_parent]._paired_info;
                        if (in(p_e, nodes_to_check)){
                            parent_assigned.push_back(in_parent);
                            //distributed_parents[i].push_back(j);
                        }
                    }
                    /*if (!only_one)
                        return empty_sol;*/
                    if (parent_assigned.size() > 0) {
                        distributed = true;
                        distributed_parents[i] = parent_assigned;
                    }
                }
                empty_sol = {distributed, distributed_parents};
                return empty_sol;
            };
            /*
             * DAG builder
             */
            std::function<void(size_t, size_t,size_t&, size_t&,
                               vector<size_t> &,vector<size_t>&, OwnNode_t , vector<float>, float)> __reach = [this, &__distribute,&__reach,
                    &traslation_map,
                    &local_graph_2,
                    &union_pei,&show,
                    &strain_freq,
                    &cautious_execution,&val_node,&val_neigh,
                    &basic_information,
                    &true_nodes,
                    &minimum_reached_by,&reached_by,&reverse_map,&frequencies_map](size_t  src, size_t target,
                                 size_t steps, size_t &branches
                    , vector<size_t> & checked
                    ,vector<size_t> & median_flow,
                                 OwnNode_t p_parent, vector<float> accumulated_flow, float freq_strain_updated)->void{
                if (Parameters::get().debug) {
                    basic_information << "( Steps: " << steps << "," << ")" << target << "(" << src << ", In_cum_flow_length: " << accumulated_flow.size()
                                      << ", In_cum_flow: " << Common::sum_vector(accumulated_flow) << ")";
                    for (auto f:accumulated_flow)
                        basic_information << " " <<f<<" ";
                }
                bool swap = false;
                size_t new_steps = 0, new_branches = 0;
                if (src != target && (in(union_pei, target))) {
                    checked[traslation_map[target]]++;
                    /*
                     * Testing flow
                     */
                    float testing_flow = abs((median_flow[median_flow.size()-1]-Common::sum_vector(accumulated_flow)));
                    // Avoid negative flows
                    testing_flow = (testing_flow < 0)?0:testing_flow;
                    if (Parameters::get().debug)
                        basic_information << " New flow expected: "<<testing_flow<<" ";
                    float flow_expected = 0.0;
                    if (STRATEGY == 0)
                    {
                        sort (median_flow.begin(), median_flow.end());
                        flow_expected = Maths::my_round(Maths::median(median_flow), true);
                    } else if (STRATEGY == 1)
                    {
                        sort (median_flow.begin(), median_flow.end());
                        flow_expected = Maths::my_round(median_flow[0], true);
                    } else if (STRATEGY == 2 || STRATEGY == 3)
                    {
                        //flow_expected = Maths::my_round(testing_flow, true);
                        flow_expected = round(testing_flow);
                    }
                    freq_strain_updated = flow_expected;
                    if (Parameters::get().debug) {
                        for (auto freq: median_flow)
                            basic_information << " " << freq << " ";
                        basic_information << " Flow expected: "<<flow_expected<<" Optimal strain freq: "<<strain_freq;
                    }
                    local_graph_2.addEdge(traslation_map[src], traslation_map[target], flow_expected, strain_freq);
                    src = target;
                    swap = true;
                }
                /*
                 * Flow from extra parents - not before
                 */
                vector<OwnNode_t> parents = getNeighbor(target, false, p_parent);
                if (parents.size() > 0)
                {
                    for (auto parent:parents)
                        accumulated_flow.push_back((float)_get_edge_frec(parent, target));
                }
                if (Parameters::get().debug)
                    basic_information << "(pre_parent - "<<p_parent<<", "<<abs(Common::sum_vector(accumulated_flow))<<")-> ";
                if (steps > D_MAX_PATH && !swap) {
                    if (Parameters::get().debug)
                        basic_information << "Maximum path allowed "<< steps<<endl;
                    return;
                }
                if (branches == D_MAX_BRANCHES && !swap) {
                    if (Parameters::get().debug)
                        basic_information << "Maximum branches allowed "<<branches<<endl;
                    return;
                }
                auto neighbors = getNeighbor(target);
                if (cautious_execution && (target == val_node || target == val_neigh)) {
                    steps += Parameters::get().kmerSize - 1;
                    new_steps += Parameters::get().kmerSize - 1;
                } else {
                    steps += getLength(target) - Parameters::get().kmerSize + 1;
                    new_steps += getLength(target) - Parameters::get().kmerSize + 1;
                }
                branches += (neighbors.size() > 1)?1:0;
                float curr_strain_freq = 0.0;
                if (STRATEGY == 2) {
                    curr_strain_freq = (float) _get_edge_frec(p_parent, target);
                    sort(accumulated_flow.begin(), accumulated_flow.end());
                    if (Parameters::get().debug) {
                        basic_information << "Current strain freq: " << curr_strain_freq<<" Previous parent: "<<p_parent<<" Target: "<<target << " "<<endl;
                        for (auto f:accumulated_flow)
                            basic_information << " "<<f<<" ";
                        basic_information << endl;
                    }
                }
                /*
                 * Get best comb
                 */
                vector<vector<float>> distributed_flow = vector<vector<float>>(neighbors.size(), vector<float>());
                if (STRATEGY == 3){
                    if (accumulated_flow.size() > 0 && neighbors.size() > 1)
                    {
                        basic_information<<endl << "Getting best combination: "<<src<<" "<<target<<" Updated Frequency: "<<freq_strain_updated<<endl;
                        vector<float> neigh_offsets;
                        auto best_comb = __distribute((float)strain_freq, target, neighbors, accumulated_flow, neigh_offsets);
                        vector<float> neigh_freqs;
                        for (auto f:_g_edges_reads[target])
                            neigh_freqs.push_back(f);
                        for (size_t i = 0; i < best_comb.size(); ++i)
                        {
                            basic_information << best_comb[i] << " In_flow "<<accumulated_flow[i]
                                <<" Best neigh: "<<neighbors[best_comb[i]]<<" with freq: "<<(float)(_g_edges_reads[target][best_comb[i]]) << ", ";
                            neigh_freqs[best_comb[i]] = ((neigh_freqs[best_comb[i]] - accumulated_flow[i] < 0)?0:neigh_freqs[best_comb[i]] - accumulated_flow[i]);
                            distributed_flow[best_comb[i]].push_back(accumulated_flow[i]);
                        }
                        basic_information << endl;
                        for (size_t f = 0; f < neigh_freqs.size(); ++f){
                            basic_information << neighbors[f] << " Freq Remaining: "<<neigh_freqs[f]<<" Flow: ";
                            for (auto flow:distributed_flow[f])
                                basic_information << " "<<flow<<", ";
                        }
                        basic_information << endl;
                    } else if ( neighbors.size() > 0)
                        distributed_flow[0] = accumulated_flow;
                }
                for (size_t i = 0; i < neighbors.size(); ++i) {
                    OwnNode_t neigh = neighbors[i];
                    /*if (checked[traslation_map[neigh]] < 1 || !in(union_pei, neigh))
                    {*/
                        /*
                         * Readapt accumulated flow from new neighs
                         */
                        /*if (Parameters::get().debug)
                            basic_information << "Number of in_flows: "<<accumulated_flow.size()<<" ";*/
                        vector<float> reaccumulated_flow;
                        if (STRATEGY == 2) {
                            reaccumulated_flow = accumulated_flow;
                            if (Parameters::get().debug)
                                basic_information << endl << "Neighbor: "<<neigh<<endl;
                            for (size_t j = 0; j < neighbors.size(); ++j) {
                                if (i == j)
                                    continue;
                                float remain_neigh_freq = (float) _g_edges_reads[target][j];
                                if (Parameters::get().debug)
                                    basic_information << endl << "      Removing flow for: "<<neighbors[j]<<" ";
                                /*if (Parameters::get().debug)
                                    basic_information << "Freq searched: "<<neigh_freq<<" ";*/
                                size_t prev_index = 0;
                                float bound = max((float) _g_edges_reads[target][j] * (float) FLOW_LIMIT_RATIO,(float) MAXIMUM_PATH_DEVIATION);
                                float ic_left = -bound,ic_right = bound;
                                if (Parameters::get().debug)
                                    basic_information << "Confidence interval: " << ic_left << " " <<
                                        ic_right<<" Remain freq: "<<remain_neigh_freq<<" "<<endl;
                                bool convergence = false;
                                while (!convergence) {
                                    convergence = true;
                                    for (int f = accumulated_flow.size() - 1; f >= 0; --f) {
                                        if (reaccumulated_flow[f] == 0)
                                            continue;
                                        float tmp_rem_flow = remain_neigh_freq - reaccumulated_flow[f];
                                        if ((tmp_rem_flow <= ic_right) & (tmp_rem_flow >= ic_left)) {
                                            if (Parameters::get().debug)
                                                basic_information << "Flow selected: " << reaccumulated_flow[f]
                                                                  << " End flow Explanation" << endl;
                                            reaccumulated_flow[f] = 0;
                                            remain_neigh_freq = 0;
                                            break;
                                        } else if (tmp_rem_flow >= ic_right) {
                                            if (Parameters::get().debug)
                                                basic_information << "Selected: " << reaccumulated_flow[f] << " " << endl;
                                            remain_neigh_freq -= reaccumulated_flow[f];
                                            reaccumulated_flow[f] = 0;
                                            convergence = false;
                                        }
                                    }
                                }
                                if (Parameters::get().debug)
                                    basic_information << "Unexplained flow: "<<remain_neigh_freq<<" Freq strain? "<<curr_strain_freq<<endl;
                                float tmp_rem_flow = remain_neigh_freq - curr_strain_freq;
                                if ((tmp_rem_flow < ic_right) & (tmp_rem_flow > ic_left))
                                    remain_neigh_freq = 0;
                                if (remain_neigh_freq > 0) {
                                    for (int f = accumulated_flow.size() - 1; f >= 0; --f)
                                    {
                                        if (reaccumulated_flow[f] > 0 && (reaccumulated_flow[f] - remain_neigh_freq) < 0) {
                                            cout << "This cannot be!" << endl;
                                            basic_information << "Check! "<<accumulated_flow[f]<<" "<<reaccumulated_flow[f]<<" "<<remain_neigh_freq<<endl;
                                        }
                                        if (reaccumulated_flow[f] > 0) {
                                            reaccumulated_flow[f] -= remain_neigh_freq;
                                            remain_neigh_freq = 0;
                                        }
                                    }
                                }
                                if (Parameters::get().debug)
                                    basic_information << "Final unexplained flow: "<<remain_neigh_freq<<endl;
                            }
                            size_t removed = 0;
                            for (size_t j = 0; j < accumulated_flow.size(); ++j) {
                                if (reaccumulated_flow[j] == 0) {
                                    reaccumulated_flow.erase(reaccumulated_flow.begin() + j - removed++);
                                    /*if (Parameters::get().debug)
                                        basic_information << "Flow removed: "<<accumulated_flow[j]<<" ";*/
                                }
                            }
                        }
                        if (STRATEGY == 3)
                            reaccumulated_flow = distributed_flow[i];
                        if (!swap) {
                            /*
                             * Once you swap the current strain you have to reset the fake_nodes
                             */
                            vector<size_t> new_median_flow = median_flow;
                            new_median_flow.push_back(_g_edges_reads[target][i]);
                            size_t local_steps = steps;
                            __reach(src, neigh, local_steps, branches, checked, new_median_flow, target, reaccumulated_flow, freq_strain_updated);
                        } else {
                            if (Parameters::get().debug)
                                basic_information << " New src: "<<neigh<< endl;
                            vector<size_t> new_median_flow(1, _g_edges_reads[target][i]);
                            __reach(src, neigh, new_steps, new_branches, checked, new_median_flow, target, reaccumulated_flow, freq_strain_updated);
                        }
                }
                if (Parameters::get().debug)
                    basic_information << "No more neighs, OUT!"<<target<<" "<<src<<endl;
            };
            /*
             * DAG builder 2.0
             */
            std::function<void(OwnNode_t,OwnNode_t,size_t,size_t,unordered_map<OwnNode_t, unordered_set<OwnNode_t>>&,
                            float,OwnNode_t,vector<float>,float,OwnNode_t,vector<OwnNode_t>, float)> 
                __dag_builder = [this, &__distribute, &__check_pe,&__dag_builder, 
            &local_graph_2, &show, &union_pei, 
            &basic_information, &true_nodes, &traslation_map](OwnNode_t src, OwnNode_t target, 
                        size_t steps, size_t branches, unordered_map<OwnNode_t, unordered_set<OwnNode_t>> & p_parents_map,
                        float last_cov_edge,OwnNode_t p_parent, vector<float> accumulated_flow, float freq_strain_updated, OwnNode_t strain_parent,
                        vector<OwnNode_t> in_parents, float original_freq)->void{
                float in_flows = Common::sum_vector(accumulated_flow);
                float new_freq = last_cov_edge;
                
                basic_information << "Source Node: "<<src<<" Target: "<<target
                    <<" Steps: "<<steps<<" Previous Parent: "<<p_parent<<" InFlows: "<<in_flows
                    <<" Strain Freq: "<<freq_strain_updated<<" Freq last Edge: "<<new_freq<<" Parent Strain: "<<strain_parent<<endl;
                basic_information << "Adding the node: "<<in(union_pei, target)<<endl;
                if (p_parents_map.find(target) == p_parents_map.end())
                    p_parents_map[target] = {p_parent};
                else
                    p_parents_map[target].emplace(p_parent);
                if ((src != target) && in(union_pei, target))
                {
                    if (STRATEGY == 3)
                    {
                        new_freq = std::min(freq_strain_updated, original_freq);//(new_freq - in_flows < 0)?1:(new_freq-inflows);
                    }
                    local_graph_2.addEdge(traslation_map[src], traslation_map[target], new_freq, last_cov_edge);
                    basic_information << "  New DAG edge from "<<src<<" to "<<target<<" with freq "<<new_freq<<endl;
                    // If I am reached as much times as in_degree I have I skip myself.    
                    /*if (checked[traslation_map[target]] >= 1)
                    {
                        basic_information << "Maximum number of traversals reached: "<<src<<" stopping traversion."<<endl;
                        return;
                    }
                    checked[traslation_map[target]] += 1;*/
                    steps = 0;
                    src = target;
                    branches = 0;
                }
                if (steps > D_MAX_PATH || branches > D_MAX_BRANCHES)
                {
                    basic_information << "Boundary reached!"<<endl;
                    return;
                }
                /*
                * Flow from extra parents - not before
                */
                if (STRATEGY == 3){
                    vector<OwnNode_t> parents = getNeighbor(target, false, p_parent);
                    if (parents.size() > 0)
                    {
                        for (auto parent:parents){
                            // Polish parents based on previous traversed parents for target
                            if (p_parents_map[target].find(parent) == p_parents_map[target].end()){
                                accumulated_flow.push_back((float)_get_edge_frec(parent, target));
                                in_parents.push_back(parent);
                            }
                        }
                    }
                }
                // Update boundary limits
                steps += (_g_nodes[target]._length - Parameters::get().kmerSize + 1);
                auto neighbors = getNeighbor(target);
                branches += (neighbors.size() > 1)?1:0;
                /*
                * Get best comb
                */
                vector<vector<float>> distributed_flow = vector<vector<float>>(neighbors.size(), vector<float>());
                vector<vector<OwnNode_t>> distributed_parents = vector<vector<OwnNode_t>>(neighbors.size(), vector<OwnNode_t>());
                // Store parent_suggestion
                OwnNode_t neigh_assigned = INF; 
                if (STRATEGY == 3)
                {
                    // Get parent suggestion  
                    if (strain_parent != target){
                        basic_information << "Get parent suggestion from: "<<strain_parent<<endl;
                        auto p_e = _g_nodes[strain_parent]._paired_info;
                        for (auto local_neigh:neighbors){
                            vector<OwnNode_t> neighs_to_check(1, local_neigh);
                            for (auto n: getNeighbor(local_neigh))
                                if (in_degree(n) == 1)
                                    neighs_to_check.push_back(n);
                            if (in(p_e, neighs_to_check))
                            {
                                if (neigh_assigned == INF)
                                    neigh_assigned = local_neigh;
                                else
                                {
                                    neigh_assigned = INF;
                                    break;
                                }
                            }
                        }
                    }
                    if (neigh_assigned != INF)
                        basic_information<<"We can use parent information, neigh suggested "<<neigh_assigned<<endl;
                    float unassigned_flow = freq_strain_updated;
                    if (accumulated_flow.size() > 0 && neighbors.size() > 1)
                    {
                        vector<float> neigh_freqs, neigh_offsets;
                        for (size_t f = 0; f < _g_edges_reads[target].size(); ++f){
                            neigh_freqs.push_back(_g_edges_reads[target][f]);
                            neigh_offsets.push_back(0);
                            if (neighbors[f] == neigh_assigned){
                                neigh_offsets[f] = min(freq_strain_updated, neigh_freqs[f]);
                                neigh_freqs[f] -= neigh_offsets[f];
                                unassigned_flow -= neigh_offsets[f];
                            }
                        }
                        basic_information << endl << "First step using paired-end information"<<endl;
                        // First try to distributed parents by paired-end information
                        auto _distributed_parents = __check_pe(p_parent,in_parents, neighbors);
                        if (_distributed_parents.first)
                        {
                            basic_information <<"There are associations"<<endl;
                            //Getting parents distribution - univoque
                            unordered_map<OwnNode_t, vector<OwnNode_t>> parent_neighs_map;
                            for (size_t w = 0; w < _distributed_parents.second.size(); ++w)
                            {
                                // Totally greedy approach
                                /*if (_distributed_parents.second[w].size() == 1)
                                {

                                } else*/ 
                                    for (auto index: _distributed_parents.second[w])
                                    {
                                        basic_information << "Parent: "<<index<<" points to: "<<neighbors[w]<<", ";
                                        if (parent_neighs_map.find(index) == parent_neighs_map.end()){
                                            parent_neighs_map[index] = vector<OwnNode_t>(1, neighbors[w]);
                                        } else {
                                            parent_neighs_map[index].push_back(neighbors[w]);
                                        }
                                        // If there is a match we give the minimum possible flow.
                                        auto iterator_index = find(in_parents.begin(), in_parents.end(), index) - in_parents.begin();
                                        accumulated_flow[iterator_index] -= MIN_FLOW_PATH;
                                        distributed_flow[w].push_back(MIN_FLOW_PATH);
                                        distributed_parents[w].push_back(index);
                                        neigh_offsets[w] += MIN_FLOW_PATH;
                                        neigh_freqs[w] -= MIN_FLOW_PATH;
                                    }
                            }
                            basic_information << "End showing associations"<<endl;
                            vector<bool> removed = vector<bool>(accumulated_flow.size(), false);
                            for (size_t w = 0; w < _distributed_parents.second.size(); ++w)
                            {
                                for (auto index: _distributed_parents.second[w])
                                {
                                    // Only 1 to 1 associations can be pre-known
                                    if (parent_neighs_map[index].size() == 1){
                                        auto iterator_index = find(in_parents.begin(), in_parents.end(), index) - in_parents.begin();
                                        float flow_from_parent = (accumulated_flow[iterator_index] - neigh_freqs[w] < 0)?accumulated_flow[iterator_index]:neigh_freqs[w];
                                        auto place = find(distributed_parents[w].begin(), distributed_parents[w].end(), index) - distributed_parents[w].begin();
                                        if (place == distributed_parents[w].size()){
                                            distributed_flow[w].push_back(flow_from_parent);
                                            distributed_parents[w].push_back(index);
                                        } else {
                                            distributed_flow[w][place]+= flow_from_parent;
                                            basic_information <<" Flow "<<distributed_flow[w][place]<<" ";
                                        }
                                        neigh_freqs[w] = neigh_freqs[w] - flow_from_parent;
                                        neigh_offsets[w] += flow_from_parent;
                                        // Change flow in parent
                                        accumulated_flow[iterator_index] -= flow_from_parent;
                                        accumulated_flow[iterator_index] *= (accumulated_flow[iterator_index] < MIN_FLOW_PATH)?0:1;
                                        basic_information << "One to one association between "<<index<<" "<<neighbors[w]<<" Offset: "<<neigh_offsets[w]
                                            <<" Remaining flow in parent: "<<accumulated_flow[iterator_index]<<" ";
                                        // Remove already explained flows
                                        if (accumulated_flow[iterator_index] < MIN_FLOW_PATH)
                                            removed[iterator_index] = true;
                                        if (neigh_freqs[w] < MIN_FLOW_PATH){
                                            neigh_offsets[w] += neigh_freqs[w];
                                            neigh_freqs[w] = 0;
                                        }
                                        basic_information << "Parent removal: "<<iterator_index<<endl;
                                    }
                                }
                            }
                            size_t offset = 0;
                            for (size_t r_i = 0; r_i < removed.size(); r_i++)
                            {
                                if (removed[r_i]){
                                    accumulated_flow.erase(accumulated_flow.begin() + r_i - offset);
                                    in_parents.erase(in_parents.begin() + r_i - offset++);
                                }
                            }
                        }
                        // Simplified distribution
                        basic_information<<endl << "Getting best combination per parent unnasigned: "<<src<<" "<<target<<" Frequency to distribute: "<<unassigned_flow<<endl;
                        basic_information<<" Accumulated flow still not distributed: "<<accumulated_flow.size()<<endl;
                        if (accumulated_flow.size() > 0){
                            size_t number_of_iteration = 0;
                            bool remain_flow = true;
                            while (remain_flow){
                                auto best_comb = __distribute(unassigned_flow, target, neighbors, accumulated_flow, neigh_offsets);
                                if (best_comb[0] == INF)
                                    break;
                                for (auto c:best_comb)
                                    basic_information << c<<", ";
                                basic_information << endl;
                                remain_flow = false;
                                for (size_t i = 0; i < best_comb.size(); ++i)
                                {
                                    basic_information << best_comb[i] << " In_flow "<<accumulated_flow[i]<<" In parent: "<<in_parents[i]
                                        <<" Best neigh: "<<neighbors[best_comb[i]]<<" with freq: "<<(float)(_g_edges_reads[target][best_comb[i]]) 
                                        <<" Sanity: "<<neigh_freqs[best_comb[i]] <<" offsets: "<<neigh_offsets[best_comb[i]]
                                        << " and already assigned "<< Common::sum_vector(distributed_flow[best_comb[i]])<<", ";
                                    // Commuatatitivity problems (it depends on the order)
                                    float flow_inserted = 0;
                                    if ((neigh_freqs[best_comb[i]] - accumulated_flow[i]) < -MIN_FLOW_PATH)
                                    {
                                        flow_inserted = neigh_freqs[best_comb[i]];
                                        accumulated_flow[i] -= neigh_freqs[best_comb[i]];
                                        neigh_offsets[best_comb[i]] += neigh_freqs[best_comb[i]];
                                        neigh_freqs[best_comb[i]] = 0;
                                    } else if ((neigh_freqs[best_comb[i]] - accumulated_flow[i]) < 0){
                                        neigh_offsets[best_comb[i]] += neigh_freqs[best_comb[i]];
                                        flow_inserted = neigh_freqs[best_comb[i]];
                                        neigh_freqs[best_comb[i]] = 0;
                                        accumulated_flow[i] = 0;
                                    } else {
                                        neigh_freqs[best_comb[i]] -= accumulated_flow[i];
                                        neigh_offsets[best_comb[i]] += accumulated_flow[i];
                                        neigh_offsets[best_comb[i]] += (neigh_freqs[best_comb[i]] < MIN_FLOW_PATH)?neigh_freqs[best_comb[i]]:0;
                                        neigh_freqs[best_comb[i]] *= (neigh_freqs[best_comb[i]] < MIN_FLOW_PATH)?0:1;
                                        flow_inserted = accumulated_flow[i];
                                        accumulated_flow[i] = 0;
                                    }
                                    if (flow_inserted != 0)
                                    {
                                        OwnNode_t w = best_comb[i];
                                        auto place = find(distributed_parents[w].begin(), distributed_parents[w].end(), in_parents[i]) - distributed_parents[w].begin();
                                        if (place == distributed_parents[w].size()){
                                            distributed_flow[w].push_back(flow_inserted);
                                            distributed_parents[w].push_back(in_parents[i]);
                                        } else {
                                            distributed_flow[w][place] += flow_inserted;
                                        }
                                        /*distributed_flow[best_comb[i]].push_back(flow_inserted);
                                        distributed_parents[best_comb[i]].push_back(in_parents[i]);*/
                                    }
                                    if (accumulated_flow[i] != 0)
                                        remain_flow = true;
                                }
                                bool all_empty = true;
                                for (size_t i = 0; i < neigh_freqs.size(); ++i)
                                {
                                    all_empty &= (neigh_freqs[i] == 0);
                                }
                                if (remain_flow)
                                    remain_flow = !(all_empty);
                                basic_information << endl;
                            }
                        }
                        for (size_t f = 0; f < neigh_freqs.size(); ++f){
                            basic_information << neighbors[f] << " Freq Remaining: "<<neigh_freqs[f]<<" Flow: ";
                            for (size_t j = 0; j < distributed_flow[f].size(); ++j)
                                basic_information << " "<<distributed_flow[f][j]<<" parent: "<<distributed_parents[f][j]<<", ";
                        }
                        basic_information << endl;
                    } else if ( neighbors.size() > 0) {
                        distributed_flow[0] = accumulated_flow;
                        distributed_parents[0] = in_parents;
                    }
                }
                for (size_t i = 0; i < neighbors.size(); ++i) {
                    basic_information << "Neighbor: "<<traslation_map[neighbors[i]]<<" "<<neighbors[i] <<endl;
                    OwnNode_t neigh = neighbors[i];
                    vector<float> reaccumulated_flow;
                    vector<OwnNode_t> strain_parents;
                    if (STRATEGY == 3){
                        reaccumulated_flow = distributed_flow[i];
                        strain_parents = distributed_parents[i];
                    }
                    float new_strain_freq = (float)_g_edges_reads[target][i] - Common::sum_vector(reaccumulated_flow);
                    /*if (neigh_assigned == neighbors[i])
                        new_strain_freq = (new_strain_freq < freq_strain_updated)?freq_strain_updated:new_strain_freq;
                    else*/
                        new_strain_freq = (new_strain_freq < MIN_FLOW_PATH)?0:new_strain_freq;
                    // Greedy as hell
                    if (new_strain_freq == 0){
                        basic_information << "No enough flow for this strain"<<endl;
                        continue;
                    }
                    if (neighbors.size() > 1)
                        strain_parent = neighbors[i];
                    __dag_builder(src, neigh, steps, branches, p_parents_map,(float)_g_edges_reads[target][i], target, 
                        reaccumulated_flow,new_strain_freq,strain_parent,strain_parents, original_freq);
                }
                basic_information << "No more neighs, OUT!"<<target<<" "<<src<<endl;
            };
            /*
             * Finding paths
             */
            if (Parameters::get().debug)
                basic_information << "Building the network"<<endl;
            for (auto pn:union_pei)
            {
                size_t steps = 0, branches = 0;
                vector<size_t> checked = vector<size_t>(union_pei.size(), 0);
                vector<size_t> median_flow, extra_in_parents;
                vector<float> extra_flow_accumulated;
                if (!cautious_execution) {
                    if (reached_by[traslation_map[pn]] == minimum_reached_by) {
                        checked[traslation_map[pn]]++;
                        if (Parameters::get().debug)
                            basic_information << "NEW LAUNCH" << endl;
                        __reach(pn, pn, steps, branches, checked, median_flow, pn, extra_flow_accumulated, strain_freq);
                    }
                } else if (parent_in_union )
                {
                    if (pn == val_node) {
                        median_flow.push_back(strain_freq);
                        cout << "Launch parent: "<<pn<<endl;
                        unordered_map<OwnNode_t, unordered_set<OwnNode_t>> p_parents_map;
                        if (Parameters::get().debug)
                            basic_information << "Launching source (parent): " << pn << endl;
                        OwnNode_t parent_strain = pn;
                        if (out_degree(pn) > 1)
                            parent_strain = val_neigh;
                        //__reach(pn, val_neigh, steps, branches, checked, median_flow, pn, extra_flow_accumulated, strain_freq);
                        __dag_builder(pn, val_neigh, steps, branches,p_parents_map, (float)strain_freq, pn, 
                            extra_flow_accumulated, strain_freq, parent_strain,extra_in_parents, (float)strain_freq);
                    }
                } else if (son_in_union) {
                    if (pn == val_neigh)
                    {
                        if (Parameters::get().debug)
                            basic_information << "Launching source (neigh): " << pn << endl;
                        __reach(pn, pn, steps, branches, checked, median_flow, pn, extra_flow_accumulated, strain_freq);
                    }
                }
            }
            /*
             * Correct flow in local graphs
             */
            string adjust_file = "debug/"+to_string(_g_nodes[i]._val)+"_"+to_string(_g_nodes[n]._val)+"_adjust_flow.txt";
            std::ofstream adjust_information;
            if (CORRECT_GRAPH)
            {
                if (Parameters::get().debug)
                {
                    string graph_file = "debug/"+to_string(_g_nodes[i]._val)+"_"+to_string(_g_nodes[n]._val)+"_graph_prereadjust_1.txt";
                    local_graph_2.export_to_graph(sequence_map, graph_file);
                    basic_information << "Readjusting nodes/edges capacities"<<endl;
                }
                if (Parameters::get().debug)
                {
                    adjust_information.open(adjust_file, std::ofstream::out);
                }
                local_graph_2.correct_graph_safe(traslation_map[node_id],adjust_information, rescoring, *(this));
            }
            /*
             * Edge and strain freq flow readjustment
             */
            if (READJUST)
            {
                if (Parameters::get().debug)
                {
                    string graph_file = "debug/"+to_string(_g_nodes[i]._val)+"_"+to_string(_g_nodes[n]._val)+"_graph_prereadjust.txt";
                    local_graph_2.export_to_graph(sequence_map, graph_file);
                    basic_information << "Readjusting nodes/edges capacities"<<endl;
                }
                local_graph_2.readjust_flow(strain_freq,min_flow, max_flow);
                // Re-correct flow
                local_graph_2.correct_graph_safe(traslation_map[node_id],adjust_information, rescoring,*(this));
                strain_freq_check = ceil((MAX_GRANULARITY_ALLOWED*100)*abs((float) max_flow)/max_flow + MAX_GRANULARITY_ALLOWED);
                strain_freq_check = strain_freq;
            }
            if (Parameters::get().debug)
                basic_information << "Maximum expected flow: "<<max_flow
                        <<" Minimum expected flow: "<<min_flow<<" Expected flow: "<<strain_freq_check<<endl;
            /*
             * Remove extra edges (cycles)
             */
            local_graph_2.post_process_cycles(reached_by, minimum_reached_by);
            /*
             * Correct the minimum
             */
            for (auto pn: union_pei)
                reached_by[traslation_map[pn]] -= minimum_reached_by;
            //local_graph_2.print(sequence_map);
            if (Parameters::get().debug)
            {
                string gfa_file = "debug/"+to_string(_g_nodes[i]._val)+"_"+to_string(_g_nodes[n]._val)+"_graph.gfa"
                , graph_file = "debug/"+to_string(_g_nodes[i]._val)+"_"+to_string(_g_nodes[n]._val)+"_graph.txt";
                local_graph_2.export_to_gfa(sequence_map,gfa_file);
                local_graph_2.export_to_graph(sequence_map, graph_file);
            }
            priority_queue<pair<size_t,vector<OwnNode_t>>>
                                       cliques;
            if (local_graph_2.vertices() > 0)
                cliques = local_graph_2.report_min_cost_flow
                                               (reached_by,strain_freq
                                                       , Parameters::get().debug, true_nodes,"debug/"+to_string(_g_nodes[i]._val)
                                                       +"_"+to_string(_g_nodes[n]._val)+"_min_cost_flow.txt", rescoring);
            /*priority_queue<pair<size_t,vector<OwnNode_t>>>
                    cliques = local_graph.report_maximal_cliques(reverse_map, *this, frequencies_map, show);*/
            vector<pair<size_t,vector<OwnNode_t>>> cliques_store;
            unordered_set<OwnNode_t> nodes_checked;
            float max_size = 0.0;
            if (Parameters::get().debug)
                paths << "Cliques size: "<<cliques.size()<<endl;
            vector<pair<size_t,pair<UG_Node,UG_Node>>> potential_nodes;
            while (!cliques.empty() && nodes_checked.size() < union_pei.size())
            {
                pair <size_t, vector<OwnNode_t>> top_click = cliques.top();
                cliques.pop();
                max_size = (max_size < top_click.second.size()) ? top_click.second.size() : max_size;
                cliques_store.push_back(top_click);
            }
            vector<bool> offset_erase;
            for (auto top_click_pair: cliques_store)
            {
                vector<OwnNode_t> top_click = top_click_pair.second;
                if (top_click.size() < SIZE_RELATION*max_size)
                    continue;
                if (show)
                    paths << "New Clique:" << endl;
                Pairedendinformation_t clique, p_info_a, p_info_b;
                bool in_a = false, in_b = false, isa = true, isb = true, parent_in = false, son_in = false;
                float rescore = 0.0;
                for (auto node: top_click)
                {
                    nodes_checked.emplace(node);
                    //UG::UG_Node real_node = local_graph[node];
                    rescore += rescoring[node];
                    UG::UG_Node real_node = local_graph_2[node];
                    clique.emplace(real_node._val);
                    if (Parameters::get().debug) {
                        paths << node<<":"<<real_node._val<<":"<<_g_nodes[real_node._val]._abundance
                            <<":"<<real_node._id<<" "<<Common::return_unitig(sequence_map, real_node._val)<<endl;
                    }
                    if (real_node._val == val_node)
                        parent_in = true;
                    if (real_node._val == val_neigh)
                        son_in = true;
                    if (node_pei.find(real_node._val) != node_pei.end()) {
                        in_a = true;
                        p_info_a.emplace(real_node._val);
                    } else {
                        isa = false;
                    }
                    if (neigh_pei.find(real_node._val) != neigh_pei.end()) {
                        in_b = true;
                        p_info_b.emplace(real_node._val);
                    } else
                        isb = false;
                }
                /*
                 * Adjust ad-hoc
                 */
                float sum_rescoring = Common::sum_vector(rescoring);
                float denom = sum_rescoring-1;
                denom = (denom == 0)?1:denom;
                rescore = 1 + ((top_click.size() - 2 > 0)?(rescore - 1)/ denom:0);
                //rescore = (top_click.size() - 2 > 0)?(rescore - 1)/ (top_click.size() - 1):1;
                if (Parameters::get().debug) {
                    paths << endl;
                    paths << "ISA: "<<isa<<" ISB: "<<isb<<endl;
                    paths << "INA: "<<in_a<<" INB: "<<in_b<<endl;
                    paths << "Is Parent: "<<parent_in<<" Is son: "<<son_in<<endl;
                }
                if (parent_in && son_in) {
                    if (Parameters::get().debug) {
                        paths << "Path forced to be added: Contains parent and son. Kept? "<<(in_a & in_b)
                            << " With flow: "<<(((float)top_click_pair.first)/max((float)0.01,rescore))<<" Rescore: "<<rescore
                            <<" Priori: "<<top_click_pair.first<<" strain freq "<<strain_freq<<" Size: "<<top_click.size()<< endl;
                    }
                    offset_erase.push_back(!(in_a & in_b));
                } else if (cautious_execution)
                {
                    if (Parameters::get().debug) {
                        paths << "Path forced to be removed: Parent or son in the union set but not in the path, cannot be a real strain"
                                << endl;
                    }
                    offset_erase.push_back(true);
                } else {
                    if (Parameters::get().debug) {
                        paths << "Standard situation" << endl;
                        paths << "Both have to provide something to the path: (parent) "<<in_a<<" (son) "<<in_b<<endl;
                    }
                        /*if ((in_a && in_b) && ((isa || isb) & !(isa && isb)) && !content) {
                            offset_erase.push_back(true);
                        } else if (in_a && in_b)
                            offset_erase.push_back(false);
                        offset_node.push_back(isa);
                        offset_neigh.push_back(isb);*/
                    if (!in_a or !in_b)
                        offset_erase.push_back(true);
                    else
                        offset_erase.push_back(false);
                }
                UG_Node parent_l(_g_nodes[i]._val, p_info_a), son_l(_g_nodes[n]._val, p_info_b);
                potential_nodes.push_back({((RESCORING)?((float)top_click_pair.first)/max((float)0.01,rescore):top_click_pair.first)
                                           ,pair<UG_Node, UG_Node>(parent_l, son_l)});
                if (Parameters::get().debug) {
                    paths << "Nodes added as potential ones."<<endl;
                    paths << parent_l.to_String();
                    paths << son_l.to_String();
                }
            }
            /*
             * Procesar las paired-end y eliminar los subconjuntos
             */
            /*for (size_t index_s = 0; index_s < potential_nodes.size(); ++index_s)
            {
                UG_Node p = potential_nodes[index_s].second.first, h = potential_nodes[index_s].second.second;
                for (size_t index_s_2 = 0; index_s_2 < potential_nodes.size(); ++index_s_2)
                {
                    UG_Node p_2 = potential_nodes[index_s_2].second.first, h_2 = potential_nodes[index_s_2].second.second;
                    if (p_2 == p && h_2 == h)
                        continue;
                    bool parent_ss = is_subset(p._paired_info,p_2._paired_info), parent_ss_2 = is_subset(p_2._paired_info, p._paired_info)
                            ,son_ss = is_subset(h._paired_info, h_2._paired_info), son_ss_2 = is_subset(h_2._paired_info,h._paired_info);
                    if (parent_ss && !parent_ss_2 && (in_parents.size() > 1))
                        offset_erase[index_s] = true;
                    if (parent_ss_2 && !parent_ss && (in_parents.size() > 1))
                        offset_erase[index_s_2] = true;
                    if (son_ss && !son_ss_2 && (neighbors.size() > 1))
                        offset_erase[index_s] = true;
                    if (son_ss_2 && !son_ss && (neighbors.size() > 1))
                        offset_erase[index_s_2] = true;
                }
            }*/
            if (offset_erase.size() > 1) {
                /*size_t num_removes = 0;
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
                }*/
                size_t num_removes = 0;
                if (Parameters::get().debug)
                    paths << "Removing cliques: ";
                for (size_t pos = 0; pos < offset_erase.size(); pos++) {
                    if (offset_erase[pos])
                    {
                        if (Parameters::get().debug)
                            paths << " "<<pos;
                        potential_nodes.erase(potential_nodes.begin() + pos - num_removes++);
                    }
                }
                if (Parameters::get().debug)
                    paths << endl;
            }
            bool empty_solution = false;
            if (potential_nodes.size() == 0)
            {
                if (Parameters::get().debug)
                {
                    paths << "No potential nodes, added complete version of both (it only happens when none (parent and son) have paired-end information)"<<endl;
                }
                //Pairedendinformation_t empty;
                UG_Node parent_l(_g_nodes[i]._val, node_pei), son_l(_g_nodes[n]._val,neigh_pei);
                potential_nodes.push_back({strain_freq_check,pair<UG_Node,UG_Node>(parent_l,son_l)});
                empty_solution = false;
            }
            if (Parameters::get().debug)
                paths << endl<<"Potential nodes"<<"("<<potential_nodes.size()<<"): "<<endl;
            size_t nodes_added = 0;
            bool only_one = (potential_nodes.size() == 1);
            for (auto s_pair: potential_nodes) {
                pair<UG_Node, UG_Node> s = s_pair.second;
                if (Parameters::get().debug)
                {
                    paths << s.first.to_String();
                    paths << s.second.to_String();
                }
                if (((strain_freq_check - s_pair.first+ MAXIMUM_PATH_DEVIATION) < 0)
                    && (nodes_added > 0))
                    break;
                nodes_added++;
                pair<OwnNode_t, OwnNode_t> index_pairs = apdbg.addNode(s.first, s.second, ((only_one)?strain_freq:s_pair.first), show, empty_solution);
                if (Parameters::get().debug){
                    if (index_pairs.first != INF) {
                        paths << "Adding edge from: "<<index_pairs.first<<" "<<index_pairs.second<< " With flow: "<<((only_one)?strain_freq:s_pair.first)<<endl;
                        apdbg.add_to_split_map(s.first._val, index_pairs.first);
                        apdbg.add_to_split_map(s.second._val, index_pairs.second);
                    } else
                        paths << "Adding connections between previous node insertions and the complete neighbor"<<endl;
                }
                /*
                 * Expected flow reached
                 */
                strain_freq_check -= s_pair.first;
            }
            if (Parameters::get().debug) {
                paths << "Finally cliques added: "<<nodes_added<<" Out-of: " << potential_nodes.size()
                    <<" Remaining flow: "<<strain_freq_check<<" Out of: "<<strain_freq<<endl;
                cout << "Finally cliques added: "<<nodes_added<<" Out-of: " << potential_nodes.size()
                      <<" Remaining flow: "<<strain_freq_check<<" Out of: "<<strain_freq<<endl;
            }
        }
    }
        bool all_processed = true;
        for (size_t i = 0; i < checked.size(); ++i)
        {
            if (!checked[i]){
                checked[i] = true;
                nodes_process.push(i);
                all_processed = false;
                break;
            }
        }
        if (all_processed)
            break;
    }
    if (Parameters::get().debug)
    {
        rejected_nodes.close();
    }
    float total_checked = 0.0;
    for (auto i:checked)
        total_checked += i;
    cout << "% nodes checked: "<<(total_checked/(float) checked.size()*100)<<"%"<<endl;
    cout << "Number of empty pairs: "<<empty_pairs<<" out of "<<_g_nodes.size()<<endl;
    //exit(1);
    /*
     * Basic polishing
     */
    apdbg.polish(true);
}

size_t DBG::out_degree(OwnNode_t node_id)
{
    return getNeighbor(node_id).size();
}

size_t DBG::in_degree(OwnNode_t node_id)
{
    return getNeighbor(node_id, false).size();
    //return _g_in_edges[node_id].size();
}

vector<vector<OwnNode_t>> DBG::export_unitigs(const vector<string> & sequence_map, bool show)
{
    vector<vector<OwnNode_t>> unitigs;
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
                repetitive_call, unitigs, local_checked, show, sequence_map);
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

void DBG::post_process_pairs(const vector<string> & sequence_map)
{
    _get_internal_stats(sequence_map);
}

void DBG::_get_internal_stats(const vector<string> & sequence_map)
{
    vector<size_t> abundances_v, paired_freqs;
    for (auto v:_g_nodes)
    {
        if (v._active) {
            float n_times = (float) sequence_map[Common::return_index(sequence_map, v._val)].length() -
                            (float) Parameters::get().kmerSize + 1;
            for (float i = 0; i < n_times; ++i)
                abundances_v.push_back(v._abundance);
            paired_freqs.push_back(v._paired_info.size());
        }
    }
    /*
     * Report histogram like
     */
    if (Parameters::get().debug)
    {
        std::ofstream histogram_abundances;
        histogram_abundances.open("stats/abundance_histogram.txt", std::ofstream::out);
        for (auto v:abundances_v)
            histogram_abundances << v << endl;
        histogram_abundances.close();
        std::ofstream histogram_freq_pairs;
        histogram_freq_pairs.open("stats/histogram_freq_pairs.txt", std::ofstream::out);
        for (auto v:paired_freqs)
            histogram_freq_pairs << v << endl;
        histogram_freq_pairs.close();
    }
    sort (abundances_v.begin(), abundances_v.end());
    _max_ab = abundances_v[abundances_v.size()-1];_min_ab = abundances_v[0];
    _quantile_L = Maths::i_quantile(abundances_v, L_QUANTILE);
    /*
     * Remove nodes from left_side
     */
    cout << "Quantile_L: "<<_quantile_L<<endl;
    abundances_v.clear();
    for (auto v:_g_nodes)
    {
        if (v._active) {
            if (v._abundance >= _quantile_L) {
                abundances_v.push_back(v._abundance);
            }
        }
    }
    sort (abundances_v.begin(), abundances_v.end());
    _quantile_H = Maths::i_quantile(abundances_v, H_QUANTILE, false);
    cout << "Quantile_H: "<<_quantile_H<<endl;
    /*
     * Both mean and median are got after "removing" the required quantiles [5,85]
     */
    float mean = 0.0;
    vector<size_t> abundances_v_post_processed;
    for (auto v:_g_nodes)
    {
        if (v._active) {
            if (v._abundance <= _quantile_H && v._abundance >= _quantile_L) {
                abundances_v_post_processed.push_back(v._abundance);
                mean += v._abundance;
            }
        }
    }
    _mean_ab = mean / (float)abundances_v_post_processed.size();
    _median_ab = abundances_v_post_processed[abundances_v_post_processed.size()/2];
}
/*
 * Private
 */
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
/*
 * Adding the extra frequency
 */
void DBG::add_read(OwnNode_t u1, size_t quantity)
{
    if (_g_edges[u1].size() > 0) {
        size_t extra_freq = floor((float) quantity / (float) _g_edges[u1].size());
        for (size_t i = 0; i < _g_edges_reads[u1].size(); ++i)
            _g_edges_reads[u1][i] += extra_freq;
    }
}

void DBG::add_read(OwnNode_t u1, OwnNode_t u2, size_t quantity, bool direction)
{
    size_t pos = _find_edge(u1,u2);
    OwnNode_t node_increment = (direction)?u1:u2;
    //_g_nodes[node_increment].increase_abundance(quantity);
    if (pos != NO_NEIGH)
    {
        /*if (u1 == 0)
            cout << "Parent: "<<u1<<" Son: "<<u2<<" frequency: "<<_g_edges_reads[u1][pos]<<endl;*/
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

vector<vector<OwnNode_t>> DBG::export_unitigs_basic(const vector <string> & sequence_map, bool show)
{
    vector<vector<vector<OwnNode_t>>> vector_unitigs(_g_nodes.size(),vector<vector<OwnNode_t>>()),
        vector_flows(_g_nodes.size(),vector<vector<OwnNode_t>>());
    vector<UG_Node> starting_points = _get_starting_nodes_basic();
    vector<bool> checked(_g_nodes.size(), false);
    /*
     * Debugging mode
     */
    std::ofstream extension_information, post_process_information;
    if (Parameters::get().debug)
    {
        extension_information.open("debug/extension_information.txt", std::ofstream::out);
        post_process_information.open("debug/post_information.txt", std::ofstream::out);
        extension_information << "Getting unitigs from APDBG"<<endl;
    }
    for (auto start_point: starting_points)
    {
        if (Parameters::get().debug)
            extension_information << "Extending unitig from: "<<start_point._id<<":"<<start_point._val
            <<": "<<Common::return_unitig(sequence_map, start_point._val)<<endl;
        bool behaviour = true;
        OwnNode_t n = INF;
        vector<OwnNode_t> unitig, flows;
        (vector_unitigs[start_point._id]).push_back(unitig);
        (vector_flows[start_point._id]).push_back(flows);
        size_t index = 0;
        _extension_basic(vector_unitigs, vector_flows, start_point, checked, extension_information,index, sequence_map,behaviour,n);
        if (Parameters::get().debug)
            extension_information << "End local extension"<<endl;
    }
    vector<vector<OwnNode_t>> final_unitigs;
    if (!Parameters::get().greedy) {
        for (size_t i = 0; i < vector_unitigs.size(); ++i) {
            for (size_t j = 0; j < vector_unitigs[i].size(); ++j) {
                vector <OwnNode_t> tmp_unitig;
                for (auto n: vector_unitigs[i][j]) {
                    if (Parameters::get().debug)
                        extension_information << _g_nodes[n]._id << ":" << _g_nodes[n]._val << " ";
                    tmp_unitig.push_back(_g_nodes[n]._val);
                }
                if (tmp_unitig.size()) {
                    final_unitigs.push_back(tmp_unitig);
                    if (Parameters::get().debug)
                        extension_information << endl;
                }
            }
        }
        if (Parameters::get().debug)
            extension_information <<"Finishing basic extension!"<< endl;
        extension_information.close();
        return final_unitigs;
    }
    return final_unitigs;
}

void DBG::_extension_basic(vector<vector<vector<OwnNode_t>>> & unitigs, vector<vector<vector<OwnNode_t>>> & flows,
        UG_Node node, vector<bool> & checked, std::ofstream & extension_information, size_t & index, const vector <string> & sequence_map,
        bool behaviour, OwnNode_t node_id_continuation)
{
    vector <OwnNode_t> neighs;
    OwnNode_t parent = node._id;
    if (node_id_continuation == INF)
    {
        if (!checked[node._id])
            neighs = getNeighbor(node._id);
        unitigs[node._id][index].push_back(node._id);
        flows[node._id][index].push_back(0.0);
        if (Parameters::get().debug)
            extension_information << "Unitig: " << node._id << ":" << node._val << "\t";
        checked[node._id] = true;
    } else {
        if (!checked[node_id_continuation])
            neighs = getNeighbor(node_id_continuation);
        if (Parameters::get().debug)
            extension_information << "Unitig: " << node._id << ":" << node._val <<"("<<
                neighs.size()<<") Incorporamos: "<<node_id_continuation<<" INDEX: "<<index<<" "<<
                unitigs[node._id][index].size()<< "\t";
        checked[node_id_continuation] = true;
        parent = node_id_continuation;
    }
    while (neighs.size() > 0)
    {
        if (neighs.size() != 1)
        {
            size_t length = 0;
            for (size_t i;  i < unitigs[node._id][index].size();++i)
            {
                auto u = unitigs[node._id][index][i];
                length += (Common::return_unitig(sequence_map,_g_nodes[u]._val).size()-Parameters::get().kmerSize + 1);
            }
            /*if (length > 500)
            {
                length = 0;
                for (size_t i;  i < unitigs[node._id][index].size();++i)
                {
                    auto u = unitigs[node._id][index][i];
                    length += (Common::return_unitig(sequence_map,_g_nodes[u]._val).size()-Parameters::get().kmerSize + 1);
                    cout << "Size: " << length <<" Unitig: "<<_g_nodes[u]._id
                        <<" VAL: "<<_g_nodes[u]._val<<" Flow: "<<flows[node._id][index][i]<<" Indegree (unitig): "<<in_degree(_g_nodes[u]._id)<<endl;
                    cout << "U: "<<Common::return_unitig(sequence_map,_g_nodes[u]._val)<<endl;
                }
                cout << "Val: " << _g_nodes[parent]._val << " Id: " << _g_nodes[parent]._id
                     << " Neighs: " << neighs.size() << " Behaviour: " << behaviour << " Pair-end size: "
                     << _g_nodes[parent]._paired_info.size()<<" Length: "<<length << endl;
                //cin.get();
            }
            else{
                cout << "You are not long enough"<<endl;
            }*/
            if (_g_nodes[parent]._paired_info.size() == 0)
            {
                vector<OwnNode_t> unitig_original = unitigs[node._id][index], original_flow = flows[node._id][index];
                for (size_t i = 0; i < neighs.size(); ++i)
                {
                    if (_g_nodes[neighs[i]]._paired_info.size() == 0)
                    {
                        if (i > 0){
                            unitigs[node._id].push_back(unitig_original);
                            flows[node._id].push_back(original_flow);
                            index++;
                        }
                        checked[parent] = true;
                        unitigs[node._id][index].push_back(neighs[i]);
                        flows[node._id][index].push_back(_g_edges_reads[node._id][i]);
                        _extension_basic(unitigs, flows, node, checked, extension_information, index, sequence_map,behaviour, neighs[i]);
                    }
                }
            }
            break;
        }
        flows[node._id][index].push_back(_g_edges_reads[parent][0]);
        parent = _g_nodes[neighs[0]]._id;
        behaviour &= (in_degree(parent) <= 1);
        size_t parent_val = _g_nodes[neighs[0]]._val;
        if (!checked[parent]){
            neighs = getNeighbor(parent);
            unitigs[node._id][index].push_back(parent);
            if (Parameters::get().debug)
                extension_information<<parent<<":"<<parent_val<<"("<<neighs.size()<<") ";
            if (getNeighbor(parent, false).size() > 1)
            {
                for (auto u: unitigs[node._id][index]){
                    cout << _g_nodes[u]._val<<" -> ";
                }
                cout << "End" << endl;
                //cin.get();
            }
            checked[parent] = true;
        } else {
            if (find(unitigs[node._id][index].begin(), unitigs[node._id][index].end(), parent) !=  unitigs[node._id][index].end()){
                cout << " Index: "<<index<<endl;
                unitigs[node._id].erase(unitigs[node._id].begin() + index);
                flows[node._id].erase(flows[node._id].begin() + index);
                if (index > 0)
                    index--;
            }
            neighs.clear();
        }
    }
    /*cout << "Index: " <<index << endl;
    if (unitigs[node._id][index].size() > 5)
    {
        size_t length = 0;
        for (size_t i;  i < unitigs[node._id][index].size();++i)
        {
            auto u = unitigs[node._id][index][i];
            length += (Common::return_unitig(sequence_map,_g_nodes[u]._val).size()-Parameters::get().kmerSize + 1);
        }
        if (length > 500)
        {
            length = 0;
            for (size_t i;  i < unitigs[node._id][index].size();++i)
            {
                auto u = unitigs[node._id][index][i];
                length += (Common::return_unitig(sequence_map,_g_nodes[u]._val).size()-Parameters::get().kmerSize + 1);
                cout << "Size: " << length <<" Unitig: "<<_g_nodes[u]._id
                    <<" VAL: "<<_g_nodes[u]._val<<" Flow: "<<flows[node._id][index][i]<<" Indegree (unitig): "<<in_degree(_g_nodes[u]._id)<<endl;
                cout << "U: "<<Common::return_unitig(sequence_map,_g_nodes[u]._val)<<endl;
            }
            cout << "Val: " << _g_nodes[parent]._val << " Id: " << _g_nodes[parent]._id
                    << " Neighs: " << neighs.size() << " Behaviour: " << behaviour << " Pair-end size: "
                    << _g_nodes[parent]._paired_info.size()<<" Length: "<<length << endl;
            cin.get();
        }
    }*/
    if (Parameters::get().debug)
        extension_information << "END"<< endl;
}

vector<DBG::UG_Node> DBG::_get_potential_sinks_basic()
{
    vector<UG_Node> potential_sinks;
    for (auto i:_g_nodes)
    {
        if (!i._active)
            continue;
        if (out_degree(i._id) == 0)
            potential_sinks.push_back(i);
    }
    cout << "Potential sinks (basic): "<<potential_sinks.size()<<" Out of: "<<_g_nodes.size()<<endl;
    return potential_sinks;
}

vector<DBG::UG_Node> DBG::_get_starting_nodes_basic()
{
    /*
     * Debugging mode
     */
    std::ofstream zero_degree;
    if (Parameters::get().debug)
    {
        zero_degree.open("debug/zero_degree.txt", std::ofstream::out);
        zero_degree << "Nodes with out_degree equal to zero"<<endl;
    }
    vector<UG_Node> suspicious_nodes;
    for (auto i:_g_nodes)
    {
        if (!i._active)
            continue;
        if (in_degree(i._id) == 0 && i._active)
            suspicious_nodes.push_back(i);
        /*else if (out_degree(i._id) > 1)
            for (auto n:getNeighbor(i._id))
                suspicious_nodes.push_back(_g_nodes[n]);*/
        if (Parameters::get().debug)
        {
                zero_degree << "Node: "<<_g_nodes[i._id]._val<<" "<<_g_edges[i._id].size()<<" "<<i._active<<endl;
        }
    }
    cout << "Starting points (basic): "<<suspicious_nodes.size()<<" Out of: "<<_g_nodes.size()<<endl;
    return suspicious_nodes;
}

void DBG::_extension(queue <UG_Node> & queue, vector<bool> & checked, vector <OwnNode_t> & unitig,
        UG_Node node, size_t & higher_out, size_t & zero, size_t & higher_in, size_t & zero_nopair, size_t &solved,
        size_t & repetitive_call, vector<vector<OwnNode_t>> & unitigs, vector<bool> & local_checked, bool show, const vector<string> & sequence_map)
{
    /*
     * Revisar punto por punto la insercion de los unitigs
     */
    checked[node._id] = true;
    local_checked[node._id] = true;
    OwnNode_t parent = NO_NEIGH;
    vector<OwnNode_t> neighs = getNeighbor(node._id);
    bool extension = true;
    if (show)
    {
        cout << "Node: "<<node._id<<" "<<node._val<<" OutDegree: "<<neighs.size()<<" Indegree: "<<in_degree(node._id)<<endl;
        Common::display_unitig(sequence_map, node._val);cout<<endl;
        cout << "Neighbors: "<<endl;
        display_container(neighs);
        cout << "Paired_end: "<<endl;
        display_container(_g_nodes[node._id]._paired_info);
    }
    while (neighs.size() == 1)
    {
        //cout << "extending path "<<endl;
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
        if (show && unitig.size() > 3)
        {
            cout << "Unitig: "<<endl;
            for (auto u:unitig)
                cout << " "<<u;
            cout << endl;
            cout << "Output = 0" << endl;
            cout << "Node: " << node._id << " " << _g_nodes[node._id]._val << endl;
            Common::display_unitig(sequence_map, _g_nodes[node._id]._val);
            cout << "Paired-end: "<<endl;
            for (auto n: _g_nodes[node._id]._paired_info)
                cout << " "<<n;
            cout << endl;
            //cin.get();
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
                checked[node._id] = true;
                for (auto i:neighs) {
                    queue.push(_g_nodes[i]);
                }
            }
        } else {
            if (!checked[node._id]) {
                queue.push(_g_nodes[node._id]);
            }
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
                Common::display_unitig(sequence_map, _g_nodes[node._id]._val);
                cout << "Paired-end: "<<endl;
                for (auto n: _g_nodes[node._id]._paired_info)
                    cout << " "<<n;
                cout << endl;
                cout << "Neighbors (out):"<<endl;
                for (auto n:neighs)
                    cout << "Neigh (val): "<<_g_nodes[n]._val<<" "<<_g_nodes[n]._id<<endl;
                cout << endl;
                //cin.get();
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
                    cout << "Unitig pre: "<<endl;
                    display_container(unitig_saved);
                    _extension(queue, checked, unitig, _g_nodes[i],
                               higher_out, zero, higher_in, zero_nopair,solved,
                               repetitive_call, unitigs, local_checked, show,sequence_map);
                    cout << "Unitig post: "<<endl;
                    display_container(unitig);
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
        if (!i._active)
            continue;
        if (in_degree(i._id) == 0)
            suspicious_nodes.push_back(i);
        else if (out_degree(i._id) > 1) {
            /*for (auto n: getNeighbor(i._id))
                suspicious_nodes.push_back(n);*/
            suspicious_nodes.push_back(i);
        } else if (in_degree(i._id) > 1)
            suspicious_nodes.push_back(i);
    }
    cout << "Starting points: "<<suspicious_nodes.size()<<" Out of: "<<_g_nodes.size()<<endl;
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
        if (steps > D_MAX_PATH)
            return;
        if (branches == D_MAX_BRANCHES)
            return;
        /*if (src != target)
            this->_reachability[src].emplace(target);*/
        steps += getLength(target)-Parameters::get().kmerSize + 1;
        auto neighbors = getNeighbor(target);
        checked.emplace(target);
        branches += (neighbors.size() > 1)?1:0;
        for (auto neigh: neighbors)
        {
            this->_reachability[src].emplace(neigh);
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
    std::ofstream stats;
    stats.open ("debug/stats.txt",std::ofstream::out );
    float neighs = 0.0, pair_info_avg = 0.0;
    size_t num_nodes_no_pair = 0, num_nodes_pair = 0, no_pair_no_neigh = 0, max_freq = 0;
    vector<size_t> histogram_abundances, histogram;
    for (auto v:_g_nodes) {
        if (!v._active)
            continue;
        histogram_abundances.push_back(v._abundance);
        histogram.push_back(_g_nodes[v._id].get_paired_information().size());
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
            max_freq = (_g_nodes[v._id].get_paired_information().size() > max_freq)?
                    _g_nodes[v._id].get_paired_information().size():max_freq;
        }
    }
    stats << "Avg neighs: "<<(neighs/_g_nodes.size())<<endl;
    stats << "Avg paired end information: "<<(pair_info_avg/num_nodes_pair)<<endl;
    stats << "Number of nodes with no paired-end information "<<num_nodes_no_pair<<" out of "<<_g_nodes.size()<<endl;
    stats << "No pair and no neighbors: "<<no_pair_no_neigh<<endl;
    stats << "Max abundance: "<<_max_ab<<" Min abundance: "<<_min_ab<<endl;
    stats << "Global LCL: "<<_quantile_L<<" Global UCL: "<<_quantile_H<<endl;
    stats << "Median abundance: "<<_median_ab<<" Mean abundance: "<<_mean_ab<<endl;
    /*
     * Export number of pairs as histogram.
     */
    std::ofstream myfile;
    myfile.open ("histogram_pairs.txt", std::ofstream::out );
    for (auto p_f:histogram)
    {
        myfile<<","<<p_f;
    }
    myfile<<endl;
    myfile.close();
    myfile.open("abundance_histogram.txt", std::ofstream::out);
    for (auto p_f:histogram_abundances)
    {
        myfile<<","<<p_f;
    }
    myfile<<endl;
    myfile.close();
}

void DBG::subsane(const vector<string> & sequence_map)
{
    /*
     * Simple adjust to correct extremely rare cases
     */
    for (size_t i = 0; i < _g_nodes.size(); ++i)
    {
        if (!_g_nodes[i]._active)
            continue;
        if (_g_edges[i].size() == 1)
        {
            size_t freq_edge = _g_edges_reads[i][0], neigh = _g_edges[i][0];
            if (_g_nodes[i]._abundance > freq_edge && _g_nodes[neigh]._abundance > freq_edge)
                _g_edges_reads[i][0] = min(_g_nodes[i]._abundance,_g_nodes[neigh]._abundance);
        }
    }
    /*
     * Set frequencies for nodes based on frequencies of the edges
     */
    for (size_t i = 0; i < _g_edges.size(); ++i)
    {
        bool show = false;
        float node_freq = _g_nodes[i]._abundance;//_g_nodes_frequency[i];
        for (size_t j = 0; j < _g_edges[i].size(); ++j)
        {
            float neigh_freq = _g_nodes[_g_edges[i][j]]._abundance;//_g_nodes_frequency[_g_edges[i][j]];
            // ¿Max o min?
            //node_freq = max(neigh_freq, node_freq);
            if (show)
            {
                cout <<"Node: "<< _g_nodes[i]._val<<" "<<node_freq<<endl;
                cout << "Neighbor: "<<_g_nodes[_g_edges[i][j]]._val<<" Frequency: "<<neigh_freq<<endl;
                cout << "Edge frequency: "<<_g_edges_reads[i][j]<<endl;
            }
            node_freq = min(neigh_freq, node_freq);
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
     * Repolish
     */
    polish();
    /*
    * Getting full reachability matrix - after polishing the matrix
    */
    _complete_reach_matrix();
}

void DBG::export_to_gfa(const vector<string> & sequence_map, string file_name)
{
    std::ofstream outfile(file_name, std::ofstream::binary);
    outfile << "H\tVN:Z:1.0"<<endl;
    for (auto v:_g_nodes) {
        if (v._active)
            outfile << "S\t" << v._id << "\t" << Common::return_unitig(sequence_map, v._val) << endl;
    }
    size_t id = 0;
    for (auto n: _g_edges)
    {
        if (_g_nodes[id]._active)
            for (auto neigh:n)
                outfile << "L\t"<<id<<"\t+\t"<<neigh<<"\t+\t0M"<<endl;
        id++;
    }
    outfile << endl;
    outfile.close();
}

void DBG::print(OwnNode_t parent, OwnNode_t son, string file_name)
{
    if (parent == INF && son == INF) {
        std::ofstream outfile(file_name, std::ofstream::binary);
        for (auto v:_g_nodes) {
            if (!v._active)
                continue;
            outfile << "Node: " << v._id << "/" << v._val;
            outfile << " Frequency: " << _g_nodes_frequency[v._id] << " - "<<v._abundance<<" - ";
            outfile << " Number of neighbors: " << _g_edges[v._id].size() << " - ";
            for (size_t i = 0; i < _g_edges[v._id].size();++i)
            {
                OwnNode_t n = _g_edges[v._id][i];
                if (!_g_nodes[n]._active)
                    continue;
                outfile << " " << n << " " << _g_edges_reads[v._id][i] << "; ";
            }
            outfile << endl;
            outfile << "Paired-end information: ";
            for (auto n: _g_nodes[v._id].get_paired_information())
                outfile << " " << n << "/" << _g_nodes[v._id].get_abundance(n);
            outfile << endl;
        }
        outfile << "End print" << endl;
        outfile.close();
    }
    if (parent == 1)
    {
        std::ofstream outfile(file_name, std::ofstream::binary);
        for (auto v:_g_nodes) {
            outfile << "Node: " << v._id<<":" << v._val << " - "<<v._abundance<<" - ";
            outfile << " Number of neighbors: " << _g_edges[v._id].size() << " - ";
            for (size_t i = 0; i < _g_edges[v._id].size();++i)
            {
                OwnNode_t n = _g_edges[v._id][i];
                if (!_g_nodes[n]._active)
                    continue;
                outfile << " " << n << " " << _g_edges_reads[v._id][i] << "; ";
            }
            outfile << endl;
            outfile << "Paired-end information: ";
            for (auto n: _g_nodes[v._id].get_paired_information())
                outfile << " " << n;
            outfile << endl;
        }
        outfile << "End graph"<<endl;
        outfile.close();
    }
}

/*
 * Undirected Graphs
 */
OwnNode_t UG::addVertex(OwnNode_t node)
{
    _g_nodes.push_back(Node(_num_vertex, node));
    _g_edges.push_back(vector<OwnNode_t>());
    _g_in_edges.push_back(vector<OwnNode_t>());
    _g_freqs.push_back(vector<size_t>());
    return _num_vertex++;
}

OwnNode_t UG::addVertexWithAbundances(OwnNode_t node, float ab)
{
    _g_nodes.push_back(Node(_num_vertex,node,ab));
    _g_edges.push_back(vector<OwnNode_t>());
    _g_in_edges.push_back(vector<OwnNode_t>());
    _g_freqs.push_back(vector<OwnNode_t>());
    return _num_vertex++;
}

bool UG::setAbundance(OwnNode_t node, size_t ab, bool id)
{
    if (id) {
        if (node < _g_nodes.size()) {
            _g_nodes[node]._abundance = ab;
            return true;
        }
        return false;
    } else {
        for (size_t i = 0; i < _g_nodes.size(); ++i) {
            if (_g_nodes[i]._val == node) {
                _g_nodes[i]._abundance = ab;
                return true;
            }
        }
        return false;
    }
}
/*
 * It returns if the node has been added or not
 */
void UG::post_process_cycles(vector<size_t> & reached_by, size_t & minimum)
{
    for (auto pos_edge:_edges_leftovers) {
        OwnNode_t i = _g_edges[pos_edge.second][pos_edge.first];
        if (reached_by[i] == minimum)
            minimum--;
        reached_by[i]--;
        _g_in_edges[i].erase(find(_g_in_edges[i].begin(),_g_in_edges[i].end(), pos_edge.second));
        _g_edges[pos_edge.second].erase(_g_edges[pos_edge.second].begin() + pos_edge.first);
        _g_freqs[pos_edge.second].erase(_g_freqs[pos_edge.second].begin() + pos_edge.first);
    }
}

void UG::readjust_flow(float freq_strain, float & min_flow, float & max_flow)
{
    /*
     * Initial Flow readjustment
     */
    for (size_t i = 0; i < _g_nodes.size(); ++i)
    {
        float node_ab = abs(_g_nodes[i]._abundance - freq_strain);
        _g_nodes[i].setAbundance(node_ab);
        max_flow = (node_ab > max_flow)?node_ab:max_flow;
        min_flow = (node_ab < min_flow)?node_ab:min_flow;
        for (size_t j = 0; j < _g_freqs[i].size(); ++j)
        {
            float edge_ab = abs(_g_freqs[i][j] - freq_strain);
            _g_freqs[i][j] = edge_ab;
            max_flow = (edge_ab > max_flow)?edge_ab:max_flow;
            min_flow = (edge_ab < min_flow)?edge_ab:min_flow;
        }
    }
    /*
     * Potential Correction
     */
    for (size_t i = 0; i < _g_freqs.size(); ++i)
    {
        for (size_t j = 0; j < _g_freqs[i].size(); ++j) {
            //_g_freqs[i][j] = (max_flow >= _g_freqs[i][j])?max_flow - _g_freqs[i][j] + min_flow:min_flow;
            //_g_freqs[i][j] = ceil((MAX_GRANULARITY_ALLOWED*100)*abs((float)_g_freqs[i][j] - max_flow)/max_flow + MAX_GRANULARITY_ALLOWED);
            _g_freqs[i][j] = ceil(freq_strain*abs((float)_g_freqs[i][j] - max_flow)/((max_flow==0)?1:max_flow) + MAX_GRANULARITY_ALLOWED);
        }
        _g_nodes[i].readjustAbundance(freq_strain, min_flow, max_flow);
    }
}
/*
 * Flow in nodes is readjusted if the in_flow and out_flow are close to each other while node abundance is higher or lower
 *       - In the middle - makes sense
 *       - Lower or higher than both - correct (even in unary paths, because the edge from A - B must be already corrected remember that)
 */
void UG::correct_graph_safe(OwnNode_t starting_point, std::ofstream & adjust_information, vector<bool> & rescoring, DBG & dbg)
{
    bool show = Parameters::get().debug;
    // Traverse the graph and store paired_end information
    unordered_set<OwnNode_t> __nodes_val;
    for (auto n: _g_nodes)
        __nodes_val.emplace(n._val);
    for (auto n: _g_nodes)
    {
        auto pe = dbg.getPairedEndInformation(n._val);
        for (auto p:pe)
        {
            if (__nodes_val.find(p) != __nodes_val.end())
            {
                if (_paired_with.find(p) == _paired_with.end())
                    _paired_with[p] = unordered_set<OwnNode_t>();
                _paired_with[p].emplace(n._val);
                if (show)
                    adjust_information << "Node: "<<n._val<<" paired with: "<<p<<endl;
            }
        }
    }
    vector<bool> reach_matrix = vector<bool>(_g_nodes.size()*_g_nodes.size(), false);
    auto __complete_reach_matrix = [this, &reach_matrix]()
    {
        size_t num_nodes = _g_nodes.size();
        for (size_t i = 0; i < _g_nodes.size(); ++i)
        {
            for (size_t j = 0; j <  _g_edges[i].size();++j)
                reach_matrix[i*num_nodes + _g_edges[i][j]] = true;
        }
        for (size_t k = 0; k < _g_nodes.size(); ++k)
        {
            for (size_t i = 0; i < _g_nodes.size(); ++i)
            {
                for (size_t j = 0; j < _g_nodes.size(); ++j)
                    if (reach_matrix[i*num_nodes + k] & reach_matrix[k*num_nodes + j])
                        reach_matrix[i*num_nodes + j] = true;
            }
        }
    };
    auto _remove_transitive_relations = [this, &reach_matrix,&adjust_information,&show](OwnNode_t starting_point)
    {
        size_t num_nodes = _g_nodes.size();
        vector<vector<bool>> nodes_to_remove(_g_nodes.size(), vector<bool>(num_nodes,false));
        queue<OwnNode_t> nodes_to_process;
        nodes_to_process.push(starting_point);
        if (show)
            adjust_information << "Number of nodes: "<<num_nodes<<endl;
        while (!nodes_to_process.empty())
        {
            OwnNode_t node = nodes_to_process.front();
            nodes_to_process.pop();
            if (show)
                adjust_information << "Node: "<<node<<endl;
            auto neighs = getNeighbors(node);
            for (size_t i = 0; i < neighs.size(); ++i)
            {
                nodes_to_process.push(neighs[i]);
                for (size_t j = 0; j < neighs.size(); ++j)
                {
                    if (reach_matrix[neighs[i] * num_nodes + neighs[j]]) {
                        if (show)
                            adjust_information << "Edge to remove: "<<node<<" - "<<neighs[j]<<endl;
                        nodes_to_remove[node][j] = true;
                    }
                }
            }
        }
        /*
         * Counting on the order of the set - study case 2081:2048 (15-ZIKA-virus)
         */
        for (size_t i = 0; i < nodes_to_remove.size(); ++i)
        {
            size_t removals = 0;
            for (size_t j = 0; j < nodes_to_remove[i].size(); ++j)
            {
                if (nodes_to_remove[i][j]) {
                    OwnNode_t neigh_pos = _g_edges[i][j - removals];
                    if (show)
                        adjust_information << "Edge removed: " << i << " - " << neigh_pos << endl;
                    _g_in_edges[neigh_pos].erase(find(_g_in_edges[neigh_pos].begin(), _g_in_edges[neigh_pos].end(), i));
                    _g_edges[i].erase(_g_edges[i].begin() + j - removals);
                    _g_freqs[i].erase(_g_freqs[i].begin() + j - removals++);
                }
            }
        }
    };
    auto check_unary_path = [this,&rescoring](float freq_ref,float  & flow,OwnNode_t  node,
                OwnNode_t & node_end, float strain_scale, float & score)
    {
        auto neighs = getNeighbors(node);
        score+= rescoring[node];
        float max_diff = abs(freq_ref - flow), length = 1;
        while (neighs.size() == 1)
        {
            length++;
            score += rescoring[neighs[0]];
            float freq = strain_scale* (float)_g_freqs[node][0];
            /*
             * Change from abs(freq - freq_ref) < min_flow to abs(freq - freq_ref) > min_flow - it shows that in the unary path there are some inconsistences
             */
            if (abs(freq - freq_ref) > max_diff)
                flow = freq;
            node_end = neighs[0];
            neighs = getNeighbors(neighs[0]);
        }
        score /= length;
    };
    auto change_freqs_unary_path = [this](size_t new_freq, OwnNode_t node, OwnNode_t last_node)
    {
        while (node != last_node)
        {
            _g_freqs[node][0] = new_freq;
            node = _g_edges[node][0];
        }
    };
    std::function<void(OwnNode_t, float, bool)> adjust_flow =
            [this, &check_unary_path, &change_freqs_unary_path, &adjust_information, &adjust_flow,&show]
            (OwnNode_t node, float strain_scale_factor, bool skip)
    {
        vector<bool> adjusted = vector<bool>(_g_nodes.size(), false);
        /*
         * nodes process - nodes still remains unprocessed
         * strain scale factor - for every strain we store the scale for that precise strain
         * prev freq - for cases where in_degree & out_degree > 1 we need to know from where are we coming
         */
        unordered_set<OwnNode_t> edges_processed;
        queue<OwnNode_t> queue_nodes_process;
        queue<float> queue_strain_scale_factor, queue_prev_freq;
        queue_nodes_process.push(node);
        queue_strain_scale_factor.push(strain_scale_factor);
        queue_prev_freq.push(0);

        while (!queue_nodes_process.empty()) {
            OwnNode_t i = queue_nodes_process.front();
            strain_scale_factor = queue_strain_scale_factor.front();
            float prev_freq = queue_prev_freq.front();
            queue_nodes_process.pop();
            queue_strain_scale_factor.pop();
            queue_prev_freq.pop();
            UG_Node n = _g_nodes[i];
            /*
             * We check the BF traversion
             */
            bool all_traversed = true;
            for (size_t j = 0; j < _g_in_edges[n._id].size(); ++j)
                all_traversed &= (adjusted[_g_in_edges[n._id][j]]);
            if (!all_traversed)
            {
                queue_nodes_process.push(i);
                queue_strain_scale_factor.push(strain_scale_factor);
                queue_prev_freq.push(prev_freq);
                continue;
            }
            adjusted[n._id] = true;
            vector<float> scale_factors = vector<float>(_g_edges[n._id].size(), strain_scale_factor);
            /*
             * Neighs_to_process stores the indexes to check
             */
            std::function<void(OwnNode_t, vector<OwnNode_t>&, float, float)> get_proper_neighs
                                   = [this](OwnNode_t node,vector<OwnNode_t> & neighbors,
                                           float prev_freq, float scale_factor)
            {
                vector<float> edge_freqs;
                for (size_t i = 0; i < neighbors.size(); ++i)
                    edge_freqs.push_back(_g_freqs[node][i]*scale_factor);
                sort(edge_freqs.begin(), edge_freqs.end());
                vector<OwnNode_t> neighs_selected;
                for (int i = neighbors.size(); i >= 0; --i)
                {

                }
            };
            vector<OwnNode_t> neighs_to_process = getNeighbors(n._id);
            /*
             * MOD!
             */
            /*if (neighs_to_process.size() > 1)
                skip = true;*/
            if (_g_in_edges[n._id].size() > 1 && _g_edges[n._id].size() > 1)
            {
                strain_scale_factor = 1;
                if (show)
                    adjust_information << "Complex case, prev freq: "<< prev_freq<<" Scale factor (switched to 1): "<<strain_scale_factor<<" OutNeighs: "
                        <<_g_edges[n._id].size()<<" InNeighs: "<<_g_in_edges[n._id].size()<<endl;
                for (size_t i = 0; i < neighs_to_process.size();++i)
                    neighs_to_process[i] = i;
            } else {
                for (size_t i = 0; i < neighs_to_process.size();++i)
                    neighs_to_process[i] = i;
            }
            if (!skip) {
                float out_flow = 0;
                for (size_t j = 0; j < _g_freqs[n._id].size(); ++j)
                    out_flow += ((float) _g_freqs[n._id][j]);
                float in_flow = 0;
                for (size_t j = 0; j < _g_in_edges[n._id].size(); ++j) {
                    OwnNode_t neigh_id = _g_in_edges[n._id][j];
                    in_flow += _g_freqs[neigh_id][find(_g_edges[neigh_id].begin(), _g_edges[neigh_id].end(), n._id) -
                                                  _g_edges[neigh_id].begin()];
                }
                if (show)
                    adjust_information << "Node: (id) " << _g_nodes[n._id]._id <<
                                       " Node: (val) " << _g_nodes[n._id]._val << " Inflow: " << in_flow
                                       << " Outflow: " << out_flow << " Scale strain: " << strain_scale_factor << endl;
                if ((out_flow > (1 + FLOW_LIMIT_RATIO) * in_flow) || (out_flow * (1 + FLOW_LIMIT_RATIO) < in_flow)) {
                    if (show) {
                        adjust_information << "Node out of flow: "
                                              "node (id): " << n._id << " Node (val): " << n._val
                                           << " Difference: " << abs((float) in_flow - (float) out_flow) << endl;
                    }
                    size_t sum_mins = 0, sum_trads = 0;
                    vector <OwnNode_t> last_nodes = vector<OwnNode_t>(neighs_to_process.size(), false),
                            min_flows_unary = vector<size_t>(neighs_to_process.size(), 0);
                    vector<float> scores = vector<float>(neighs_to_process.size(), 0.0);
                    for (size_t j = 0; j < neighs_to_process.size(); ++j) {
                        OwnNode_t neigh = _g_edges[n._id][neighs_to_process[j]];
                        OwnNode_t last_node = _g_nodes[neigh]._id;
                        float flow = strain_scale_factor * (float) _g_freqs[n._id][neighs_to_process[j]], score = 0.0;
                        adjust_information << "Checking: " << last_node << " flow: " << flow << endl;
                        check_unary_path(in_flow, flow, _g_nodes[neigh]._id, last_node,
                                         strain_scale_factor, score);
                        adjust_information << "New flow: " << flow << " Until: " << last_node << " Edge flow: "
                                           << neigh << " With score: " << score << endl;
                        scores[j] = score;
                        min_flows_unary[j] = abs(in_flow - flow);
                        last_nodes[j] = last_node;
                        sum_mins += min_flows_unary[j];
                    }
                    /*
                     * Standarized scoring
                     * Change: 
                     */
                    float sum_scores = Common::sum_vector(scores);
                    for (size_t j = 0; j < scores.size(); ++j)
                        /*scores[j] = (sum_scores == 0)?1:scores[j]/sum_scores;*/
                        scores[j] = (scores[j] == 0)?0:1;
                    /*
                     * Testing non-linear behaviour
                     */
                    for (size_t j = 0; j < min_flows_unary.size(); ++j) {
                        min_flows_unary[j] = pow(sum_mins - min_flows_unary[j], 1) * scores[j];
                        sum_trads += min_flows_unary[j];
                        adjust_information << " Val ori: " << _g_freqs[n._id][neighs_to_process[j]] << " Val trad: " << min_flows_unary[j]
                                           << endl;
                    }
                    adjust_information << "Sum trads: " << sum_trads << endl;
                    for (size_t j = 0; j < min_flows_unary.size(); ++j) {
                        OwnNode_t neigh = neighs_to_process[j];
                        // Removed round
                        size_t new_flow = round((min_flows_unary[j] == 0) ?
                                          (sum_trads == 0) ? ((float) in_flow / (float) min_flows_unary.size())
                                                           : MIN_FLOW_PATH :
                                          ((float) min_flows_unary[j] / (float) sum_trads) * (float) in_flow);
                        scale_factors[j] = ((float) new_flow / (float) _g_freqs[n._id][neigh]);
                        adjust_information << "Previous flow: " << _g_freqs[n._id][neigh] << endl;
                        _g_freqs[n._id][neigh] = new_flow;
                        //change_freqs_unary_path(new_flow, _g_nodes[_g_edges[n._id][j]]._id, last_nodes[j]);
                        if (show) {
                            adjust_information << "Flow changed for edge (val1,val2): " << _g_nodes[n._id]._val << " "
                                               << _g_nodes[neigh]._val
                                               << " New freq flow: " << new_flow << " Previous flow: "
                                               << new_flow / scale_factors[j] << endl;
                        }
                    }
                }
            }
            skip = false;
            for (size_t j = 0; j < neighs_to_process.size(); ++j) {
                /*
                 * We only add an edge once
                 */
                OwnNode_t edge = ((n._id << 15) | (_g_edges[n._id][neighs_to_process[j]]));
                if (show)
                    adjust_information<<"EDGE: "<<edge<<endl;
                if (edges_processed.find(edge) == edges_processed.end()) {
                    if (show)
                        adjust_information << " Adding neigh: " << _g_edges[n._id][neighs_to_process[j]]
                                           << " With scale: " << scale_factors[j] << " With freq: "
                                           << _g_freqs[n._id][neighs_to_process[j]] << endl;
                    queue_nodes_process.push(_g_edges[n._id][neighs_to_process[j]]);
                    queue_strain_scale_factor.push(scale_factors[j]);
                    queue_prev_freq.push(_g_freqs[n._id][neighs_to_process[j]]);
                    edges_processed.emplace(edge);
                }
            }
        }
    };
    if (show)
    {
        adjust_information << "First step - Correcting transitive relations in the DAG."<<endl;
    }
    __complete_reach_matrix();
    _remove_transitive_relations(starting_point);
    if (show) {
        adjust_information << "Adjusting flow discrepancies in the pre-network flow graph" << endl;
        adjust_information << "Adjusting flow from: " << starting_point << endl;
    }
    adjust_flow(starting_point,1.0, true);
    if (show)
        adjust_information << "Adjusting flow in nodes. "<<endl;
    for (size_t i = 0; i < _g_nodes.size(); ++i) {
        Node n = _g_nodes[i];
        /*
         * Isolated nodes does not make sense
         */
        if (in_degree(n._id) == 0 && out_degree(n._id) == 0)
        {
            _g_nodes[n._id]._abundance = 1;
            adjust_information << "Isolated node turn abundance to 1"<<endl;
            continue;
        }
        /*
         * For in-between nodes
         */
        if (_g_freqs[i].size() && _g_in_edges[i].size()) {
            size_t out_flow = 0;
            for (size_t j = 0; j < _g_freqs[i].size(); ++j)
                out_flow += _g_freqs[i][j];
            size_t in_flow = 0;
            for (size_t j = 0; j < _g_in_edges[i].size(); ++j) {
                OwnNode_t neigh_id = _g_in_edges[i][j];
                in_flow += _g_freqs[neigh_id][find(_g_edges[neigh_id].begin(), _g_edges[neigh_id].end(), n._id) -
                                              _g_edges[neigh_id].begin()];
            }
            if ((n._abundance > out_flow && n._abundance > in_flow) ||
                (n._abundance < out_flow && n._abundance < in_flow)) {
                _g_nodes[i]._abundance = min(in_flow, out_flow);
                if (show)
                    adjust_information << "Node: (id) " << _g_nodes[i]._id << " (val) " << _g_nodes[i]._val
                                       << "CORRECTION! " << _g_nodes[i]._abundance << endl;
            }
        }
        /*
         * For sinks
         */
        if (_g_freqs[i].size() == 0 && _g_in_edges[i].size())
        {
            size_t in_flow = 0;
            for (size_t j = 0; j < _g_in_edges[i].size(); ++j) {
                OwnNode_t neigh_id = _g_in_edges[i][j];
                in_flow += _g_freqs[neigh_id][find(_g_edges[neigh_id].begin(), _g_edges[neigh_id].end(), n._id) -
                                              _g_edges[neigh_id].begin()];
            }
            if ((n._abundance > in_flow) || (n._abundance < in_flow)) {
                _g_nodes[i]._abundance = in_flow;
                if (show)
                    adjust_information << "Node: (id) " << _g_nodes[i]._id << " (val) " << _g_nodes[i]._val
                                       << "CORRECTION! " << _g_nodes[i]._abundance << endl;
            }
        }
    }
    if (show)
        adjust_information << "End flow adjustment."<<endl;
}

void UG::correct_graph(std::ofstream & adjust_information)
{
    if (Parameters::get().debug)
        adjust_information << "Adjusting flow discrepancies pre-network flow graph"<<endl;
    for (size_t i = 0; i < _g_nodes.size(); ++i)
    {
        if (_g_freqs[i].size())
        {
            Node n = _g_nodes[i];
            size_t out_flow = 0;
            for (size_t j = 0; j < _g_freqs[i].size(); ++j)
                out_flow += _g_freqs[i][j];
            if (Parameters::get().debug) {
                adjust_information << "Node: " << n._id << ":" << n._val << "-" << n._abundance << endl;
                adjust_information << "OutFlow: " << out_flow << " MaxAllowed: "
                                   << ((float) out_flow * (1 + CORRECT_RATIO)) << endl;
            }
            if ((n._abundance > (float) out_flow * (1 + CORRECT_RATIO)) || (n._abundance*(1 + CORRECT_RATIO) < (float) out_flow))
            {
                if (Parameters::get().debug)
                    adjust_information << "Abundance from "<<n._abundance<<" to "<<out_flow<<endl;
                _g_nodes[i]._abundance = out_flow;
            }
        } else {
            Node n = _g_nodes[i];
            size_t in_flow = 0;
            for (size_t j = 0; j < _g_in_edges[i].size();++j)
            {
                OwnNode_t neigh_id = _g_in_edges[i][j];
                in_flow += _g_freqs[neigh_id][find(_g_edges[neigh_id].begin(), _g_edges[neigh_id].end(),n._id)-_g_edges[neigh_id].begin()];
            }
            if (((n._abundance > (float) in_flow * (1 + CORRECT_RATIO)) || (n._abundance*(1 + CORRECT_RATIO) < (float) in_flow)) & in_flow != 0)
            {
                if (Parameters::get().debug)
                    adjust_information << "Abundance from "<<n._abundance<<" to "<<in_flow<<endl;
                _g_nodes[i]._abundance = in_flow;
            }
        }
    }
    if (Parameters::get().debug)
        adjust_information << "End adjust."<<endl;
}

/*
    Pensar algo para este optimal cost.
*/
void UG::addEdge(OwnNode_t i, OwnNode_t j, size_t cost, size_t optimal_cost)
{
    vector<OwnNode_t>::iterator pos_j_in_i = find(_g_edges[i].begin(), _g_edges[i].end(), j);
    if ( pos_j_in_i == _g_edges[i].end())
    {
        _num_edges++;
        if (!_directed) {
            _g_edges[i].push_back(j);
            _g_edges[j].push_back(i);
            _g_freqs[i].push_back(cost);
            _g_freqs[j].push_back(cost);
        } else {
            /*
             * Check frequencies - keep the highest one when A -> B and B -> A (if A reachs B but previously B had reached A)
             */
            bool remove = false, add = true;
            size_t index = 0;
            vector<OwnNode_t >::iterator pos = std::find (_g_edges[j].begin(), _g_edges[j].end(), i);
            if (pos != _g_edges[j].end()) {
                index = std::distance( _g_edges[j].begin(), pos);
                remove = (_g_freqs[j][index] < cost);
                if (remove) {
                    _edges_leftovers.push_back({index, j});
                } else
                    add = false;
            }
            if (add)
            {
                _g_edges[i].push_back(j);
                _g_freqs[i].push_back(cost);
                _g_in_edges[j].push_back(i);
            }
        }
    } else {
        /*
         * When j has already been added to i
         */
        if (STRATEGY == 3)
        {
            _g_freqs[i][distance(_g_edges[i].begin(), pos_j_in_i)] = cost;
        } else {
            size_t index = std::distance(_g_edges[i].begin(), pos_j_in_i);
            float dif_actual = abs((float)_g_freqs[i][index] - (float) optimal_cost),
                dif_potencial = abs((float) cost - (float) optimal_cost), dif_add = abs((float)_g_freqs[i][index] + (float) cost - (float) optimal_cost);
            float min_diff = min(min(dif_actual, dif_potencial), dif_add);
            if (min_diff == dif_potencial)
                _g_freqs[i][index] = cost;
            else if (min_diff == dif_add)
                _g_freqs[i][index] += cost;
        }
    }
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

size_t UG::getFreqs(OwnNode_t node, size_t index)
{
    return _g_freqs[node][index];
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

priority_queue<pair<size_t,vector<OwnNode_t>>> UG::report_min_cost_flow(const vector<size_t> & reached_by,
        size_t strain_freq,bool show, unordered_map<size_t, bool> true_nodes, string file_name, vector<bool> & rescoring)
{
    /*
     * Debugging files
     */
    std::ofstream nw_file;
    if (show)
        nw_file.open (file_name,std::ofstream::out);
    /*
     * Define maximal allowed granularity
     */
    float max_granularity = std::max(std::min((float)MAX_GRANULARITY_ALLOWED, (float)strain_freq),(float)1.0);
    /*
     * Result paths
     */
    priority_queue<pair<size_t,vector<OwnNode_t>>> result_paths;
    unordered_map<OwnNode_t, vector<float>> pre_flow_map;
    /*
     * Truncate demanded flow
     */
    //strain_freq = Maths::my_round(strain_freq, true);
    Graph l_g, solution_g;
    /*
     * Solution network containers
     */
    Graph::ArcMap<Capacity> available_flow(solution_g);
    Graph::ArcMap<Graph::Arc> offsetArcToSolutionArc(l_g);
    Graph::ArcMap<int> arcDirection(l_g);
    /*
     * Supply map
     */
    Graph::NodeMap<float> supplyMap(l_g);
    vector<Graph::Node> out_nodes;
    vector<std::pair<Graph::Node, Graph::Node>> map_in_out;
    Graph::NodeMap<OwnNode_t> translation_map(l_g), solution_to_ids(solution_g);
    vector<Graph::Node> offset_to_solution(_num_vertex+2);
    Graph::Node source = l_g.addNode(), target = l_g.addNode(),
                s_source = l_g.addNode(), s_target = l_g.addNode();
    /*
     * Extra_ids
     */
    size_t extra_id = reached_by.size();
    size_t id_source = extra_id++, id_target = extra_id++, id_s_source = extra_id++, id_s_target = extra_id++;
    // Build the source/sink and synthetic sources and sinks - to avoid infinite and insufficient
    translation_map[source] = id_source; translation_map[target] = id_target;
    //translation_map[g_source] = INF; translation_map[g_target] = INF;
    translation_map[s_source] = id_s_source; translation_map[s_target] = id_s_target;
    if (show)
    {
        nw_file << "Global source: "<<id_source<<" Global target: "<<id_target<<endl;
        nw_file << "Start source: "<<id_s_source<<" Start target: "<<id_s_target<<endl;
    }
    /*
     * Out flow and in flow
     */
    unordered_map<OwnNode_t, float> in_flow, out_flow;
    in_flow[id_source] = 0.0;out_flow[id_source] = 0.0;
    in_flow[id_target] = 0.0;out_flow[id_target] = 0.0;
    /*
     * Arc capacities and weights
     */
    vector<Graph::Arc> arc_vector;
    ArcMap<Capacity> capacities(l_g), lowerMap(l_g);
    ArcMap<Weight> weights(l_g);
    ArcMap<bool> true_arcs(l_g);
    /*
     * Set of sources and targets
     */
    unordered_set<OwnNode_t> sources_set, targets_set;
    float flow_from_target = 0;
    for (auto n:_g_nodes)
    {
        if (_g_in_edges[n._id].size() == 0) {
            sources_set.emplace(n._id);
        }
        if (_g_edges[n._id].size() == 0) {
            flow_from_target += _g_nodes[n._id]._abundance;
            targets_set.emplace(n._id);
        }
    }
    for (auto n:_g_nodes)
    {
        if (show)
            nw_file << "Node: "<<n._id<<" Val: "<<n._val<<endl;
        /*
         * Initialize flow map
         */
        pre_flow_map[n._id] = vector<float>(_g_nodes.size(), 0.0);
        /*
         * Splitting nodes
         */
        Graph::Node node = solution_g.addNode();
        offset_to_solution[n._id] = node;
        solution_to_ids[node] = n._id;
        Graph::Node in_node = l_g.addNode(), out_node = l_g.addNode();
        map_in_out.push_back(pair<Graph::Node,Graph::Node>(in_node, out_node));
        translation_map[in_node] = n._id;translation_map[out_node] = n._id;
        out_nodes.push_back(out_node);
        /*
         * Fill supply map
         */
        supplyMap[in_node] = 0.0;supplyMap[out_node] = 0.0;
        in_flow[n._id] = 0.0;out_flow[n._id] = 0.0;
        //Arc forward - using split into several edges
        float pair_abundance =(READJUST)?((float) n._abundance):(abs((float)n._abundance - strain_freq));
        pair_abundance = (pair_abundance <= MAX_GRANULARITY_ALLOWED)?(MAX_GRANULARITY_ALLOWED+1):pair_abundance;
        float arc_increment = ceil((float) pair_abundance / (float) max_granularity),
            capacity = ceil((float) n._abundance / (float) max_granularity);
        if (show)
            nw_file << "Arc_increment: "<<arc_increment<<" Abundance: "<<n._abundance<<" Abundance corrected: "<<pair_abundance<<endl;
        float freq = arc_increment * max_granularity;
        for (float i = arc_increment; i <= freq; i += arc_increment)
        {
            Graph::Arc a = l_g.addArc(in_node, out_node);
            lowerMap[a] = 0;capacities[a] = capacity;true_arcs[a] = false;
            weights[a] = (i == arc_increment)?i/freq:pow((i-arc_increment),2)/freq;
            if (show)
                nw_file <<"Capacity: "<<capacities[a]<<" Weight: "<<weights[a]<<endl;
            arc_vector.push_back(a);
        }
        //Arc backward - using split into several edges
        for (float i = arc_increment; i <= freq; i += arc_increment)
        {
            Graph::Arc a = l_g.addArc(out_node, in_node);
            lowerMap[a] = 0;capacities[a] = capacity;true_arcs[a] = false;
            weights[a] = (i == arc_increment)?i/freq:pow((i-arc_increment),2)/freq;
            arc_vector.push_back(a);
        }
    }
    if (show)
        nw_file <<"End nodes arcs"<<endl;
    /*
     * Edge from sink to source
     */
    for (size_t i = 0; i < _g_edges.size(); ++i)
    {
        if (show)
            nw_file << "Node: "<<_g_nodes[i]._id<<endl;
        pair<Graph::Node, Graph::Node> n_1 = map_in_out[i];
        for (size_t j = 0; j < _g_edges[i].size(); ++j)
        {
            if (show)
                nw_file << "   Neighs (id): "<<_g_nodes[_g_edges[i][j]]._id<<" Val: "<<_g_nodes[_g_edges[i][j]]._val<<endl;
            auto n = _g_edges[i][j];
            pair<Graph::Node, Graph::Node> n_2 = map_in_out[n];

            float freq =(READJUST)?(float) _g_freqs[i][j]:(abs((float) _g_freqs[i][j] - strain_freq));
            // 2 is just a hot fix
            freq = (freq <= MAX_GRANULARITY_ALLOWED)?(MAX_GRANULARITY_ALLOWED+1):freq;
            /*
             * Solution Arc
             */
            Graph::Arc solution_arc = solution_g.addArc(offset_to_solution[i], offset_to_solution[n]);
            available_flow[solution_arc] = _g_freqs[i][j];
            float arc_increment = ceil((float) freq / (float) max_granularity),
                capacity = ceil((float) _g_freqs[i][j] / (float) max_granularity);
            if (show)
                nw_file << "Arc_increment: "<<arc_increment<<" Frequency: "<<freq<<endl;
            freq = arc_increment*max_granularity;
            // Forward flow
            for (float z = arc_increment; z <= freq; z += arc_increment)
            {
                Graph::Arc a = l_g.addArc(n_1.second, n_2.first);
                offsetArcToSolutionArc[a] = solution_arc;
                arcDirection[a] = 1;
                lowerMap[a] = 0;capacities[a] = capacity;true_arcs[a] = true;
                weights[a] = ((z == arc_increment)?z/freq:pow((z-arc_increment),2)/freq)+SECURE_OFFSET*j;
                if (show)
                    nw_file <<"Capacities: "<<capacities[a]<< " Weight: "<<weights[a]<<endl;
                arc_vector.push_back(a);
            }
            // Backward flow
            for (float z = arc_increment; z <= freq; z += arc_increment)
            {
                Graph::Arc a = l_g.addArc(n_2.first, n_1.second);
                offsetArcToSolutionArc[a] = solution_arc;
                arcDirection[a] = -1;
                lowerMap[a] = 0;capacities[a] = capacity;true_arcs[a] = true;
                weights[a] = ((z == arc_increment)?z/freq:pow((z-arc_increment),2)/freq)+SECURE_OFFSET*j;
                arc_vector.push_back(a);
            }
            // Update in and out flow
            out_flow[n] += _g_freqs[i][j];
            in_flow[i] += _g_freqs[i][j];
        }
        //float in_distribution = (float) Maths::my_round(in_flow[i], true) / (float) _g_edges[i].size();
        float in_distribution = round((float) in_flow[i] / (float) _g_edges[i].size());
        size_t diff = std::abs((float) _g_nodes[i]._abundance - (float) in_flow[i]);
        // Edge from new source to previous sources null cost and infinite capacity
        if (_g_in_edges[i].size() == 0)
        {
            //solution_g.addArc(solution_source,offset_to_solution[i]);
            // Forward Flow
            Graph::Arc a = l_g.addArc(source, n_1.first);
            lowerMap[a] = 0;capacities[a] = INF;weights[a] = 0;true_arcs[a] = false;
            arc_vector.push_back(a);
            // Backward Flow
            a = l_g.addArc(n_1.first, source);
            lowerMap[a] = 0; capacities[a] = INF; weights[a] = 0;true_arcs[a] = false;
            arc_vector.push_back(a);
        }
        // Edge from previous targets to new target with null cost and infinite capacity
        if (_g_edges[i].size() == 0 /*|| diff > in_distribution*/)
        {
            //solution_g.addArc(offset_to_solution[i], solution_target);
            // Forward Flow
            Graph::Arc a = l_g.addArc(n_1.second, target);
            lowerMap[a] = 0; capacities[a] = INF;weights[a] = 0;true_arcs[a] = false;
            arc_vector.push_back(a);
            // Backward Flow
            a = l_g.addArc(target, n_1.second);
            lowerMap[a] = 0; capacities[a] = INF; weights[a] = 0; true_arcs[a] = false;
            arc_vector.push_back(a);
        }
        if (show)
            nw_file << "Node: "<<_g_nodes[i]._val<<" Id: "<<i<<" Reached by: "<<_g_in_edges[i].size()<<" Neighs: "<<_g_edges[i].size()<<endl;
    }
    if (show)
        nw_file <<"End edges arcs"<<endl;
    /*
     * Check in_flow to potential sources:
     */
    /*for (auto o_flow: out_flow)
    {
        if (o_flow.first == id_source || o_flow.first == id_target)
            continue;
        float out_distribution = (float) Maths::my_round(o_flow.second) / (float) _g_edges[o_flow.first].size();
        size_t diff = std::abs((float) _g_nodes[o_flow.first]._abundance - (float) out_flow[o_flow.first]);
        if (diff > out_distribution && sources_set.find(o_flow.first) == sources_set.end())
        {
            pair<Graph::Node, Graph::Node> n_1 = map_in_out[o_flow.first];
            Graph::Arc a = l_g.addArc(source,n_1.first);
            capacities[a] = INF;
            weights[a] = 0;//1 / sqrt(_g_nodes[i]._abundance);
            arc_vector.push_back(a);
            sources_set.emplace(o_flow.first);
        }
    }*/
    /*
     * Sources and targets
     */
    if (show)
    {
        nw_file << "Sources set: " << endl;
        nw_file << set_to_string(sources_set);
        nw_file << "Targets set: " << endl;
        nw_file << set_to_string(targets_set);
    }
    /*
     * Readjust flow: remember - n_in -> n_out; flow_in goes to n_out while flow_out goes to n_in
     */
    float supply_source_s = 0.0, supply_sink_s = 0.0;
    for (auto k_v: in_flow)
    {
        if (k_v.first == id_source || k_v.first == id_target)
            continue;
        float ab = _g_nodes[k_v.first]._abundance;
        float exogenous_flow_in = k_v.second - ab, exogenous_flow_out = out_flow[k_v.first] - ab;
        if (show)
        {
            nw_file << "Node: "<<_g_nodes[k_v.first]._val<<" In flow: "<<k_v.second<<" Out flow: "<<out_flow[k_v.first]<<" Id: "<<k_v.first<<endl;
            nw_file << "Exogenous flow in: " << exogenous_flow_in << " Exogenous flow out: " << exogenous_flow_out << " Abundance: "<<ab<< endl;
        }
        if (exogenous_flow_in < 0)
        {
            if (show)
                nw_file << "    From s* to OutNode "<<(-exogenous_flow_in)<<endl;
            Graph::Arc a = l_g.addArc(s_source, map_in_out[k_v.first].second);true_arcs[a] = false;
            lowerMap[a] = -exogenous_flow_in;capacities[a] = - exogenous_flow_in;weights[a] = 0;
            supply_source_s += (-exogenous_flow_in);arc_vector.push_back(a);
        } else if (exogenous_flow_in > 0)
        {
            if (show)
                nw_file << "    From OutNode to t* "<<exogenous_flow_in<<endl;
            Graph::Arc a = l_g.addArc(map_in_out[k_v.first].second, s_target);true_arcs[a] = false;
            lowerMap[a] = exogenous_flow_in;capacities[a] = exogenous_flow_in;weights[a] = 0;
            supply_sink_s += exogenous_flow_in;arc_vector.push_back(a);
        }
        if (exogenous_flow_out < 0)
        {
            if (show)
                nw_file << "    From InNode to t* "<<(-exogenous_flow_out)<<endl;
            Graph::Arc a = l_g.addArc(map_in_out[k_v.first].first, s_target);true_arcs[a] = false;
            lowerMap[a] = -exogenous_flow_out;capacities[a] = -exogenous_flow_out;weights[a] = 0;
            supply_sink_s += (-exogenous_flow_out);arc_vector.push_back(a);
        } else if (exogenous_flow_out > 0)
        {
            if (show)
                nw_file << "    From s* to InNode "<<exogenous_flow_out<<endl;
            Graph::Arc a = l_g.addArc(s_source, map_in_out[k_v.first].first);true_arcs[a] = false;
            lowerMap[a] = exogenous_flow_out;capacities[a] = exogenous_flow_out;weights[a] = 0;
            supply_source_s += exogenous_flow_out;arc_vector.push_back(a);
        }
    }
    /*
     * Exogenous flow from sources to s_star and target to s_star
     */
    /*Graph::Arc a = l_g.addArc(source, s_target);
    lowerMap[a] = in_flow[id_source];capacities[a] = in_flow[id_source];weights[a] = 0;true_arcs[a] = false;
    supply_sink_s += in_flow[id_source];
    a = l_g.addArc(s_source, target);
    lowerMap[a] = out_flow[id_target];capacities[a] = out_flow[id_target];weights[a] = 0;true_arcs[a] = false;
    supply_source_s += out_flow[id_target];*/
    if (show)
        nw_file <<"End sources supply: (supply source*)"<<supply_source_s<<" (supply sink*)"<<supply_sink_s<<endl;
    /*
     * Fill the supplyMap with the flow in the start_source/sink
     */
    supplyMap[s_source] = supply_source_s;supplyMap[s_target] = -supply_sink_s;
    if (show) {
        nw_file << "From star source: "<<endl;
        for (Graph::OutArcIt a(l_g, s_source); a != INVALID; ++a)
        {
            nw_file << "From: "<<translation_map[l_g.source(a)]<<" to "<<translation_map[l_g.target(a)]<<" "<<capacities[a]<<endl;
        }
        nw_file << "To sink star: "<<endl;
        for (Graph::InArcIt a(l_g, s_target); a != INVALID; ++a)
        {
            nw_file << "From: "<<translation_map[l_g.source(a)]<<" to "<<translation_map[l_g.target(a)]<<" "<<capacities[a]<<endl;
        }
    }
    /*
     * Edge from sink to source
     */
    Graph::Arc a = l_g.addArc(target, source);
    lowerMap[a] = 0;capacities[a] = INF;weights[a] = 0;
    arc_vector.push_back(a);

    NS ns(l_g);
    if (show)
        nw_file << "Demanded flow: "<<strain_freq<<endl;
    ns.costMap(weights).lowerMap(lowerMap).upperMap(capacities).supplyMap(supplyMap);
    ArcMap<Capacity> flows(l_g);
    NS::ProblemType status = ns.run(NS::BLOCK_SEARCH);
    switch (status) {
        case NS::INFEASIBLE:
            cerr << "Infeasible flow" << endl;
            exit(1);
            break;
        case NS::OPTIMAL:
            ns.flowMap(flows);
            if (show)
                nw_file << "Cost: "<<ns.totalCost()<<endl;
            size_t num_edges = 0;
            for (Graph::ArcIt a(l_g); a != INVALID; ++a)
            {
                if (show)
                    nw_file << "Source: "<<translation_map[l_g.source(a)]<<" to: "
                        <<translation_map[l_g.target(a)]<<" Flow: "<<ns.flow(a)
                                <<" Solution edge: "<<true_arcs[a]<<" Direccion: "<<arcDirection[a]<<endl;
                if (true_arcs[a])
                    available_flow[offsetArcToSolutionArc[a]] += (ns.flow(a) * arcDirection[a]);
                num_edges++;
            }
            if (show)
                nw_file<<"Sanity edges: "<<num_edges<<" "<<arc_vector.size()<<endl;
            /*
             * Flow from solution source to sources
             */
            if (show)
                nw_file<<"Available flow:"<<endl;
            for (Graph::ArcIt n(solution_g); n != INVALID; ++n)
                if (show)
                    nw_file << solution_to_ids[solution_g.source(n)] <<" -> "<<solution_to_ids[solution_g.target(n)]<<" "<<available_flow[n]<<endl;
            if (show)
                nw_file << endl;
            if (show)
                nw_file << "End available flow"<<endl;
            std::function<void(Graph::Node, vector<Graph::Arc>&, OwnNode_t&, float&)> __flow_to_path =
                    [this,& __flow_to_path,& solution_g, &available_flow, &targets_set, &nw_file,&solution_to_ids](Graph::Node src,
                            vector<Graph::Arc> & cur_path, OwnNode_t anchor_point, float & it_max_flow)->void {
                float max_flow = -INF, max_flow_pointed = -INF;
                Graph::Arc selected, selected_by_pointed;
                /*
                 * Path from the maximal available flow
                 */
                for (Graph::OutArcIt a(solution_g, src); a != INVALID; ++a)
                {
                   
                    if (max_flow < available_flow[a]){
                        max_flow = available_flow[a];
                        selected = a;
                    }
                    if (anchor_point != INF){
                        OwnNode_t val = _g_nodes[solution_to_ids[solution_g.target(a)]]._val;
                        if ((_paired_with[val].find(anchor_point) != _paired_with[val].end()) & (max_flow_pointed < available_flow[a])){
                            max_flow_pointed = available_flow[a];
                            selected_by_pointed = a;
                        }
                    }
                }
                if (max_flow == -INF)
                    return;
                if (max_flow_pointed != -INF)
                {
                    max_flow = max_flow_pointed;
                    selected = selected_by_pointed;
                }
                if (max_flow != 0)
                    it_max_flow = (max_flow < it_max_flow) ? max_flow : it_max_flow;
                if (targets_set.find(solution_to_ids[solution_g.target(selected)]) != targets_set.end())
                {
                    cur_path.push_back(selected);
                    return;
                } 
                Node target_node = _g_nodes[solution_to_ids[solution_g.target(selected)]];
                if (_g_in_edges[target_node._id].size() > 1)
                    anchor_point = _g_nodes[solution_to_ids[solution_g.source(selected)]]._val;
                cur_path.push_back(selected);
                __flow_to_path(solution_g.target(selected),cur_path, anchor_point, it_max_flow);
            };
            /*
             * Flow to paths
             */
            unordered_set<OwnNode_t> sources_added;
            vector<bool> traversed_node(_g_nodes.size(),false);
            size_t total = _g_nodes.size();
            while(true)
            {
                vector <Graph::Arc> cur_path;
                OwnNode_t anchor_point = INF;
                Graph::Arc selected;
                float flow = -INF;
                nw_file << "Getting source with highest flow (swag :P): "<<endl;
                for (auto s: sources_set)
                {
                    Graph::Node source_node = offset_to_solution[s];
                    Graph::OutArcIt a(solution_g, source_node);
                    /*if (a == INVALID && sources_added.find(s) == sources_added.end())
                    {
                        vector <OwnNode_t> real_path(1,s);
                        result_paths.push({0,real_path});
                        sources_added.emplace(s);
                        if (show)
                            nw_file << "Source added: "<<s<<endl;
                    }*/
                    if (show)
                        nw_file << " Source: "<<s<<endl;
                    for (Graph::OutArcIt a(solution_g, source_node); a != INVALID; ++a) {
                        float study_flow = available_flow[a];
                        if (show)
                            nw_file << " -> "<< solution_to_ids[solution_g.target(a)]<<"   Flow:  "<<study_flow<<endl;
                        if (study_flow > flow) {
                            flow = study_flow;
                            selected = a;
                        }
                    }
                }
                if (flow == 0 || flow == -INF)
                    break;
                if (show) {
                    nw_file << "And the winner is... " << flow << " From: " << solution_to_ids[solution_g.source(selected)]<<" "
                         << solution_to_ids[solution_g.target(selected)] << "¿Real? "<<true_nodes[solution_to_ids[solution_g.source(selected)]]<< endl;
                }
                bool store = true_nodes[solution_to_ids[solution_g.source(selected)]];
                if (!store)
                    nw_file << "Not real path - avoiding storing the path. STORE: "<<store<<endl;
                else
                    nw_file << "Path from real node. STORE: "<<store<<endl;
                cur_path.push_back(selected);
                __flow_to_path(solution_g.target(selected), cur_path, anchor_point, flow);
                if (show)
                    nw_file << "Max flow: " << flow << " anchor point: "<<anchor_point<< endl;
                vector <OwnNode_t> real_path;
                float rescore = 0.0;
                for (size_t i_p = 0; i_p < cur_path.size(); i_p++) {
                    OwnNode_t solution_node_id_source = solution_to_ids[solution_g.source(cur_path[i_p])]
                            , solution_node_id_target = solution_to_ids[solution_g.target(cur_path[i_p])];
                    if (!traversed_node[solution_node_id_source])
                    {
                        total--;
                        traversed_node[solution_node_id_source] = true;
                    }
                    if (!traversed_node[solution_node_id_target])
                    {
                        total--;
                        traversed_node[solution_node_id_target] = true;
                    }
                    available_flow[cur_path[i_p]] = (available_flow[cur_path[i_p]] - flow < 0) ? 0 :
                                                    available_flow[cur_path[i_p]] - flow;
                    if (show) {
                        nw_file << "From: " << solution_node_id_source << " to "
                             << solution_node_id_target<<" "<<rescoring[solution_node_id_source]<<" "<<store << endl;
                    }
                    if (store)
                    {
                        rescore += rescoring[solution_node_id_source];
                        real_path.push_back(solution_node_id_source);
                        if (i_p == cur_path.size() - 1) {
                            rescore += rescoring[solution_node_id_target];
                            real_path.push_back(solution_node_id_target);
                        }
                    }
                }
                float sum_rescoring = Common::sum_vector(rescoring);
                float denom = sum_rescoring-1;
                denom = (denom == 0)?1:denom;
                // Change rescore 1 + rescore -> we give extrapoints to more legit paths but we not penalize others
                rescore = 1 + (((real_path.size() - 2) > 0)?(rescore - 1)/denom:0);
                //rescore = ((real_path.size() - 2) > 0)?(rescore - 1)/(real_path.size() - 1):1;
                if (show) {
                    if (RESCORING)
                        nw_file << "New path with " <<flow<<" Rescore: "<<rescore<<" "<<(denom)<<" With Rescore: "
                                <<ceil(flow*rescore)<<" Denom: "<<denom << " flow Path size: "<< real_path.size()<< endl;
                    else
                        nw_file << "New path with " << flow << " flow." << endl;
                }
                flow = (RESCORING)?ceil(flow*rescore):flow;
                if (flow >= MIN_FLOW_PATH)
                {
                    result_paths.push(pair<size_t, vector < OwnNode_t >>(flow, real_path));
                }
                /*if (total == 0)
                {
                    if (show)
                        nw_file << "Set of paths is already a path cover - Finishing path traversal"<<endl;
                    break;
                }*/
            }
            if (show)
                nw_file << "End MCF"<<endl;
            break;
        case NS::UNBOUNDED:
            cerr << "infinite flow" << endl;
            exit(1);
            break;
        default:
            break;
    }
    return result_paths;
}

/*
 * Get maximal flow for compute min_cost_flow
 */
float DBG::to_max_flow_solution()
{
    using namespace lemon;
    typedef ListDigraph Graph;
    typedef int Capacity;

    Graph fn;
    Graph::NodeMap<OwnNode_t> translation_map(fn);
    vector<Graph::Node> id_to_fn = vector<Graph::Node>(_g_nodes.size() + 2);
    Graph::ArcMap<Capacity> available_flow(fn);
    /*
     * Adding nodes
     */
    Graph::Node s = fn.addNode(), t = fn.addNode();
    size_t ids = _g_nodes.size();
    id_to_fn[ids] = s;
    translation_map[s] = ids++;
    id_to_fn[ids] = t;
    translation_map[t] = ids;
    /*
     * Solution network containers
     */
    Graph::ArcMap<int> capacity(fn);
    Graph::Arc a = fn.addArc(t,s);
    capacity[a] = INF;
    for (size_t i = 0; i < _g_nodes.size(); ++i)
    {
        Graph::Node n = fn.addNode();
        translation_map[n] = _g_nodes[i]._id;
        id_to_fn[_g_nodes[i]._id] = n;
    }
    /*
     * Edges from s to sources and from sinks to t with inf capacity
     */
    vector<UG_Node> starting_points = _get_starting_nodes_basic(), sink_nodes = _get_potential_sinks_basic();
    for (auto start_point: starting_points)
    {
        OwnNode_t id = start_point._id;
        Graph::Arc a = fn.addArc(s, id_to_fn[id]);
        capacity[a] = INF;
    }
    for (auto sink_point: sink_nodes)
    {
        OwnNode_t id = sink_point._id;
        Graph::Arc a = fn.addArc(id_to_fn[id], t);
        capacity[a] = INF;
    }
    /*
     * Remaining edges
     */
    size_t node_max = 0, neigh_max = 0, cap = 0;
    for (size_t i = 0; i < _g_edges.size(); ++i)
    {
        for (size_t j = 0; j < _g_edges[i].size(); ++j)
        {
            Graph::Arc a = fn.addArc(id_to_fn[_g_nodes[i]._id], id_to_fn[_g_edges[i][j]]);
            capacity[a] = (int) _g_edges_reads[i][j];
            if (capacity[a] > cap){
                node_max = _g_nodes[i]._val;
                neigh_max = _g_nodes[_g_edges[_g_nodes[i]._id][j]]._val;
                cap = capacity[a];
            }
        }
    }
    cout << "Check: "<<cap << " " <<node_max<<" "<<neigh_max<<endl;
    /*
     * Preflow
     */
    Preflow<Graph> preflow (fn, capacity, s, t);
    Graph::ArcMap<Capacity> flows(fn);
    preflow.run();

    /*
     * Flow to paths
     */
    cout << "Maximum flow by edmonds karp: " << preflow.flowValue() << endl;
    //cin.get();
    return preflow.flowValue();
}

/*
 * Get flow from precomputed max_flow based on min-cost flow (easy approach)
 */
priority_queue<pair<size_t,vector<OwnNode_t>>> DBG::get_min_cost_flow_paths(float pre_flow)
{
    Graph fn;
    Graph::NodeMap<OwnNode_t> translation_map(fn);
    vector<Graph::Node> id_to_fn = vector<Graph::Node>(_g_nodes.size() + 2);
    Graph::ArcMap<Capacity> available_flow(fn);
    /*
     * Adding nodes
     */
    Graph::Node s = fn.addNode(), t = fn.addNode();
    size_t ids = _g_nodes.size();
    id_to_fn[ids] = s;
    translation_map[s] = ids++;
    id_to_fn[ids] = t;
    translation_map[t] = ids;
    /*
     * Solution network containers
     */
    Graph::ArcMap<Capacity> capacity(fn);
    Graph::ArcMap<Weight> weights(fn);
    for (size_t i = 0; i < _g_nodes.size(); ++i)
    {
        Graph::Node n = fn.addNode();
        translation_map[n] = _g_nodes[i]._id;
        id_to_fn[_g_nodes[i]._id] = n;
    }
    /*
     * Edges from s to sources and from sinks to t with inf capacity
     */
    vector<UG_Node> sources_nodes = _get_starting_nodes_basic(), sink_nodes = _get_potential_sinks_basic();
    for (auto start_point: sources_nodes)
    {
        OwnNode_t id = start_point._id;
        Graph::Arc a = fn.addArc(s, id_to_fn[id]);
        capacity[a] = INF;weights[a] = 0.0;
    }
    for (auto sink_point: sink_nodes)
    {
        OwnNode_t id = sink_point._id;
        Graph::Arc a = fn.addArc(id_to_fn[id], t);
        capacity[a] = INF;weights[a] = 0.0;
    }
    /*
     * Rest edges
     */
    for (size_t i = 0; i < _g_edges.size(); ++i)
    {
        /*
         * Forward edge
         */
        for (size_t j = 0; j < _g_edges[i].size(); ++j)
        {
            Graph::Arc a = fn.addArc(id_to_fn[_g_nodes[i]._id], id_to_fn[_g_edges[i][j]]);
            capacity[a] = _g_edges_reads[i][j]; weights[a] = (1/((_g_edges_reads[i][j] == 0)?0.001:_g_edges_reads[i][j]));
        }
    }

    cout << "Distributing the flow: "<<pre_flow<<endl;
    NS ns(fn);
    ns.costMap(weights).upperMap(capacity).stSupply(s, t, pre_flow);

    ArcMap<Capacity> flows(fn);
    NS::ProblemType status = ns.run();
    switch (status) {
        case NS::INFEASIBLE:
            cerr << "insufficient flow" << endl;
            break;
        case NS::OPTIMAL:
            unordered_set<OwnNode_t> targets_set, sources_set;
            priority_queue<pair<size_t,vector<OwnNode_t>>> result_paths;
            for (auto t: sink_nodes)
                targets_set.emplace(t._id);
            for (auto s: sources_nodes)
                sources_set.emplace(s._id);
            std::function<void(Graph::Node, vector<Graph::Arc>&, float&)> __flow_to_path =
                    [this,& __flow_to_path,& fn, &available_flow, &targets_set,&translation_map]
                        (Graph::Node src,vector<Graph::Arc> & cur_path, float & it_max_flow)->void {
                        float max_flow = -INF;
                        Graph::Arc selected;
                        /*
                         * Path from the maximal available flow
                         */
                        for (Graph::OutArcIt a(fn, src); a != INVALID; ++a)
                        {
                            if (max_flow < available_flow[a]){
                                max_flow = available_flow[a];
                                selected = a;
                            }
                        }
                        if (max_flow == -INF)
                            return;
                        if (max_flow != 0)
                            it_max_flow = (max_flow < it_max_flow) ? max_flow : it_max_flow;
                        if (targets_set.find(translation_map[fn.target(selected)]) != targets_set.end())
                        {
                            cur_path.push_back(selected);
                            return;
                        }
                        cur_path.push_back(selected);
                        __flow_to_path(fn.target(selected),cur_path, it_max_flow);
                    };
            ns.flowMap(flows);
            cout << "cost=" << ns.totalCost() << endl;
            cout << "Filling flow: "<<endl;
            for (Graph::ArcIt a(fn); a != INVALID; ++a)
                available_flow[a] = ns.flow(a);
            cout << "Checking flow: "<<endl;
            unordered_set<OwnNode_t> sources_added;
            vector<bool> traversed_node(_g_nodes.size(),false);
            size_t total = _g_nodes.size();
            while(true)
            {
                vector <Graph::Arc> cur_path;
                Graph::Arc selected;
                float flow = -INF;
                cout << "Getting source with highest flow: "<<endl;
                for (auto s: sources_set)
                {
                    Graph::Node source_node = id_to_fn[s];
                    Graph::OutArcIt a(fn, source_node);
                    if (a == INVALID && sources_added.find(s) == sources_added.end())
                    {
                        vector <OwnNode_t> real_path(1,s);
                        result_paths.push({0,real_path});
                        sources_added.emplace(s);
                        cout << "Source added: "<<s<<endl;
                    }
                    for (Graph::OutArcIt a(fn, source_node); a != INVALID; ++a) {
                        float study_flow = available_flow[a];
                        if (study_flow > flow) {
                            flow = study_flow;
                            selected = a;
                        }
                    }
                }
                if (flow == 0 || flow == -INF) {
                    cout << "No more path available!"<<endl;
                    break;
                }
                cout << "And the winner is... " << flow << " From: " << translation_map[fn.source(selected)]<<" "
                            << translation_map[fn.target(selected)]<< endl;
                cur_path.push_back(selected);
                __flow_to_path(fn.target(selected), cur_path, flow);
                cout << "Max flow: " << flow << endl;
                vector <OwnNode_t> real_path;
                for (size_t i_p = 0; i_p < cur_path.size(); i_p++) {
                    OwnNode_t solution_node_id_source = translation_map[fn.source(cur_path[i_p])]
                    , solution_node_id_target = translation_map[fn.target(cur_path[i_p])];
                    if (!traversed_node[solution_node_id_source])
                    {
                        total--;
                        traversed_node[solution_node_id_source] = true;
                    }
                    if (!traversed_node[solution_node_id_target])
                    {
                        total--;
                        traversed_node[solution_node_id_target] = true;
                    }
                    if (available_flow[cur_path[i_p]] - flow < 0)
                    {
                        cout << "Something goes wrong!"<<endl;
                        exit(1);
                    }
                    available_flow[cur_path[i_p]] = (available_flow[cur_path[i_p]] - flow < 0) ? 0 :
                                                    available_flow[cur_path[i_p]] - flow;
                    cout << "From: " << solution_node_id_source << " to "<< solution_node_id_target<< endl;
                    if (solution_node_id_source != ids - 1)
                        real_path.push_back(_g_nodes[solution_node_id_source]._val);
                    if (i_p == cur_path.size() - 1 && solution_node_id_target != ids)
                        real_path.push_back(_g_nodes[solution_node_id_target]._val);
                }
                if (flow >= MIN_FLOW_PATH)
                {
                    result_paths.push(pair<size_t, vector < OwnNode_t >>((flow == INF)?_g_nodes[translation_map[fn.target(cur_path[0])]]._abundance:flow, real_path));
                }
            }
            cout << "Results path: "<<result_paths.size()<<endl;
            return result_paths;
            break;
        case NS::UNBOUNDED:
            cerr << "infinite flow" << endl;
            break;
        default:
            break;
    }
}

/*
 * MCP with no flow in a graph
 * */
priority_queue<pair<size_t,vector<OwnNode_t>>> DBG::solve_std_mcp(const vector<string> & sequence_map)
{
    /*
     * Debugging files
     */
    std::ofstream nw_file;
    if (Parameters::get().debug)
        nw_file.open ("graphs/mcp_apdb_log.txt",std::ofstream::out);
    /*
     * Flow network
     */
    Graph fn;
    Graph::NodeMap<OwnNode_t> translation_map(fn);
    vector<Graph::Node> id_to_fn = vector<Graph::Node>(_g_nodes.size() + 4);
    Graph::ArcMap<Capacity> available_flow(fn);
    /*
     * When no flow is required
     */
    Graph::NodeMap<float> supplyMap(fn);
    Graph::ArcMap<bool> trueArc(fn);
    Graph::ArcMap<Graph::Arc> backward_to_forward(fn), forward_to_backward(fn);
    /*
     * Adding nodes
     */
    Graph::Node s = fn.addNode(), t = fn.addNode();
    Graph::Node s_star = fn.addNode(), t_star = fn.addNode();
    size_t ids = _g_nodes.size();
    id_to_fn[ids] = s;translation_map[s] = ids++;
    id_to_fn[ids] = t;translation_map[t] = ids++;
    id_to_fn[ids] = s_star;translation_map[s_star] = ids++;
    id_to_fn[ids] = t_star;translation_map[t_star] = ids;
    nw_file << "Source global "<<(ids - 3)<<" Target global "<<(ids - 2)<<" Source star "<<(ids - 1)<<" Target star "<<(ids)<<endl;
    /*
     * Solution network containers
     */
    Graph::ArcMap<Capacity> capacity(fn), lowerMap(fn);
    Graph::ArcMap<Weight> weights(fn);
    for (size_t i = 0; i < _g_nodes.size(); ++i)
    {
        Graph::Node n = fn.addNode();
        supplyMap[n] = 0;
        translation_map[n] = _g_nodes[i]._id;
        id_to_fn[_g_nodes[i]._id] = n;
    }
    supplyMap[s] = 0;
    supplyMap[t] = 0;
    /*
     * Edges from s to sources and from sinks to t with inf capacity
     */
    vector<UG_Node> sources_nodes = _get_starting_nodes_basic(), sink_nodes = _get_potential_sinks_basic();
    vector<OwnNode_t> sources_id_nodes, sink_id_nodes;
    for (auto start_point: sources_nodes)
    {
        OwnNode_t id = start_point._id;
        sources_id_nodes.emplace_back(id);
        if (Parameters::get().debug)
            nw_file << "Source: "<<id<<endl;
        Graph::Arc a = fn.addArc(s, id_to_fn[id]);
        capacity[a] = INF;weights[a] = 0.0;lowerMap[a] = 0;
        trueArc[a] = false;
    }
    for (auto sink_point: sink_nodes)
    {
        OwnNode_t id = sink_point._id;
        sink_id_nodes.emplace_back(id);
        if (Parameters::get().debug)
            nw_file << "Sink: "<<id<<endl;
        Graph::Arc a = fn.addArc(id_to_fn[id], t);
        capacity[a] = INF;weights[a] = 0.0;lowerMap[a] = 0;
        trueArc[a] = true;
    }
    /*
     * Rest edges
     */
    for (size_t i = 0; i < _g_edges.size(); ++i)
    {
        vector<Graph::Arc> arcs;
        nw_file << "Node: "<<i<<endl;
        float length_unitig = (float) Common::return_unitig(sequence_map, _g_nodes[i]._val).size();
        /*
         * Forward edge
         */
        for (size_t j = 0; j < _g_edges[i].size(); ++j)
        {
            Graph::Node src = id_to_fn[_g_nodes[i]._id], tg = id_to_fn[_g_edges[i][j]];
            Graph::Arc a = fn.addArc(src, tg);
            trueArc[a] = true;arcs.push_back(a);lowerMap[a] = 0;
            capacity[a] = _g_edges_reads[i][j]; 
            //weights[a] = (1/((_g_edges_reads[i][j] == 0)?0.001:_g_edges_reads[i][j]));
            float weight = 0;
            // Check extreme cases.
            if (_g_edges[i].size() > 1 || in_degree(_g_edges[i][j]) > 1 || in_degree(i) > 1 || in_degree(i) == 0 )
                weight = 1/_g_edges_reads[i][j];
            // Lets prioritize removing some flow in the end rather than increase on the rest of the strain
            if (out_degree(_g_edges[i][j]) == 0)
                weight = 2/_g_edges_reads[i][j];
            weights[a] = weight;
            nw_file << "neigh: "<<_g_edges[i][j]<<" Flow: "<<capacity[a]<<" Weight: "<<weights[a]<<endl;
            /*
             * Supply map
             */
            supplyMap[src] += capacity[a];
            supplyMap[tg] -= capacity[a];
            nw_file << "Supply: "<<i<<" "<<supplyMap[src]<<endl;
            nw_file << "Supply: "<<_g_edges[i][j]<<" "<<supplyMap[tg]<<endl;
        }
        /*
         * Backward edge
         */
        for (size_t j = 0; j < _g_edges[i].size(); ++j)
        {
            Graph::Node tg = id_to_fn[_g_nodes[i]._id], src = id_to_fn[_g_edges[i][j]];
            Graph::Arc a = fn.addArc(src, tg);
            trueArc[a] = false;backward_to_forward[a] = arcs[j];lowerMap[a] = 0;
            forward_to_backward[arcs[j]] = a;
            capacity[a] = _g_edges_reads[i][j]; 
            //weights[a] = (1/((_g_edges_reads[i][j] == 0)?0.001:_g_edges_reads[i][j]));
            float weight = 0;
            if (_g_edges[i].size() > 1 || in_degree(_g_edges[i][j]) > 1 || in_degree(i) > 1 || in_degree(i) == 0 )
                weight = 1/_g_edges_reads[i][j];
            if (out_degree(_g_edges[i][j]) == 0)
                weight = 2/_g_edges_reads[i][j];
            weights[a] = weight;
        }
    }
    /*
     * Exogenous flow
     */
    for (size_t i = 0; i < _g_nodes.size(); ++i)
    {
        Graph::Node src = id_to_fn[_g_nodes[i]._id];
        if (supplyMap[src] > 0)
        {
            Graph::Arc a = fn.addArc(src, t_star);
            capacity[a] = supplyMap[src];weights[a] = 0;
            supplyMap[src] = 0;
            trueArc[a] = false;lowerMap[a] = 0;
            supplyMap[t_star] -= capacity[a];
            if (Parameters::get().debug)
                nw_file << "Edge from: "<<_g_nodes[i]._id<<" Val: "<<_g_nodes[i]._val<<" to target_star with flow "<<capacity[a]
                << " Is source? "<<(find(sources_id_nodes.begin(), sources_id_nodes.end(),_g_nodes[i]._id) != sources_id_nodes.end())<<endl;
        } else if (supplyMap[src] < 0)
        {
            Graph::Arc a = fn.addArc(s_star, src);
            capacity[a] = -supplyMap[src]; weights[a] = 0;
            supplyMap[src] = 0;
            trueArc[a] = false;lowerMap[a] = 0;
            supplyMap[s_star] += capacity[a];
            if (Parameters::get().debug)
                nw_file << "Edge from source_star to: "<<_g_nodes[i]._id<<" Val: "<<_g_nodes[i]._val<<" with flow "<<capacity[a]<<
                " Is target? "<<(find(sink_id_nodes.begin(),sink_id_nodes.end(),_g_nodes[i]._id) != sink_id_nodes.end())<<endl;
        }
    }
    if (Parameters::get().debug)
        nw_file << "Source star flow: "<<supplyMap[s_star]<<" Sink star flow: "<<supplyMap[t_star]<<endl;
    /*
     * Arc from sink to source
     */
    Graph::Arc a = fn.addArc(t, s);
    capacity[a] = INF; weights[a] = 0;
    trueArc[a] = false;lowerMap[a] = 0;
    /*
     * Build the network flow
     */
    NS ns(fn);
    ns.costMap(weights).lowerMap(lowerMap).upperMap(capacity).supplyMap(supplyMap);

    ArcMap<Capacity> flows(fn);
    NS::ProblemType status = ns.run(NS::BLOCK_SEARCH);
    switch (status) {
        case NS::INFEASIBLE:
            cerr << "insufficient flow" << endl;
            break;
        case NS::OPTIMAL:
            unordered_set<OwnNode_t> targets_set, sources_set;
            priority_queue<pair<size_t,vector<OwnNode_t>>> result_paths;
            for (auto t: sink_nodes)
                targets_set.emplace(t._id);
            for (auto s: sources_nodes)
                sources_set.emplace(s._id);
            std::function<void(Graph::Node, vector<Graph::Arc>&, float&)> __flow_to_path =
                    [this,& __flow_to_path,& fn, &available_flow,
                     &targets_set,&translation_map, &trueArc]
                            (Graph::Node src,vector<Graph::Arc> & cur_path, float & it_max_flow)->void {
                        float max_flow = -INF, min_diff = INF;
                        Graph::Arc selected;
                        /*
                         * Neigh selection based on the minimum distance to the current it_max_flow
                         * Or max_flow still deciding
                         */
                        for (Graph::OutArcIt a(fn, src); a != INVALID; ++a)
                        {
                            if (!trueArc[a])
                                continue;
                            /*if (min_diff > abs(available_flow[a] - it_max_flow) & (available_flow[a] != 0)){
                                max_flow = available_flow[a];
                                min_diff = abs(available_flow[a] - it_max_flow);
                                selected = a;
                            }*/
                            if (available_flow[a] > max_flow)
                            {
                                max_flow = available_flow[a];
                                selected = a;
                            }
                        }
                        if (max_flow == -INF)
                            return;
                        if (max_flow != 0)
                            it_max_flow = (max_flow < it_max_flow) ? max_flow : it_max_flow;
                        if (targets_set.find(translation_map[fn.target(selected)]) != targets_set.end())
                        {
                            cur_path.push_back(selected);
                            return;
                        }
                        cur_path.push_back(selected);
                        __flow_to_path(fn.target(selected),cur_path, it_max_flow);
                    };
            ns.flowMap(flows);
            nw_file << "cost=" << ns.totalCost() << endl;
            nw_file << "Filling flow: "<<endl;
            for (Graph::ArcIt a(fn); a != INVALID; ++a)
            {
                    available_flow[a] = (trueArc[a])?capacity[a]:0;
            }
            /*
             * Just for testing purposes
             */
            nw_file << "Checking flow conservative rule"<<endl;
            for (Graph::NodeIt i(fn); i != INVALID; ++i)
            {
                nw_file << "Node: "<<translation_map[i]<<endl;
                nw_file << "Required flow: "<<supplyMap[i]<<endl;
                for (Graph::OutArcIt edge(fn, i); edge != INVALID; ++edge)
                {
                    nw_file << "Outflow node: "<<flows[edge]<<endl;
                    nw_file << "Capacity: "<<capacity[edge]<<endl;
                    nw_file << "Source: "<<translation_map[fn.source(edge)]<<" Target: "<<translation_map[fn.target(edge)]<<endl;
                }
                for (Graph::InArcIt edge(fn, i); edge != INVALID; ++edge)
                {
                    nw_file << "Inflow node: "<<flows[edge]<<endl;
                    nw_file << "Capacity: "<<capacity[edge] <<endl;
                    nw_file << "Source: "<<translation_map[fn.source(edge)]<<" Target: "<<translation_map[fn.target(edge)]<<endl;
                }
            }
            for (Graph::ArcIt a(fn); a != INVALID; ++a)
            {
                if (trueArc[a]){
                    available_flow[a] += flows[a];
                    if (Parameters::get().debug) {
                            nw_file << "TRUE ARC"<<endl;
                            nw_file << "Source: " << translation_map[fn.source(a)] << " Target: "
                                    << translation_map[fn.target(a)] << " Capacity: " <<capacity[a] 
                                    << " Available flow: " << available_flow[a] << " Flow (forward): "
                                    << flows[a]<<endl;
                            nw_file << "Flow (backward): "<<flows[forward_to_backward[a]]<<" Source: "<<translation_map[fn.source(forward_to_backward[a])]
                                    << " Target: "<<translation_map[fn.target(forward_to_backward[a])]<<endl;
                    }
                } else if (fn.target(a) != t_star && fn.source(a) != s_star && fn.source(a) != s) {
                    available_flow[backward_to_forward[a]] -= flows[a];
                    /*if (Parameters::get().debug) {
                        if (!trueArc[a]) {
                            nw_file << "Source: " << translation_map[fn.source(a)] << " Target: "
                                    << translation_map[fn.target(a)]
                                    << " Available flow: " << available_flow[backward_to_forward[a]] << " Flows: "
                                    << flows[a]
                                    << " Flows(forward): " << flows[backward_to_forward[a]] << endl;
                            nw_file << "Source: " << translation_map[fn.source(backward_to_forward[a])]
                                    << " Target: " << translation_map[fn.target(backward_to_forward[a])] << endl;
                        }
                    }*/
                } else {
                    /*if (Parameters::get().debug) {
                        if (!trueArc[a]) {
                            nw_file << "From exceptional point: "<<endl;
                            nw_file << "Source: " << translation_map[fn.source(a)] << " Target: "
                                    << translation_map[fn.target(a)]
                                    << " Available flow: " << available_flow[backward_to_forward[a]] << " Flows: "
                                    << flows[a]
                                    << " Flows(forward): " << flows[backward_to_forward[a]] << endl;
                            nw_file << "Source: " << translation_map[fn.source(backward_to_forward[a])]
                                    << " Target: " << translation_map[fn.target(backward_to_forward[a])] << endl;
                        }
                    }*/
                }
            }
            unordered_set<OwnNode_t> sources_added;
            vector<bool> traversed_node(_g_nodes.size(),false);
            size_t total = _g_nodes.size();
            while(true)
            {
                vector <Graph::Arc> cur_path;
                Graph::Arc selected;
                float flow = -INF;
                nw_file << "Getting source with highest flow: "<<endl;
                for (auto s: sources_set)
                {
                    Graph::Node source_node = id_to_fn[s];
                    Graph::OutArcIt a(fn, source_node);
                    if (a == INVALID && sources_added.find(s) == sources_added.end())
                    {
                        vector <OwnNode_t> real_path(1,s);
                        result_paths.push({0,real_path});
                        sources_added.emplace(s);
                        nw_file << "Source added: "<<s<<endl;
                    }
                    for (Graph::OutArcIt a(fn, source_node); a != INVALID; ++a) {
                        if (!trueArc[a])
                            continue;
                        float study_flow = available_flow[a];
                        if (study_flow > flow) {
                            flow = study_flow;
                            selected = a;
                        }
                    }
                }
                if (flow == 0 || flow == -INF) {
                    nw_file << "No more path available!"<<endl;
                    break;
                }
                nw_file << "And the winner is... " << flow << " From: " << translation_map[fn.source(selected)]<<" "
                     << translation_map[fn.target(selected)]<< endl;
                cur_path.push_back(selected);
                __flow_to_path(fn.target(selected), cur_path, flow);
                nw_file << "Max flow: " << flow << endl;
                vector <OwnNode_t> real_path;
                for (size_t i_p = 0; i_p < cur_path.size(); i_p++) {
                    OwnNode_t solution_node_id_source = translation_map[fn.source(cur_path[i_p])]
                    , solution_node_id_target = translation_map[fn.target(cur_path[i_p])];
                    nw_file << "From: " << solution_node_id_source << " to "<< solution_node_id_target<<" Cur flow: "<<available_flow[cur_path[i_p]] << endl;
                    if (!traversed_node[solution_node_id_source])
                    {
                        total--;
                        traversed_node[solution_node_id_source] = true;
                    }
                    if (!traversed_node[solution_node_id_target])
                    {
                        total--;
                        traversed_node[solution_node_id_target] = true;
                    }
                    if (available_flow[cur_path[i_p]] - flow < 0)
                    {
                        nw_file << "Something goes wrong!" << available_flow[cur_path[i_p]]<<" "<<flow<<endl;
                        exit(1);
                    }
                    available_flow[cur_path[i_p]] = (available_flow[cur_path[i_p]] - flow < 0) ? 0 :
                                                    available_flow[cur_path[i_p]] - flow;
                    if (solution_node_id_source != (ids - 3))
                        real_path.push_back(_g_nodes[solution_node_id_source]._val);
                    if (i_p == cur_path.size() - 1 && solution_node_id_target != (ids-2))
                        real_path.push_back(_g_nodes[solution_node_id_target]._val);
                }
                if (flow >= MIN_FLOW_PATH)
                {
                    result_paths.push(pair<size_t, vector < OwnNode_t >>((flow == INF)?_g_nodes[translation_map[fn.source(cur_path[0])]]._abundance:flow, real_path));
                }
            }
            nw_file << "Results path: "<<result_paths.size()<<endl;
            return result_paths;
            break;
        case NS::UNBOUNDED:
            cerr << "infinite flow" << endl;
            break;
        default:
            break;
    }
}

/*
 * MCP Example
 */
void DBG::mcp_example()
{
    using namespace lemon;
    using namespace std;


    using Weight = int;
    using Capacity = int;

    using Graph = SmartDigraph;

    using Node = Graph::Node;
    using Arc = Graph::Arc;
    using ArcMap = SmartDigraph::ArcMap<int>;

    using NS = NetworkSimplex<SmartDigraph, Capacity, Weight>;

    Graph g;
    /*
     * u - 5 -> v - 3 -> x
     */
    Node u = g.addNode();
    Node v = g.addNode();
    Node x = g.addNode();
    Node s = g.addNode(), s_s = g.addNode(), t = g.addNode(), t_s = g.addNode();
    Arc a = g.addArc(u, v), a_r = g.addArc(v,u);
    Arc b = g.addArc(v,x), b_r = g.addArc(x,u);
    Arc uts = g.addArc(u,t_s), ssv = g.addArc(s_s,v), ssb = g.addArc(s_s,x);
    Arc su = g.addArc(s,u), xt = g.addArc(x,t);
    Arc ts = g.addArc(t,s);

    ArcMap weights(g);
    weights[a] = 5;weights[a_r] = 5;
    weights[b] = 3;weights[b_r] = 3;
    weights[ssb] = 0;weights[uts] = 0; weights[ssv] = 0;
    weights[su] = 0;weights[xt] = 0;
    weights[ts] = 0;

    ArcMap capacities(g);
    capacities[a] = 5;capacities[a_r] = 5;
    capacities[b] = 3;capacities[b_r] = 3;
    capacities[ssb] = 3;capacities[uts] = 5; capacities[ssv] = 2;
    capacities[su] = INF;capacities[xt] = INF;
    capacities[ts] = INF;

    SmartDigraph::NodeMap<int> supply(g);
    supply[s_s] = 5;supply[t_s] = -5;

    NS ns(g);
    ns.costMap(weights).upperMap(capacities).supplyMap(supply);

    ArcMap flows(g);
    NS::ProblemType status = ns.run();
    switch (status) {
    case NS::INFEASIBLE:
        cout << "insufficient flow" << endl;
        break;
    case NS::OPTIMAL:
        ns.flowMap(flows);
        cout << "flow[a]=" << ns.flow(a) << " flow[b]=" << flows[b] << endl;
        cout << "flow[ar]=" << ns.flow(a_r) << " flow[b_r]=" << flows[b_r] << endl;
        cout << "flow[ssb]=" << ns.flow(ssb) << " flow[uts]=" << flows[uts] << endl;
        cout << "flow[xt]=" << ns.flow(xt) << endl;
        cout << "flow[ts]=" << ns.flow(ts) << endl;
        cout << "cost=" << ns.totalCost() << endl; 
        break;
    case NS::UNBOUNDED:
        cout << "infinite flow" << endl;
        break;
    default:
        break;
    }

    return 0;
}