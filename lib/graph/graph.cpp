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
                for (size_t i = 0; i < results.size(); ++i)
                {
                    if (i == results.size() - 1)
                    {
                        float abundance;
                        istringstream iss(results[i]);
                        iss >> abundance;
                        _g_nodes[node].set_abundance(abundance);
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
                    }
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
    size_t empty_pairs = 0;
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
    for (size_t i = 0; i < _num_nodes; ++i)
    {
        size_t val_node = _g_nodes[i]._val;
        vector<size_t> local_histogram = vector<size_t>(10000);
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
            if (getAbundance(pn) < _quantile_L || getAbundance(pn) > _quantile_H)
                continue;
            node_pei.emplace(pn);
        }
        UG_Node parent(_g_nodes[i]._val, node_pei);
        vector<OwnNode_t> neighbors = getNeighbor(i);
        if (neighbors.size() == 0)
        {
            apdbg.addNode(parent, show);
        }
        //if (show)
            cout << " Node: "<<i<<" "<<_g_nodes[i]._val<<" Neighbors: "<<neighbors.size()<<" InNeighbors: "<<in_degree(i)<<endl;
        if (node_pei.empty())
            empty_pairs++;
        for (size_t n_i = 0; n_i < neighbors.size(); ++n_i)
        {
            auto n = neighbors[n_i];
            auto in_parents = getNeighbor(_g_nodes[n]._id, false);
            size_t val_neigh = _g_nodes[n]._val;
            size_t strain_freq = Maths::my_round(_g_edges_reads[i][n_i], true);
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
            Pairedendinformation_t neigh_pei;
            /*
             * Polish node pairs based on abundance
             */
            for (auto pn: getPairedEndInformation(n))
            {
                if (getAbundance(pn) < _quantile_L || getAbundance(pn) > _quantile_H)
                {
                    rejected_pairs << "Kmer: "<<_g_nodes[pn]._val<<" "<<((size_t)_g_nodes[pn]._abundance)<<" "
                                   <<Common::return_unitig(sequence_map, _g_nodes[pn]._val)<<endl;
                    continue;
                }
                neigh_pei.emplace(pn);
            }
            Pairedendinformation_t neigh_pei_copy = neigh_pei, new_node_pei = node_pei;

            bool content = (neigh_pei.size() > node_pei.size())?is_subset(node_pei,neigh_pei):is_subset(neigh_pei,node_pei);
            if (Parameters::get().debug) {
                basic_information << "Node: "<<_g_nodes[i]._val<<" ";
                basic_information << Common::return_unitig(sequence_map, _g_nodes[i]._val);
                basic_information<<" "<<out_degree(_g_nodes[i]._id)<<endl;
                basic_information << "Size Node (paired-end): " << node_pei.size() << " Neighbors: " << getNeighbor(i).size() << " Indegree: "<<in_degree(_g_nodes[n]._id)<<endl;
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
            * Weird/Strange behaviour?
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
            for (auto pn: union_pei)
            {
                local_histogram[std::min((size_t) 999,getAbundance(n,pn) + getAbundance(i,pn))]++;
                size_t id = local_graph.addVertex(pn);
                float node_abundance = (float) Maths::my_round(getAbundance(pn), true);
                if (Parameters::get().debug)
                    basic_information << "Pair: "<<pn<<" "<<getAbundance(pn)<<" "<<abs(node_abundance - strain_freq)<<" "<<node_abundance<<endl;
                if (READJUST) {
                    node_abundance = abs(node_abundance - strain_freq);
                    max_flow = (max_flow < node_abundance) ? node_abundance : max_flow;
                    min_flow = (min_flow > node_abundance) ? node_abundance : min_flow;
                }
                local_graph_2.addVertexWithAbundances(pn, node_abundance);
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
             * DAG builder
             */
            std::function<void(size_t, size_t,size_t&, size_t&,
                               unordered_set<OwnNode_t> &,
                                       vector<size_t>&, float &, float &, vector<pair<OwnNode_t,size_t>>&)> __reach = [this, &__reach,
                    &traslation_map,
                    &local_graph_2,
                    &union_pei,&show,
                    &strain_freq,
                    &cautious_execution,&val_node,&val_neigh,
                    &basic_information,
                    &true_nodes,
                    &minimum_reached_by,&reached_by,&reverse_map,&frequencies_map](size_t  src, size_t target,
                                 size_t steps, size_t &branches
                    , unordered_set<OwnNode_t> & checked
                    ,vector<size_t> & median_flow,
                                 float & max_flow,
                                 float & min_flow,
                                 vector<pair<OwnNode_t, size_t>> & in_nodes)->void{
                if (Parameters::get().debug)
                    basic_information <<"("<<steps<<","<<in_nodes.size()<<")"<<target<<"("<<src<<","<< median_flow.size()<<")->";
                bool swap = false;
                size_t new_steps = 0, new_branches = 0;
                if (src != target && (in(union_pei, target))) {
                    checked.emplace(target);
                    sort (median_flow.begin(), median_flow.end());
                    //float flow_expected = Maths::my_round(Maths::median(median_flow), true);
                    float flow_expected = Maths::my_round(median_flow[0], true);
                    if (READJUST) {
                        flow_expected = abs(strain_freq - flow_expected);
                        max_flow = (max_flow < flow_expected) ? flow_expected : max_flow;
                        min_flow = (min_flow > flow_expected) ? flow_expected : min_flow;
                    }
                    local_graph_2.addEdge(traslation_map[src], traslation_map[target], flow_expected);
                    for (auto fake_node:in_nodes) {
                        float flow_expected = Maths::my_round(fake_node.second, true);
                        float node_abundance = (float) Maths::my_round(getAbundance(fake_node.first), true);
                        if (READJUST) {
                            node_abundance = abs(node_abundance - strain_freq);
                            flow_expected = abs(strain_freq - flow_expected);
                            float max_comparator = (flow_expected > node_abundance)?flow_expected:node_abundance
                                , min_comparator = (flow_expected < node_abundance)?flow_expected:node_abundance;
                            max_flow = (max_flow < max_comparator) ? max_comparator : max_flow;
                            min_flow = (min_flow > min_comparator) ? min_comparator : min_flow;
                        }
                        if (traslation_map.find(fake_node.first) == traslation_map.end()) {
                            if (Parameters::get().debug)
                                basic_information << "{extra parent "<<fake_node.first<<" son "<<target<<" "<<flow_expected<<" "<<fake_node.second<<"}";
                            size_t id = local_graph_2.addVertexWithAbundances(fake_node.first,node_abundance);
                            reached_by.push_back(minimum_reached_by);
                            true_nodes[id] = false;
                            traslation_map[fake_node.first] = id;
                            frequencies_map[id] = node_abundance;
                            reverse_map[id] = fake_node.first;
                            local_graph_2.addEdge(id,traslation_map[target], flow_expected);
                        } else {
                            if (Parameters::get().debug)
                                basic_information << "{extra parent "<<fake_node.first<<" son "<<target<<" "<<flow_expected<<" "<<fake_node.second<<"}";
                            local_graph_2.addEdge(traslation_map[fake_node.first], traslation_map[target],
                                                  flow_expected);
                        }
                    }
                    src = target;
                    swap = true;
                }
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
                for (size_t i = 0; i < neighbors.size(); ++i) {
                    OwnNode_t neigh = neighbors[i];
                    if (checked.find(neigh) == checked.end())
                    {
                        if (!swap) {
                            /*
                             * Once you swap the current strain you have to reset the fake_nodes
                             */
                            vector<pair<OwnNode_t,size_t>> new_in_nodes = in_nodes;
                            vector<size_t> new_median_flow = median_flow;
                            vector<OwnNode_t> parents = getNeighbor(neigh, false, target);
                            if (neighbors.size() > 1 || parents.size() > 1)
                            {
                                new_in_nodes.clear();
                                for (auto parent:parents)
                                    new_in_nodes.push_back({parent, _get_edge_frec(parent, neigh)});
                            }
                            new_median_flow.push_back(_g_edges_reads[target][i]);
                            size_t local_steps = steps;
                            __reach(src, neigh, local_steps, branches, checked, new_median_flow, max_flow, min_flow, new_in_nodes);
                        } else {
                            if (Parameters::get().debug)
                                basic_information << endl;
                            vector<pair<OwnNode_t,size_t>> new_in_nodes;
                            vector<size_t> new_median_flow(1, _g_edges_reads[target][i]);
                            __reach(src, neigh, new_steps, new_branches, checked, new_median_flow, max_flow, min_flow, new_in_nodes);
                        }
                    } else {
                        if (Parameters::get().debug)
                            basic_information<<"Node of interest already extended, joining both "<<src<<"->"<<neigh<<endl;
                        sort (median_flow.begin(), median_flow.end());
                        //float flow_expected = Maths::my_round(Maths::median(median_flow), true);
                        float flow_expected = Maths::my_round(median_flow[0], true);
                        if (READJUST) {
                            flow_expected = abs(strain_freq - flow_expected);
                            max_flow = (max_flow < flow_expected) ? flow_expected : max_flow;
                            min_flow = (min_flow > flow_expected) ? flow_expected : min_flow;
                        }
                        local_graph_2.addEdge(traslation_map[src], traslation_map[neigh], flow_expected);
                        for (auto fake_node:in_nodes) {
                            float flow_expected = Maths::my_round(fake_node.second, true);
                            float node_abundance = (float) Maths::my_round(getAbundance(fake_node.first), true);
                            if (READJUST) {
                                node_abundance = abs(node_abundance - strain_freq);
                                flow_expected = abs(strain_freq - flow_expected);
                                float max_comparator = (flow_expected > node_abundance)?flow_expected:node_abundance
                                    , min_comparator = (flow_expected < node_abundance)?flow_expected:node_abundance;
                                max_flow = (max_flow < max_comparator) ? max_comparator : max_flow;
                                min_flow = (min_flow > min_comparator) ? min_comparator : min_flow;
                            }
                            if (traslation_map.find(fake_node.first) == traslation_map.end()) {
                                if (Parameters::get().debug)
                                    basic_information << "{extra parent "<<fake_node.first<<" son "<<neigh<<" "<<flow_expected<<"}";
                                size_t id = local_graph_2.addVertexWithAbundances(fake_node.first,
                                                                                  node_abundance);
                                reached_by.push_back(minimum_reached_by);
                                true_nodes[id] = false;
                                traslation_map[fake_node.first] = id;
                                frequencies_map[id] = node_abundance;
                                reverse_map[id] = fake_node.first;
                                local_graph_2.addEdge(id,traslation_map[neigh], flow_expected);
                            } else {
                                if (Parameters::get().debug)
                                    basic_information << "{extra parent "<<fake_node.first<<" son "<<neigh<<" "<<flow_expected<<"}";
                                local_graph_2.addEdge(traslation_map[fake_node.first], traslation_map[neigh],
                                                      flow_expected);
                            }
                        }
                    }
                }
                if (Parameters::get().debug)
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
                unordered_set <OwnNode_t> checked;
                vector<size_t> median_flow;
                vector<pair<OwnNode_t, size_t>> in_nodes;
                if (!cautious_execution) {
                    if (reached_by[traslation_map[pn]] == minimum_reached_by) {
                        checked.emplace(pn);
                        if (Parameters::get().debug)
                            basic_information << "NEW LAUNCH" << endl;
                        __reach(pn, pn, steps, branches, checked, median_flow, max_flow, min_flow, in_nodes);
                    }
                } else if (parent_in_union )
                {
                    if (pn == val_node) {
                        if (Parameters::get().debug)
                            basic_information << "Launching source (parent): " << pn << endl;
                        __reach(pn, pn, steps, branches, checked, median_flow, max_flow, min_flow,in_nodes);
                    }
                } else if (son_in_union) {
                    if (pn == val_neigh)
                    {
                        if (Parameters::get().debug)
                            basic_information << "Launching source (neigh): " << pn << endl;
                        __reach(pn, pn, steps, branches, checked, median_flow, max_flow, min_flow, in_nodes);
                    }
                }
            }
            /*
             * Edge and strain freq flow readjust
             */
            if (READJUST)
            {
                if (Parameters::get().debug)
                {
                    string graph_file = "debug/"+to_string(_g_nodes[i]._val)+"_"+to_string(_g_nodes[n]._val)+"_graph_prereadjust.txt";
                    local_graph_2.export_to_graph(sequence_map, graph_file);
                    basic_information << "Readjusting nodes/edges capacities"<<endl;
                }
                local_graph_2.readjust_flow(min_flow, max_flow);
                strain_freq_check = ceil((MAX_GRANULARITY_ALLOWED*100)*abs((float) max_flow)/max_flow + MAX_GRANULARITY_ALLOWED);
            }
            if (Parameters::get().debug)
                basic_information << "Maximum expected flow: "<<max_flow
                        <<" Minimum expected flow: "<<min_flow<<" Expected flow: "<<strain_freq_check<<endl;
            if (CORRECT_GRAPH)
            {
                string adjust_file = "debug/"+to_string(_g_nodes[i]._val)+"_"+to_string(_g_nodes[n]._val)+"_adjust_flow.txt";
                std::ofstream adjust_information;
                if (Parameters::get().debug)
                {
                    adjust_information.open(adjust_file, std::ofstream::out);
                }
                local_graph_2.correct_graph(adjust_information);
            }
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
                                               (reached_by, _g_edges_reads[i][n_i]
                                                       , Parameters::get().debug, true_nodes,"debug/"+to_string(_g_nodes[i]._val)
                                                       +"_"+to_string(_g_nodes[n]._val)+"_min_cost_flow.txt");
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
                if (top_click.size() <  SIZE_RELATION*max_size)
                    continue;
                if (show)
                    paths << "New Clique:" << endl;
                Pairedendinformation_t clique, p_info_a, p_info_b;
                bool in_a = false, in_b = false, isa = true, isb = true, parent_in = false, son_in = false;
                for (auto node: top_click)
                {
                    nodes_checked.emplace(node);
                    //UG::UG_Node real_node = local_graph[node];
                    UG::UG_Node real_node = local_graph_2[node];
                    clique.emplace(real_node._val);
                    if (Parameters::get().debug) {
                        paths << node<<":"<<real_node._val<<":"<<_g_nodes[real_node._val]._abundance
                            <<":"<<real_node._id<<" "<<Common::return_unitig(sequence_map, real_node._val)<<" ";
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
                if (Parameters::get().debug) {
                    paths << endl;
                    paths << "ISA: "<<isa<<" ISB: "<<isb<<endl;
                    paths << "INA: "<<in_a<<" INB: "<<in_b<<endl;
                    paths << "Is Parent: "<<parent_in<<" Is son: "<<son_in<<endl;
                }
                if (parent_in && son_in) {
                    if (Parameters::get().debug) {
                        paths << "Path forced to be added: Contains parent and son" << endl;
                    }
                    offset_erase.push_back(false);
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
                potential_nodes.push_back({top_click_pair.first,pair<UG_Node, UG_Node>(parent_l, son_l)});
                if (Parameters::get().debug) {
                    paths << "Nodes added as potential ones."<<endl;
                    paths << parent_l.to_String();
                    paths << son_l.to_String();
                }
            }
            /*
             * Procesar las paired-end y eliminar los subconjuntos
             */
            for (size_t index_s = 0; index_s < potential_nodes.size(); ++index_s)
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
            }
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
            if (potential_nodes.size() == 0)
            {
                if (Parameters::get().debug)
                {
                    paths << "No potential nodes, added empty version of both (it only happens when none (parent and son) have paired-end information)"<<endl;
                }
                Pairedendinformation_t empty;
                UG_Node parent_l(_g_nodes[i]._val, empty), son_l(_g_nodes[n]._val,empty);
                potential_nodes.push_back({strain_freq_check,pair<UG_Node,UG_Node>(parent_l,son_l)});
            }
            if (Parameters::get().debug)
                paths << endl<<"Potential nodes: "<<endl;
            for (auto s_pair: potential_nodes) {
                pair<UG_Node, UG_Node> s = s_pair.second;
                if (Parameters::get().debug)
                {
                    paths << s.first.to_String();
                    paths << s.second.to_String();
                }
                pair<OwnNode_t, OwnNode_t> index_pairs = apdbg.addNode(s.first, s.second, show);
                if (Parameters::get().debug)
                    paths << "Adding edge from: "<<index_pairs.first<<" "<<index_pairs.second<<endl;
                /*
                 * Expected flow reached
                 */
                strain_freq_check -= s_pair.first;
                if (strain_freq_check < MIN_FLOW_PATH)
                {
                    strain_freq_check += s_pair.first;
                    break;
                }
            }

            if (Parameters::get().debug) {
                paths << "Finally cliques to add: " << potential_nodes.size() <<" Remaining flow: "<<strain_freq_check<<endl;
            }
        }
        /*
         * Crear directorios antes!
         */
        //string histogram_output = "stats/histograms_nodes/histogram_node"+to_string(i)+".txt";
        //write_histogram(local_histogram, histogram_output);
    }
    if (Parameters::get().debug)
    {
        rejected_nodes.close();
    }
    cout << "Number of empty pairs: "<<empty_pairs<<" out of "<<_g_nodes.size()<<endl;
    //exit(1);
}

size_t DBG::out_degree(OwnNode_t node_id)
{
    return _g_edges[node_id].size();
}

size_t DBG::in_degree(OwnNode_t node_id)
{
    return _g_in_edges[node_id].size();
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
    _get_internal_stats();
    for (size_t i = 0; i < _g_nodes.size(); ++i)
        _g_nodes[i].post_process_pairs(sequence_map);
}

void DBG::_get_internal_stats()
{
    vector<size_t> abundances_v;
    for (auto v:_g_nodes)
    {
        abundances_v.push_back(v._abundance);
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
        if (v._abundance >= _quantile_L) {
            abundances_v.push_back(v._abundance);
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
        if (v._abundance <= _quantile_H && v._abundance >= _quantile_L) {
            abundances_v_post_processed.push_back(v._abundance);
            mean += v._abundance;
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
 * REVISAR
 */
void DBG::add_read(OwnNode_t u1, size_t quantity)
{
    _g_nodes[u1].increase_abundance(quantity);
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
    vector<vector<OwnNode_t>> vector_unitigs(_g_nodes.size(),vector<OwnNode_t>());
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
            extension_information << "Extending unitig from: "<<start_point._id<<":"<<start_point._val<<": "<<Common::return_unitig(sequence_map, start_point._val)<<endl;
        _extension_basic(vector_unitigs, start_point, checked, extension_information);
        if (Parameters::get().debug)
            extension_information << "End extension"<<endl;
    }
    vector<vector<OwnNode_t>> final_unitigs;
    if (!Parameters::get().greedy) {
        for (size_t i = 0; i < vector_unitigs.size(); ++i) {
            vector <OwnNode_t> tmp_unitig;
            for (auto n: vector_unitigs[i]) {
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
        if (Parameters::get().debug)
            extension_information <<"Finishing basic extension!"<< endl;
        extension_information.close();
        return final_unitigs;
    }
    if (Parameters::get().debug)
        post_process_information << "Post-process unitigs"<<endl;
    checked = vector<bool>(vector_unitigs.size(), false);
    std::function<void(size_t, size_t)> __greedy_extension =
            [this, &__greedy_extension,
             &vector_unitigs, &post_process_information,
             &final_unitigs, &checked, &sequence_map](size_t unitig_id, size_t num_call)->void
            {
                checked[unitig_id] = true;
                vector<OwnNode_t> unitig = vector_unitigs[unitig_id];
                if (unitig.size() == 0) {
                    post_process_information<<"Unitig size = 0"<<endl;
                    return;
                }
                OwnNode_t last_node = unitig[unitig.size()-1];
                if (Parameters::get().debug)
                {
                    post_process_information << "Unitig number: " << unitig_id<< " "
                                             << Common::return_unitig(sequence_map, _g_nodes[unitig_id]._val) << endl;
                    post_process_information << "ends with: " << _g_nodes[last_node]._id << ":"<< _g_nodes[last_node]._val<<" "
                                             << Common::return_unitig(sequence_map, _g_nodes[last_node]._val) << endl;
                }
                vector<OwnNode_t> neighs = getNeighbor(last_node);
                if (neighs.size() > 1)
                {
                    OwnNode_t pivot = _g_nodes[neighs[0]]._val;
                    bool same_flag_neigh = true;
                    for (auto n: neighs)
                        same_flag_neigh &= (_g_nodes[n]._val == pivot);
                    if (same_flag_neigh) {
                        for (auto n:neighs) {
                            if (Parameters::get().debug)
                                post_process_information << "continues by "<<_g_nodes[n]._id<<":"<<_g_nodes[n]._val
                                                      << Common::return_unitig(sequence_map, _g_nodes[n]._val)
                                                      << endl;
                            Op_Container::join_containers(vector_unitigs[unitig_id], vector_unitigs[n]);
                            if (num_call == 1) {
                                if (Parameters::get().debug)
                                    post_process_information << "Adding final unitig extension: "<< num_call<<" Unitig_id: "<<unitig_id<<endl;
                                final_unitigs.push_back(vector_unitigs[unitig_id]);
                            }
                        }
                    }
                }
            };
    /*
     * Look for predecessors - check!
     */
    map<OwnNode_t, vector<OwnNode_t>> predecessors;
    for (auto s:starting_points)
    {
        if (Parameters::get().debug)
            extension_information <<"Starting point: "<<s._id <<" Unitig length: "<<vector_unitigs[s._id].size()<<endl;
        OwnNode_t last_node_id = vector_unitigs[s._id][vector_unitigs[s._id].size() - 1];
        for (auto neigh: getNeighbor(last_node_id))
        {
            if (Parameters::get().debug)
                extension_information << " one way "<<_g_nodes[neigh]._id<<":"<<_g_nodes[neigh]._val<<" ";
            if (predecessors.find(last_node_id) == predecessors.end())
                predecessors[neigh] = vector<OwnNode_t>(s._id);
            else
                predecessors[neigh].push_back(s._id);
        }
        if (Parameters::get().debug)
            extension_information<<endl;
    }
    /*
     * Built stack for launchings
     */
    for (auto s:starting_points)
    {
        OwnNode_t node = s._id;
        stack<OwnNode_t> launch_nodes, assay_nodes;
        if (predecessors.find(node) == predecessors.end())
            launch_nodes.push(node);
        else
            assay_nodes.push(node);
        while (!assay_nodes.empty())
        {
            node = assay_nodes.top();
            for (auto n: predecessors[node]) {
                if (predecessors.find(n) == predecessors.end())
                    launch_nodes.push(n);
                else
                    assay_nodes.push(n);
            }
            assay_nodes.pop();
        }

        while (!launch_nodes.empty())
        {
            node = launch_nodes.top();
            if (!checked[node]) {
                if (Parameters::get().debug)
                    post_process_information << "Extending unitig from " << node <<" "
                        << Common::return_unitig(sequence_map,_g_nodes[node]._val)<<" Starting at: "<<Common::return_unitig(sequence_map, _g_nodes[s._id]._val)<< endl;
                __greedy_extension(node, 1);
            }
            launch_nodes.pop();
        }
    }
    /*
     * Translate/Plot the final unitigs
     */
    if (Parameters::get().debug)
        extension_information<<"Final unitigs:"<<endl;
    for (size_t i = 0; i < final_unitigs.size(); ++i)
    {
        vector<OwnNode_t> tmp_unitig;
        for (auto n: final_unitigs[i])
        {
            if (Parameters::get().debug)
                extension_information << _g_nodes[n]._id << ":" << _g_nodes[n]._val << " ";
            tmp_unitig.push_back(_g_nodes[n]._val);
        }
        final_unitigs[i] = tmp_unitig;
        if (Parameters::get().debug)
            extension_information << endl;
    }
    extension_information.close();
    return final_unitigs;
}

void DBG::_extension_basic(vector<vector<OwnNode_t>> & unitigs, UG_Node node, vector<bool> & checked, std::ofstream & extension_information)
{
    /*
     * Greedy method - continue until outdegree node > 1
     */
    vector<OwnNode_t> neighs = getNeighbor(node._id);
    OwnNode_t parent = NO_NEIGH;
    unitigs[node._id].push_back(node._id);
    if (Parameters::get().debug)
        extension_information << "Unitig: "<<node._id<<":"<<node._val<<"\t";
    checked[node._id] = true;
    while (neighs.size() > 0)
    {
        if (neighs.size() != 1)
        {
            if (unitigs[node._id].size() > 1)
            {
                cout << "Val: "<<_g_nodes[parent]._val<<" Id: "<<_g_nodes[parent]._id<<" Neighs: "<<neighs.size()<<endl;
                //cin.get();
            }
            break;
        }
        parent = _g_nodes[neighs[0]]._id;
        size_t parent_val = _g_nodes[neighs[0]]._val;
        neighs = getNeighbor(parent);
        unitigs[node._id].push_back(parent);
        if (Parameters::get().debug)
            extension_information<<parent<<":"<<parent_val<<"("<<neighs.size()<<") ";
        checked[parent] = true;
    }
    if (Parameters::get().debug)
        extension_information << endl;
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
        if (in_degree(i._id) == 0)
            suspicious_nodes.push_back(i);
        else if (out_degree(i._id) > 1)
            for (auto n:getNeighbor(i._id))
                suspicious_nodes.push_back(_g_nodes[n]);
        if (Parameters::get().debug)
        {
            if (out_degree(i._id) == 0)
                zero_degree << "Node: "<<_g_nodes[i._id]._val<<endl;
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
        if (steps > MAX_PATH)
            return;
        if (branches == MAX_BRANCHES)
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
     * Set frequencies for nodes based on frequencies of the edges
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
    /*for (size_t i = 0; i < _g_nodes_frequency.size()/2; ++i)
    {
        _g_nodes[i]._abundance /= (sequence_map[i].size() - Parameters::get().kmerSize +1);
        _g_nodes[i+sequence_map.size()]._abundance /= (sequence_map[i].size() - Parameters::get().kmerSize +1);
    }*/
    /*
    * Getting full reachability matrix - after polishing the matrix
    */
    _complete_reach_matrix();
}

void DBG::export_to_gfa(const vector<string> & sequence_map, string file_name)
{
    std::ofstream outfile(file_name, std::ofstream::binary);
    outfile << "H\tVN:Z:1.0"<<endl;
    for (auto v:_g_nodes)
        outfile << "S\t"<<v._id<<"\t"<<Common::return_unitig(sequence_map, v._val)<<endl;
    size_t id = 0;
    for (auto n: _g_edges)
    {
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
            outfile << "Node: " << v._id << "/" << v._val;
            outfile << " Frequency: " << _g_nodes_frequency[v._id] << " - "<<v._abundance<<" - ";
            outfile << " Number of neighbors: " << _g_edges[v._id].size() << " - ";
            for (size_t i = 0; i < _g_edges[v._id].size();++i)
            {
                OwnNode_t n = _g_edges[v._id][i];
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
            for (auto n:_g_edges[v._id])
                outfile << " " << n;
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

void UG::readjust_flow(float min_flow, float max_flow)
{
    for (size_t i = 0; i < _g_freqs.size(); ++i)
    {
        for (size_t j = 0; j < _g_freqs[i].size(); ++j) {
            //_g_freqs[i][j] = (max_flow >= _g_freqs[i][j])?max_flow - _g_freqs[i][j] + min_flow:min_flow;
            _g_freqs[i][j] = ceil((MAX_GRANULARITY_ALLOWED*100)*abs((float)_g_freqs[i][j] - max_flow)/max_flow + MAX_GRANULARITY_ALLOWED);
        }
        _g_nodes[i].readjustAbundance(min_flow, max_flow);
    }
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

void UG::addEdge(OwnNode_t i, OwnNode_t j, size_t cost)
{
    _num_edges++;
    vector<OwnNode_t>::iterator pos_j_in_i = find(_g_edges[i].begin(), _g_edges[i].end(), j);
    if ( pos_j_in_i == _g_edges[i].end())
    {
        if (!_directed) {
            _g_edges[i].push_back(j);
            _g_edges[j].push_back(i);
            _g_freqs[i].push_back(cost);
            _g_freqs[j].push_back(cost);
        } else {
            /*
             * Check frequencies - keep the highest one
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
        size_t index = std::distance(_g_edges[i].begin(), pos_j_in_i);
        if (_g_freqs[i][index] > cost)
            _g_freqs[i][index] = cost;
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
        size_t max_flow,
        bool show, unordered_map<size_t, bool> true_nodes, string file_name = "")
{
    /*
     * CORREGIR EL FLUJO CABEZON
     */
    std::ofstream nw_file;
    if (show)
        nw_file.open (file_name,std::ofstream::out);
    /*
     * Define maximal allowed granularity
     */
    float max_granularity = std::max(std::min((float)MAX_GRANULARITY_ALLOWED, (float)max_flow),(float)1.0);
    /*
     * Result paths
     */
    priority_queue<pair<size_t,vector<OwnNode_t>>> result_paths;
    unordered_map<OwnNode_t, vector<float>> pre_flow_map;
    /*
     * Truncate demanded flow
     */
    max_flow = Maths::my_round(max_flow, true);
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
    float max_capacity = 0, flow_from_target = 0;
    for (auto n:_g_nodes)
    {
        if (_g_in_edges[n._id].size() == 0) {
            max_capacity += _g_nodes[n._id]._abundance;
            sources_set.emplace(n._id);
        }
        if (_g_edges[n._id].size() == 0) {
            flow_from_target += _g_nodes[n._id]._abundance;
            targets_set.emplace(n._id);
        }
    }
    /*
     * Decide if extra_sources/extra_targets and needed
     */
    bool extra_source = (max_flow > max_capacity), extra_target = (max_flow > flow_from_target);
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
        float arc_increment = ceil((float) n._abundance / (float) max_granularity);
        if (show)
            nw_file << "Arc_increment: "<<arc_increment<<" Abundance: "<<n._abundance<<endl;
        float freq = arc_increment * max_granularity;
        for (float i = arc_increment; i <= freq; i += arc_increment)
        {
            Graph::Arc a = l_g.addArc(in_node, out_node);
            lowerMap[a] = 0;capacities[a] = arc_increment;true_arcs[a] = false;
            weights[a] = (i == arc_increment)?i/freq:pow((i-arc_increment),2)/freq;
            if (show)
                nw_file << " Weight: "<<weights[a]<<endl;
            arc_vector.push_back(a);
        }
        //Arc backward - using split into several edges
        for (float i = arc_increment; i <= freq; i += arc_increment)
        {
            Graph::Arc a = l_g.addArc(out_node, in_node);
            lowerMap[a] = 0;capacities[a] = arc_increment;true_arcs[a] = false;
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
            nw_file << "Node: "<<i<<endl;
        pair<Graph::Node, Graph::Node> n_1 = map_in_out[i];
        for (size_t j = 0; j < _g_edges[i].size(); ++j)
        {
            if (show)
                nw_file << "   Neighs: "<<j<<endl;
            auto n = _g_edges[i][j];
            pair<Graph::Node, Graph::Node> n_2 = map_in_out[n];

            float pivot_flow_tmp = (float) max_flow / (float)_g_edges[i].size();
            float freq = (float) _g_freqs[i][j];
            freq = (freq == 0)?1:freq;
            /*
             * Solution Arc
             */
            Graph::Arc solution_arc = solution_g.addArc(offset_to_solution[i], offset_to_solution[n]);
            available_flow[solution_arc] = freq;
            float arc_increment = ceil((float) freq / (float) max_granularity);
            freq = arc_increment*max_granularity;
            // Forward flow
            for (float i = arc_increment; i <= freq; i += arc_increment)
            {
                Graph::Arc a = l_g.addArc(n_1.second, n_2.first);
                offsetArcToSolutionArc[a] = solution_arc;
                arcDirection[a] = 1;
                lowerMap[a] = 0;capacities[a] = arc_increment;true_arcs[a] = true;
                weights[a] = (i == arc_increment)?i/freq:pow((i-arc_increment),2)/freq;
                if (show)
                    nw_file << " Weight: "<<weights[a]<<endl;
                arc_vector.push_back(a);
            }
            // Backward flow
            for (float i = arc_increment; i <= freq; i += arc_increment)
            {
                Graph::Arc a = l_g.addArc(n_2.first, n_1.second);
                offsetArcToSolutionArc[a] = solution_arc;
                arcDirection[a] = -1;
                lowerMap[a] = 0;capacities[a] = arc_increment;true_arcs[a] = true;
                weights[a] = (i == arc_increment)?i/freq:pow((i-arc_increment),2)/freq;
                if (show)
                    nw_file << " Weight: "<<weights[a]<<endl;
                arc_vector.push_back(a);
            }
            // Update in and out flow
            out_flow[n] += freq;
            in_flow[i] += freq;
        }
        float in_distribution = (float) Maths::my_round(in_flow[i], true) / (float) _g_edges[i].size();
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
            supply_source_s += (-exogenous_flow_in);
        } else if (exogenous_flow_in > 0)
        {
            if (show)
                nw_file << "    From OutNode to t* "<<exogenous_flow_in<<endl;
            Graph::Arc a = l_g.addArc(map_in_out[k_v.first].second, s_target);true_arcs[a] = false;
            lowerMap[a] = exogenous_flow_in;capacities[a] = exogenous_flow_in;weights[a] = 0;
            supply_sink_s += exogenous_flow_in;
        }
        if (exogenous_flow_out < 0)
        {
            if (show)
                nw_file << "    From InNode to t* "<<(-exogenous_flow_out)<<endl;
            Graph::Arc a = l_g.addArc(map_in_out[k_v.first].first, s_target);true_arcs[a] = false;
            lowerMap[a] = -exogenous_flow_out;capacities[a] = -exogenous_flow_out;weights[a] = 0;
            supply_sink_s += (-exogenous_flow_out);
        } else if (exogenous_flow_out > 0)
        {
            if (show)
                nw_file << "    From s* to InNode "<<exogenous_flow_out<<endl;
            Graph::Arc a = l_g.addArc(s_source, map_in_out[k_v.first].first);true_arcs[a] = false;
            lowerMap[a] = exogenous_flow_out;capacities[a] = exogenous_flow_out;weights[a] = 0;
            supply_source_s += exogenous_flow_out;
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
        for (Graph::OutArcIt a(l_g, s_source); a != INVALID; ++a)
        {
            if (show)
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
    max_capacity = Maths::my_round(max_capacity, true);
    max_flow = (max_flow > max_capacity)? max_flow:max_capacity;
    if (show)
        nw_file << "Demanded flow: "<<max_flow<<endl;
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
            for (Graph::ArcIt a(l_g); a != INVALID; ++a)
            {
                if (show)
                    nw_file << "Source: "<<translation_map[l_g.source(a)]<<" to: "
                        <<translation_map[l_g.target(a)]<<" Flow: "<<ns.flow(a)<<" Solution edge: "<<true_arcs[a]<<endl;
                if (true_arcs[a])
                    available_flow[offsetArcToSolutionArc[a]] += (ns.flow(a) * arcDirection[a]);
            }
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
            std::function<void(Graph::Node, vector<Graph::Arc>&, float&)> __flow_to_path =
                    [this,& __flow_to_path,& solution_g, &available_flow, &targets_set, &nw_file,&solution_to_ids](Graph::Node src,
                            vector<Graph::Arc> & cur_path, float & it_max_flow)->void {
                float max_flow = -INF;
                Graph::Arc selected;
                /*
                 * Path from the maximal available flow
                 */
                for (Graph::OutArcIt a(solution_g, src); a != INVALID; ++a)
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
                if (targets_set.find(solution_to_ids[solution_g.target(selected)]) != targets_set.end())
                {
                    cur_path.push_back(selected);
                    return;
                }
                cur_path.push_back(selected);
                __flow_to_path(solution_g.target(selected),cur_path, it_max_flow);
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
                Graph::Arc selected;
                float flow = -INF;
                nw_file << "Getting source with highest flow (swag :P): "<<endl;
                for (auto s: sources_set)
                {
                    Graph::Node source_node = offset_to_solution[s];
                    Graph::OutArcIt a(solution_g, source_node);
                    if (a == INVALID && sources_added.find(s) == sources_added.end())
                    {
                        vector <OwnNode_t> real_path(1,s);
                        result_paths.push({0,real_path});
                        sources_added.emplace(s);
                        if (show)
                            nw_file << "Source added: "<<s<<endl;
                    }
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
                __flow_to_path(solution_g.target(selected), cur_path, flow);
                if (show)
                    nw_file << "Max flow: " << flow << endl;
                vector <OwnNode_t> real_path;
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
                             << solution_node_id_target << endl;
                    }
                    if (store)
                    {
                        real_path.push_back(solution_node_id_source);
                        if (i_p == cur_path.size() - 1)
                            real_path.push_back(solution_node_id_target);
                    }
                }
                if (flow >= MIN_FLOW_PATH)
                {
                    if (show)
                        nw_file << "New path with "<<flow<<" flow."<<endl;
                    result_paths.push(pair<size_t, vector < OwnNode_t >>(flow, real_path));
                }
                if (total == 0)
                {
                    if (show)
                        nw_file << "Set of paths is already a path cover - Finishing path traversal"<<endl;
                    break;
                }
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