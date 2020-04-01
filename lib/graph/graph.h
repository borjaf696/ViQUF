//
// Created by borja on 8/11/19.
//
#pragma once
#ifndef GATB_TRIAL_GRAPH_H
#define GATB_TRIAL_GRAPH_H
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <set>
#include <functional>
#include <algorithm>
#include <cmath>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <stdlib.h>
#ifndef VIADBG_TRIAL_EXTRA_H
    #include "../extra/extra.h"
#endif
#include <lemon/smart_graph.h>
#include <lemon/network_simplex.h>

// MaxPATH = 3 - MAX_BRANCHES 6
#define INF 9999999
#define NO_NEIGH 9999999
#define MAX_PATH 1000
#define MAX_BRANCHES 1000
//Maximum distance for directed graph
#define D_MAX_PATH 50
#define D_MAX_BRANCHES 50
#define COMPLETE 0
#define CLIQUE_LIMIT 2
#define SIZE_RELATION 0.2

#define STD_LIMIT 2
#define RATE 0.5
#define MIN_RATE 0.25

#define RATIO 10
#define RATIO_SIMILARITY 50

#define FIRSTLASTLENGTH 100

using namespace lemon;
using namespace std;

typedef size_t OwnNode_t;

typedef unordered_set<OwnNode_t> Pairedendinformation_t;
template<typename T>
bool is_subset(unordered_set<T> s1, unordered_set<T> s2)
{
    if (s1.empty() | s2.empty())
        return false;
    for (auto t: s1)
    {
        if (s2.find(t)==s2.end())
            return false;
    }
    return true;
}
template<typename T>
unordered_set<T> intersect(unordered_set<T> s1, unordered_set<T> s2)
{
    unordered_set<T> intersection;
    if (s1.empty() | s2.empty())
        return intersection;
    for (auto t:s1)
    {
        if (s2.find(t) != s2.end())
            intersection.emplace(t);
    }
    return intersection;
}
template<typename T>
unordered_set<T> union_set(unordered_set<T> s1, unordered_set<T> s2)
{
    unordered_set<T> union_s = s1;
    for (auto s: s2)
        union_s.emplace(s);
    return union_s;
}
template<typename T>
unordered_set<T> union_sets(vector<unordered_set<T>> v)
{
    unordered_set<T> union_s = v[0];
    for (auto s: v)
        for (auto e:s)
            union_s.emplace(e);
    return union_s;
}
template <typename T>
unordered_set<T> minus_set(unordered_set<T> s1, unordered_set<T> s2)
{
    unordered_set<T> minus_s = s1;
    for (auto s: s2)
        minus_s.erase(s);
    return minus_s;
}
template<typename T>
void display_set(unordered_set<T> s1)
{
    for (auto s: s1)
        cout << " "<<s;
    cout <<endl;
}
template<typename T>
bool in(unordered_set<T> s1, T el)
{
    return (s1.find(el)!=s1.end());
}


/*
 * future strategy pattern
template <typename T>
class Graph
{
public:
    Graph(){}
    virtual T getNeighbor(T, bool);
    virtual unordered_set<T> getPairedEndInformation(T);
    virtual size_t vertices();
    virtual size_t edges();
};
*/
class DBG
{
public:
    struct UG_Node
    {
        UG_Node(OwnNode_t val):_val(val),_abundance(0)
        {
            _paired_info = Pairedendinformation_t();
            _map_abundance = unordered_map<OwnNode_t,size_t>();
        }

        UG_Node(size_t id, OwnNode_t val):_val(val), _id(id),_abundance(0)
        {
            _paired_info = Pairedendinformation_t();
            _map_abundance = unordered_map<OwnNode_t,size_t>();
        }

        UG_Node(OwnNode_t val, Pairedendinformation_t paired_info):_val(val),_abundance(0),_paired_info(paired_info)
        {
            _map_abundance = unordered_map<OwnNode_t,size_t>();
        }

        UG_Node& operator=(const UG_Node & node)
        {
            _val = node._val;
            _id = node._id;
            _length = node._length;
            _abundance = node._abundance;
            _paired_info = node._paired_info;
            _map_abundance = node._map_abundance;
            return *this;
        }

        bool operator==(const UG_Node & other)
        {
            return (other._val == _val) && (other._paired_info == _paired_info);
        }

        void setId(size_t id)
        {
            _id = id;
        }

        void setLength(size_t length)
        {
            _length = length;
        }

        void add_paired_information(OwnNode_t val)
        {
            _paired_info.emplace(val);
            if (_map_abundance.find(val) != _map_abundance.end())
                _map_abundance[val]++;
            else
                _map_abundance[val] = 1;
            _mean_abundance++;
        }

        void add_paired_information(Pairedendinformation_t pairs)
        {
            _paired_info = pairs;
        }

        Pairedendinformation_t get_paired_information()
        {
            return _paired_info;
        }

        size_t get_abundance(OwnNode_t node)
        {
            if (_map_abundance.find(node) == _map_abundance.end())
                return 0;
            return _map_abundance[node];
        }

        void print_info()
        {
            cout << "Val: "<<_val<<" Id: "<<_id<<" Abundance: "<<_abundance<<endl;
            for (auto i: _paired_info)
                cout << " "<<i;
            cout << endl;
        }

        /*
         * Pendiente de revision profunda - parece que cumple por ahora
         */
        void post_process_pairs(const vector<string> & sequence_map)
        {
            // For quasis - 3
            // For meta - 2
            bool show = (_val == NO_NEIGH);
            size_t nodo = (_val >= sequence_map.size())?_val- sequence_map.size():_val;
            //cout << "Node: "<<_val<<" Seq: "<<sequence_map[nodo]<<endl;
            float pivot = 0.0, pivot_less = INF;
            /*for (auto p:_map_abundance)
            {
                if (p.second > pivot)
                    pivot = (float) p.second;
                if (p.second < pivot_less)
                    pivot_less = (float) p.second;
            }*/
            _mean_abundance = 0.0;
            for (auto p:_map_abundance) {
                size_t place = (p.first >= sequence_map.size())?p.first - sequence_map.size():p.first;
                _mean_abundance += ((float) p.second / (float) sequence_map.at(place).length());
            }
            _mean_abundance /= _map_abundance.size();
            pivot = _mean_abundance*RATE;pivot_less = _mean_abundance*MIN_RATE;
            float tmp_diff = 0.0, num_pairs = 0.0;
            _mean_abundance = 0.0;
            for (auto p:_map_abundance) {
                size_t place = (p.first >= sequence_map.size())?p.first - sequence_map.size():p.first;
                float local_score = ((float) p.second / (float) sequence_map.at(place).length());
                if (local_score >= pivot || local_score <= pivot_less)
                    continue;
                _mean_abundance += ((float) p.second / (float) sequence_map.at(place).length());
                num_pairs++;
            }
            _mean_abundance /= num_pairs;
            num_pairs = 0.0;
            for (auto p: _map_abundance) {
                size_t place = (p.first >= sequence_map.size())?p.first - sequence_map.size():p.first;
                float local_score = ((float) p.second / (float) sequence_map.at(place).length());
                if (show)
                {
                    cout <<"Nodo: "<<p.first<<" "<<p.second<<" "<<local_score<<endl;
                }
                if (local_score >= pivot || local_score <= pivot_less)
                    continue;
                tmp_diff += pow(_mean_abundance - local_score, 2);
                num_pairs++;
            }
            _std = sqrt(tmp_diff / (num_pairs-1));
            /*
             * Process _paired_info
             */
            float LCL = _mean_abundance-STD_LIMIT*_std, UCL = _mean_abundance+STD_LIMIT*_std;
            if (show)
            {
                cout << "Pivote: "<<pivot<<" Pivote less: "<<pivot_less<<" Mean: "<<_mean_abundance<<" SD: "<<_std<<endl;
                cout << "LCL: "<<LCL<<endl;
            }
            for (auto p:_map_abundance)
            {
                size_t place = (p.first >= sequence_map.size())?p.first - sequence_map.size():p.first;
                float local_score = ((float) p.second / (float) sequence_map.at(place).length());
                if (local_score < LCL)
                {
                    //cout << "Remove: "<<p.first<<" "<<p.second<<endl;
                    _paired_info.erase(p.first);
                }else
                {
                    if (show)
                        cout << "Remove: "<<p.first<<" "<<p.second<<endl;
                }
            }
        }
        void increase_abundance(size_t quantity)
        {
            _abundance+=quantity;
        }

        OwnNode_t _val;
        size_t _id, _length;
        float _abundance;
        Pairedendinformation_t _paired_info;
        unordered_map<OwnNode_t,size_t> _map_abundance;
        float _mean_abundance = 0.0, _std = 0.0;
    };
    /*
     * Constructors
     */
    DBG():_num_nodes(0),_numedges(0){}
    DBG(char *);
    DBG(const DBG&);

    /*
     * Methods modify graph
     */
    void addNode(UG_Node node, bool show = false)
    {
        vector<OwnNode_t> index_1 = _get_posible_pos(node._val);
        OwnNode_t i_1 = NO_NEIGH;
        /*
         * Check index positions
         */
        for (auto i: index_1)
        {
            if (_g_nodes[i] == node)
            {
                i_1 = i;
                break;
            }
        }
        if (i_1 == NO_NEIGH)
            index_1.clear();
        if (index_1.size() == 0)
        {
            i_1 = _num_nodes;
            _map_pos[node._val].push_back(i_1);
            node.setId(_num_nodes++);
            _g_nodes.push_back(node);
            _g_edges.push_back(vector<OwnNode_t>());
            _g_edges_reads.push_back(vector<size_t>());
            _g_in_edges.push_back(vector<OwnNode_t>());
            if (show)
                cout << node._val<<" "<<_g_nodes[_num_nodes-1]._val<<" "<<node._id<<" "<<_g_nodes[_num_nodes-1]._id<<endl;
        }
    }

    void addNode(UG_Node node_parent, UG_Node node_son)
    {
        vector<OwnNode_t> index_1 = _get_posible_pos(node_parent._val),index_2 = _get_posible_pos(node_son._val);
        OwnNode_t i_1 = NO_NEIGH, i_2 = NO_NEIGH;
        /*
         * Check index positions
         */
        for (auto i: index_1)
        {
            if (_g_nodes[i] == node_parent)
            {
                i_1 = _g_nodes[i]._id;
                break;
            }
            /*
             * Revisar - subconjuntos no siempre son iguales los cliques
             */
            /*if (!(_g_nodes[i] == node_parent) && (is_subset(_g_nodes[i]._paired_info, node_parent._paired_info)
                || is_subset(node_parent._paired_info, _g_nodes[i]._paired_info)))
            {
                Pairedendinformation_t insert =  node_parent._paired_info, study = _g_nodes[i]._paired_info;
                double ratio1 = ((insert.size()*100)/(insert.size() + study.size())),
                    ratio2 = ((study.size()*100)/(insert.size()+study.size()));
                if (ratio1 >= RATIO_SIMILARITY || ratio2 >= RATIO_SIMILARITY) {
                    _g_nodes[i]._paired_info = (insert.size() > study.size())?insert:study;
                    i_1 = _g_nodes[i]._id;
                    break;
                }
            }*/
        }
        for (auto i: index_2)
        {
            if (_g_nodes[i] == node_son)
            {
                i_2 = _g_nodes[i]._id;
                break;
            }
            /*if (!(_g_nodes[i] == node_son) && (is_subset(_g_nodes[i]._paired_info, node_son._paired_info)
                                                  || is_subset(node_son._paired_info, _g_nodes[i]._paired_info)))
            {
                Pairedendinformation_t insert =  node_son._paired_info, study = _g_nodes[i]._paired_info;
                double ratio1 = ((insert.size()*100)/(insert.size() + study.size())),
                        ratio2 = ((study.size()*100)/(insert.size()+study.size()));
                if (ratio1 >= RATIO_SIMILARITY || ratio2 >= RATIO_SIMILARITY) {
                    _g_nodes[i]._paired_info = (insert.size() > study.size())?insert:study;
                    i_2 = _g_nodes[i]._id;
                    break;
                }
            }*/
        }
        /*if (node_parent._val == 1011 || node_son._val == 1011)
        {
            cout <<" Caso 1011: "<<endl;
            cout << "Valor indice: "<<i_1<<endl;
            if (i_1 != NO_NEIGH)
            {
                cout <<"Indice encontrado: "<<i_1<<endl;
            } else {
                cout << "Indice usado: "<<_num_nodes<<endl;
            }
        }*/
        if (i_1 == NO_NEIGH)
            index_1.clear();
        if (i_2 == NO_NEIGH)
            index_2.clear();
        if (index_1.size() == 0)
        {
            i_1 = _num_nodes;
            _map_pos[node_parent._val].push_back(i_1);
            node_parent.setId(_num_nodes++);
            _g_nodes.push_back(node_parent);
            _g_edges.push_back(vector<OwnNode_t>());
            _g_edges_reads.push_back(vector<size_t>());
            _g_in_edges.push_back(vector<OwnNode_t>());
            _numedges++;
        }
        if (index_2.size() == 0)
        {
            i_2 = _num_nodes;
            _map_pos[node_son._val].push_back(i_2);
            node_son.setId(_num_nodes++);
            _g_nodes.push_back(node_son);
            _g_edges.push_back(vector<OwnNode_t>());
            _g_edges_reads.push_back(vector<size_t>());
            _g_in_edges.push_back(vector<OwnNode_t>());
        }
        _g_edges[i_1].push_back(i_2);
        _g_edges_reads[i_1].push_back(0);
        _g_in_edges[i_2].push_back(i_1);
        //cout <<" Adding edge from "<<i_1<<"/"<<_g_nodes[i_1]._val<<" to "<<i_2<<"/"<<_g_nodes[i_2]._val<<endl;
        // Lets check the information
        //print();
    }

    /*
     * Operators
     */
    DBG& operator=(const DBG && dbg)
    {
        _g_nodes = dbg._g_nodes;
        _g_edges = dbg._g_edges;
        _g_in_edges = dbg._g_in_edges;
        _g_edges_reads = dbg._g_edges_reads;
        _numedges = dbg._numedges;
        _num_nodes = dbg._num_nodes;
        return *this;
    }

    DBG& operator=(const DBG & dbg)
    {
        _g_nodes = dbg._g_nodes;
        _g_edges = dbg._g_edges;
        _g_in_edges = dbg._g_in_edges;
        _g_edges_reads = dbg._g_edges_reads;
        _numedges = dbg._numedges;
        _num_nodes = dbg._num_nodes;
        return *this;
    }

    size_t operator[](size_t index)
    {
        if (index > _g_nodes.size())
        {
            cout << "Array index out of bound, exiting";
            exit(0);
        }
        return index;
    }

    /*
     * Methods graph
     */
    vector<size_t> getNeighbor(OwnNode_t, bool = true, size_t = NO_NEIGH);
    void add_read(OwnNode_t, size_t);
    void add_read(OwnNode_t, OwnNode_t, size_t, bool);
    void set_length(OwnNode_t, size_t);
    void subsane(const vector<string>&);
    /*
     * PairedEnd
     */
    size_t getLength(OwnNode_t);
    float getFreqNode(OwnNode_t);
    float getAbundance(OwnNode_t);
    Pairedendinformation_t getPairedEndInformation(OwnNode_t);
    size_t getAbundance(OwnNode_t, OwnNode_t);
    void addPair(OwnNode_t,OwnNode_t, bool = false, bool = false);
    size_t addPair(OwnNode_t, Pairedendinformation_t);
    void build_process_cliques(DBG &,const vector<string> &, size_t);
    /*
     * Basics
     */
    size_t vertices();
    size_t edges();
    size_t out_degree(OwnNode_t);
    size_t in_degree(OwnNode_t);
    /*
     * Unitigs
     */
    vector<vector<OwnNode_t>> export_unitigs(bool = false);
    /*
     * Show information
     */
    string toString();
    void stats();
    void print();
    void post_process_pairs(const vector<string> &);
private:
    size_t _find_edge(OwnNode_t, OwnNode_t, bool = true);
    vector<OwnNode_t> _get_posible_pos(OwnNode_t node_val)
    {
        if (_map_pos.find(node_val) != _map_pos.end())
            return _map_pos[node_val];
        _map_pos[node_val] = vector<OwnNode_t>();
        return vector<OwnNode_t>();
    }
    void _complete_reach_matrix();
    bool _reach_e(OwnNode_t, OwnNode_t);
    bool _reach(OwnNode_t, OwnNode_t);
    void _print_reachability();
    /*
     * Extension
     */
    OwnNode_t _get_possible_extension(OwnNode_t);
    vector<UG_Node> _get_starting_nodes();
    void _extension(queue<UG_Node>&, vector<bool>&, vector<OwnNode_t>&, UG_Node,
            size_t&, size_t&,size_t&, size_t&, size_t&,size_t&,vector<vector<OwnNode_t>>&,
                    vector<bool>&, bool = false);
    vector<UG_Node> _g_nodes;
    vector<float>_g_nodes_frequency;
    vector<vector<OwnNode_t>> _g_edges, _g_in_edges;
    vector<vector<size_t>> _g_edges_reads;
    vector<unordered_set<OwnNode_t>> _reachability;
    size_t _numedges, _num_nodes;
    unordered_map<OwnNode_t, vector<OwnNode_t>> _map_pos;
};

class UG
{
public:
    /*
    * Build the environment to min cost flow
    */
    using Weight = int;
    using Capacity = int;

    using Graph = SmartDigraph;

    template<typename ValueType>
    using ArcMap = SmartDigraph::ArcMap<ValueType>;

    using NS = NetworkSimplex<SmartDigraph, Capacity, Weight>;

    struct UG_Node
    {
        UG_Node(size_t id, OwnNode_t val):_val(val), _id(id){}
        UG_Node(size_t id, OwnNode_t val, float abundance):_val(val), _id(id), _abundance(abundance){}
        UG_Node& operator=(const UG_Node & node)
        {
            _val = node._val;
            _id = node._id;
            _abundance = node._abundance;
            return *this;
        }
        OwnNode_t _val;
        size_t _id;
        float _abundance;
    };
    typedef UG_Node Node;
    UG():_num_edges(0), _num_vertex(0), _directed(false){}
    UG(bool directed):_num_edges(0),_num_vertex(0),_directed(directed){}
    /*
     * Operators
     */
    UG_Node operator[](size_t index)
    {
        if (index > _g_nodes.size())
        {
            cout << "Array index out of bound, exiting";
            exit(0);
        }
        return _g_nodes[index];
    }
    /*
     * Basic methods
     */
    OwnNode_t addVertex(OwnNode_t);
    OwnNode_t addVertexWithAbundances(OwnNode_t, float);
    bool setAbundance(OwnNode_t, size_t);
    void addEdge(OwnNode_t, OwnNode_t);
    size_t edges()
    {
        return _num_edges;
    }
    size_t vertices()
    {
        return _num_vertex;
    }
    /*
     * Neighbors - return directly the place you should look at
     */
    unordered_set<OwnNode_t > getNeighbors(OwnNode_t i);
    vector<Node> getNodeNeighbors(OwnNode_t);
    /*
     * Cliques + NF methods
     */
    priority_queue<pair<size_t,vector<OwnNode_t>>>
            report_maximal_cliques(unordered_map<size_t,size_t>, DBG &, unordered_map<size_t,float> &, bool );
    void report_min_cost_flow(const vector<size_t> &, bool);

    size_t degree(OwnNode_t i)
    {
        if (i >= _g_edges.size())
            throw std::invalid_argument( "Out of bounds");
        return _g_edges[i].size();
    }

    void print()
    {
        for (auto n:_g_nodes)
        {
            cout << "Id: "<<n._id<<" Val: "<<n._val<<" Abundance: "<<n._abundance<<endl<<" Neighs: ";
            for (auto neigh: _g_edges[n._id])
                cout << " "<<neigh;
            cout << endl;
        }
    }
private:
    vector<Node> _g_nodes;
    vector<unordered_set<OwnNode_t >> _g_edges;
    size_t _num_edges, _num_vertex;
    bool _directed;
};
#endif //GATB_TRIAL_GRAPH_H

