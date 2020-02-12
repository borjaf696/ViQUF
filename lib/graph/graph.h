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
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>

// MaxPATH = 3 - MAX_BRANCHES 6
#define NO_NEIGH 9999999
#define MAX_PATH 20
#define MAX_BRANCHES 30
#define COMPLETE false

#define CLIQUE_LIMIT 3

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
template<typename T>
void display_set(unordered_set<T> s1)
{
    for (auto s: s1)
        cout << " "<<s;
    cout <<endl;
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
        UG_Node(OwnNode_t val):_val(val)
        {
            _paired_info = Pairedendinformation_t();
            _map_abundance = unordered_map<OwnNode_t,size_t>();
        }

        UG_Node(size_t id, OwnNode_t val):_val(val), _id(id)
        {
            _paired_info = Pairedendinformation_t();
            _map_abundance = unordered_map<OwnNode_t,size_t>();
        }

        UG_Node(OwnNode_t val, Pairedendinformation_t paired_info):_val(val),_paired_info(paired_info)
        {
            _map_abundance = unordered_map<OwnNode_t,size_t>();
        }

        UG_Node& operator=(const UG_Node & node)
        {
            _val = node._val;
            _id = node._id;
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

        void add_paired_information(OwnNode_t val)
        {
            _paired_info.emplace(val);
            if (_map_abundance.find(val) != _map_abundance.end())
                _map_abundance[val]++;
            else
                _map_abundance[val] = 1;
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

        OwnNode_t _val;
        size_t _id;
        Pairedendinformation_t _paired_info;
        unordered_map<OwnNode_t,size_t> _map_abundance;
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
        }
        for (auto i: index_2)
        {
            if (_g_nodes[i] == node_son)
            {
                i_2 = _g_nodes[i]._id;
                break;
            }
        }
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
            _g_in_edges.push_back(vector<OwnNode_t>());
        }
        _g_edges[i_1].push_back(i_2);
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
        _numedges = dbg._numedges;
        _num_nodes = dbg._num_nodes;
        return *this;
    }

    DBG& operator=(const DBG & dbg)
    {
        _g_nodes = dbg._g_nodes;
        _g_edges = dbg._g_edges;
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
    vector<size_t> getNeighbor(OwnNode_t, bool = true);
    /*
     * PairedEnd
     */
    Pairedendinformation_t getPairedEndInformation(OwnNode_t);
    size_t getAbundance(OwnNode_t, OwnNode_t);
    void addPair(OwnNode_t,OwnNode_t);
    size_t addPair(OwnNode_t, Pairedendinformation_t);
    void build_process_cliques(DBG &);
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
    vector<vector<OwnNode_t>> export_unitigs();
    /*
     * Show information
     */
    string toString();
    void stats();
    void print();
private:
    vector<OwnNode_t> _get_posible_pos(OwnNode_t node_val)
    {
        if (_map_pos.find(node_val) != _map_pos.end())
            return _map_pos[node_val];
        _map_pos[node_val] = vector<OwnNode_t>();
        return vector<OwnNode_t>();
    }
    void _complete_reach_matrix();
    bool _reach(OwnNode_t, OwnNode_t);
    void _print_reachability();
    void _subsane();
    /*
     * Extension
     */
    OwnNode_t _get_possible_extension(OwnNode_t);
    vector<UG_Node> _get_starting_nodes();
    void _extension(queue<UG_Node>&, vector<bool>&, vector<OwnNode_t>&, UG_Node,
            size_t&, size_t&,size_t&, size_t&, size_t&,vector<vector<OwnNode_t>>&,vector<bool>&);
    vector<UG_Node> _g_nodes;
    vector<vector<OwnNode_t>> _g_edges, _g_in_edges;
    vector<unordered_set<OwnNode_t>> _reachability;
    size_t _numedges, _num_nodes;
    unordered_map<OwnNode_t, vector<OwnNode_t>> _map_pos;
};

class UG
{
public:
    struct UG_Node
    {
        UG_Node(size_t id, OwnNode_t val):_val(val), _id(id){}
        UG_Node(size_t id, OwnNode_t val, size_t abundance):_val(val), _id(id), _abundance(abundance){}
        UG_Node& operator=(const UG_Node & node)
        {
            _val = node._val;
            _id = node._id;
            _abundance = node._abundance;
            return *this;
        }
        OwnNode_t _val;
        size_t _id, _abundance;
    };
    typedef UG_Node Node;
    UG():_num_edges(0), _num_vertex(0){}
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
    OwnNode_t addVertexWithAbundances(OwnNode_t, size_t);
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
     * Others
     */
    priority_queue<pair<size_t,vector<OwnNode_t>>> report_maximal_cliques(unordered_map<size_t,size_t>, DBG &);

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
};
#endif //GATB_TRIAL_GRAPH_H

