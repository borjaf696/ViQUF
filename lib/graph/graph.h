//
// Created by borja on 8/11/19.
//

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

#define NO_NEIGH 9999999
#define MAX_PATH 15
#define MAX_BRANCHES 10

#define CLIQUE_LIMIT 2

using namespace std;
typedef size_t OwnNode_t;
typedef unordered_set<size_t> Pairedendinformation_t;
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
     * Basic methods
     */
    void addVertex(OwnNode_t);
    void addVertexWithAbundances(OwnNode_t, size_t);
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
    vector<OwnNode_t > getNeighbors(OwnNode_t i);
    vector<Node> getNodeNeighbors(OwnNode_t);
    /*
     * Others
     */
    priority_queue<pair<size_t,vector<OwnNode_t>>> report_maximal_cliques();

    size_t degree(OwnNode_t i)
    {
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
    vector<vector<OwnNode_t >> _g_edges;
    size_t _num_edges, _num_vertex;
};
class DBG
{
public:
    /*
     * Constructors
     */
    DBG(){}
    DBG(char *);
    DBG(const DBG&);

    /*
     * Operators
     */

    DBG& operator=(const DBG && dbg)
    {
        _g = dbg._g;
        _pairedEndInformation = dbg._pairedEndInformation;
        _numedges = dbg._numedges;
        return *this;
    }

    DBG& operator=(const DBG & dbg)
    {
        _g = dbg._g;
        _pairedEndInformation = dbg._pairedEndInformation;
        _numedges = dbg._numedges;
        return *this;
    }

    size_t operator[](size_t index)
    {
        if (index > _g.size())
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
    void addPair(OwnNode_t,OwnNode_t);
    size_t addPair(OwnNode_t, Pairedendinformation_t);
    void build_process_cliques();
    /*
     * Basics
     */
    size_t vertices();
    size_t edges();

    /*
     * Show information
     */
    string toString();
    void print();
private:
    void _complete_reach_matrix();
    vector<vector<OwnNode_t>> _g;
    vector<bool> _reachibility;
    vector<unordered_set<OwnNode_t>> _pairedEndInformation;
    vector<unordered_map<OwnNode_t,size_t>> _abundance_map_node;
    size_t _numedges, _num_nodes;
};

#endif //GATB_TRIAL_GRAPH_H

