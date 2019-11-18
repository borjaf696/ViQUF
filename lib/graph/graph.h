//
// Created by borja on 8/11/19.
//

#ifndef GATB_TRIAL_GRAPH_H
#define GATB_TRIAL_GRAPH_H
#include <unordered_set>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

#define NO_NEIGH 9999999

using namespace std;

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
    /*
     * Methods graph
     */
    vector<size_t> getNeighbor(size_t, bool = true);

    /*
     * PairedEnd
     */
    unordered_set<size_t> getPairedEndInformation(size_t);
    void addPair(size_t,size_t);
    size_t addPair(size_t, unordered_set<size_t>);

    /*
     * Basics
     */
    size_t vertices();
    size_t edges();
private:
    vector<vector<size_t>> _g;
    vector<unordered_set<size_t>> _pairedEndInformation;
    size_t _numedges;
};

#endif //GATB_TRIAL_GRAPH_H

