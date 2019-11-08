//
// Created by borja on 8/11/19.
//
#ifndef GATB_TRIAL_GRAPH_H
#include "graph.h"
#endif

DBG::DBG(char * file):_numedges(0)
{

}

size_t DBG::edges()
{
    return _numedges;
}


size_t DBG::vertices()
{
    return _g.size();
}

size_t DBG::getNeighbor(size_t v , bool dir)
{
    if (dir)
        return _g[v];
}

size_t DBG::getPairedEndInformation(size_t v)
{
    return _pairedEndInformation[v];
}

void DBG::addPair(size_t v, size_t p)
{
    _pairedEndInformation[v].emplace(p);
}

void DBG::addPair(size_t v, unordered_set <size_t> pairInformation)
{
    _pairedEndInformation[v] = pairInformation;
    return pairInformation.size();
}