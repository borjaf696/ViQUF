//
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
    strs.push_back( txt.substr( initialPos, std::min( pos, txt.size() ) - initialPos + 1 ) );

    return strs.size();
}

DBG::DBG(const DBG & dbg):_g(dbg._g),_numedges(dbg._numedges),_pairedEndInformation(dbg._pairedEndInformation)
{}

DBG::DBG(char * file):_numedges(0)
{
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
            size_t numNodes;
            iss >> numNodes;
            cout << "Number of nodes: "<<numNodes<<endl;
            _g = vector<vector<size_t>>(numNodes, vector<size_t>());
            _pairedEndInformation = vector<unordered_set<size_t>>(numNodes, unordered_set<size_t>());
        }
        vector<string> results;
        split(line,results,' ');
        istringstream iss(results[0]);
        size_t node;
        iss >> node;
        if (results.size() > 1)
        {
            line = results[1];
            results.clear();
            split(line, results,',');
            for (auto s: results)
            {
                size_t neighbor;
                iss = istringstream(s);
                iss >> neighbor;
                _g[node].push_back(neighbor);
            }
        }
    }

    inFile.close();
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