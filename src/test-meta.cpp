//
// Created by borja on 29/10/19.
//
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/bit_vectors.hpp>
#include <fstream>
#include <vector>
#include <chrono>
#include <map>
#include <string>
#include <queue>
#include <stack>
#include <set>
#ifndef GATB_TRIAL_GRAPH_H
    #include "../lib/graph/graph.h"
#endif
using namespace std;
using namespace sdsl;

typedef csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> FMIndex;
typedef rank_support_v<1> RankOnes;
typedef select_support_mcl<1> SelectOnes;

void _buildFmIndex(FMIndex & fm, char * unitigs)
{
    string fileOut = "tmp/index";
    /*
     * Construct the fmindex + storing into a file
     */
    cout << "Building FMIndex... This may take a while..."<<endl;
    auto start = chrono::steady_clock::now();
    construct(fm, unitigs, 1); // generate index
    store_to_file(fm, fileOut);

    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
    cout << "Index construction complete, index requires " << size_in_mega_bytes(fm) << " MiB." << endl;
}

void _buildBitMap(bit_vector & bDolars, RankOnes & bDolarRank, SelectOnes  & bDolarSelect, char * dolars)
{
    ifstream dolarsFile(dolars);
    string line;
    size_t count = 0;
    while (getline(dolarsFile, line)){
        std::istringstream iss(line);
        if (!(iss >> count)) { break; }
    }
    dolarsFile.close();
    bDolars = bit_vector(count+1);
    dolarsFile.open(dolars);
    while (getline(dolarsFile, line))
    {
        std::istringstream iss(line);
        int a;
        if (!(iss >> a)) { break; }
        bDolars[a] = 1;
    }
    dolarsFile.close();
    bDolarSelect = select_support_mcl<1>(&bDolars);
    bDolarRank = rank_support_v<1>(&bDolars);
    cout << "Rank: "<<bDolarRank(count)<<endl;
    //cout << "Select: "<<bDolarSelect(1)<<" Rank: "<<bDolarRank(bDolarSelect(1))<<endl;
}

void _traverseReads(char * file,const FMIndex & fm, const RankOnes & rank, const SelectOnes & select, size_t kmerSize)
{
    auto start = chrono::steady_clock::now();
    IBank* inputBank = Bank::open (file);
    Iterator<Sequence>* itSeq = inputBank->iterator();

    for (itSeq->first(); !itSeq->isDone(); itSeq->next())
    {
        /*
         * KmerIteration
         */
        Data data ((char*)itSeq->item().toString().c_str());
        Kmer<256>::ModelDirect model (kmerSize);
        Kmer<256>::ModelDirect::Iterator it (model);
        it.setData (data);
        for (it.first(); !it.isDone(); it.next())
        {
            //cout << "kmer " << model.toString(it->value()) << ",  value " << it->value() << endl;
            string query = model.toString(it->value());
            size_t occs = sdsl::count(fm, query.begin(), query.end());
            if (occs != 0)
            {
                auto locations = locate(fm, query.begin(), query.begin()+query.size());
                sort(locations.begin(), locations.end());
                for (auto l:locations)
                {
                    size_t remaining = (select(rank(l)+1)-l);
                    size_t count = 0;
                    /*
                     * What ends first? Query or Read
                     */
                    while (count++ < remaining && !it.isDone())
                        it.next();
                }
            }
        }
    }
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
}

void _buildGraph(char * file)
{
    cout << "Build DBG graph from file"<<endl;
}

int main (int argc, char* argv[])
{
    char * unitigs = argv[1], * dolars = argv[2], * file = argv[3], * graphFile = argv[5];
    size_t kmerSize = atoi(argv[4]);
    FMIndex fmIndex;
    bit_vector bDolars;
    RankOnes bDolarRank;
    SelectOnes bDolarSelect;

    _buildGraph(graphFile);
    _buildBitMap(bDolars, bDolarRank, bDolarSelect, dolars);
    _buildFmIndex(fmIndex, unitigs);
    _traverseReads(file, fmIndex,bDolarRank, bDolarSelect, kmerSize);
}
