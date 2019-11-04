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

using namespace std;
using namespace sdsl;

typedef csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> FMIndex;
typedef rank_support_v<1> RankOnes;
void _buildFmIndex(FMIndex & fm, char * unitigs)
{
    string fileOut = "tmp/index";
    /*
     * Construct the fmindex + storing into a file
     */
    auto start = chrono::steady_clock::now();
    construct(fm, unitigs, 1); // generate index
    store_to_file(fm, fileOut);

    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
    cout << "Index construction complete, index requires " << size_in_mega_bytes(fm) << " MiB." << endl;
}

void _buildBitMap(bit_vector & bDolars, RankOnes & bDolarRank, char * dolars)
{
    bit_vector();
    ifstream dolarsFile(dolars);
    string line;
    while (getline(dolarsFile, line))
    {
        std::istringstream iss(line);
        int a;
        if (!(iss >> a)) { break; }
    }
}

void _traverseReads(char * file,const FMIndex & fm, size_t kmerSize)
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
        Kmer<128>::ModelCanonical model (kmerSize);
        Kmer<128>::ModelCanonical::Iterator it (model);
        it.setData (data);
        for (it.first(); !it.isDone(); it.next())
        {
            //cout << "kmer " << model.toString(it->value()) << ",  value " << it->value() << endl;
            string query = model.toString(it->value());
            size_t occs = sdsl::count(fm, query.begin(), query.end());
        }
    }
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
}

int main (int argc, char* argv[])
{
    char * unitigs = argv[1], * dolars = argv[2], * file = argv[3];
    size_t kmerSize = atoi(argv[4]);
    FMIndex fmIndex;
    bit_vector bDolars;
    RankOnes bDolarRank;

    _buildBitMap(bDolars, bDolarRank, dolars);
    _buildFmIndex(fmIndex, unitigs);
    _traverseReads(file, fmIndex, kmerSize);
}
