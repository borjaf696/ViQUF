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

#define SPAN 128
#define ILLUMINA true
#define VERSION 2

using namespace std;
using namespace sdsl;

typedef csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> FMIndex;
typedef rank_support_v<1> RankOnes;
typedef select_support_mcl<1> SelectOnes;
struct stored_info
{
    stored_info(){}
    stored_info(size_t left,size_t right, size_t unitig):_pos(left), _remaining(right), _unitig(unitig){}
    stored_info& operator=(const stored_info & info)
    {
        _pos = info._pos;
        _remaining = info._remaining;
        _unitig = info._unitig;
        return *this;
    }
    size_t _pos, _remaining, _unitig = 0;
};

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
    while (getline(dolarsFile, line))
    {
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
    cout << "Select: "<<bDolarSelect(1)<<endl;
    cout << "Rank: "<<bDolarRank(bDolarSelect(1)+1)<<endl;
    //cout << "Select: "<<bDolarSelect(1)<<" Rank: "<<bDolarRank(bDolarSelect(1))<<endl;
}

/*
 * Hash version
 */

void _buildHash(unordered_map<Kmer<SPAN>::Type,stored_info> & kmer_map,
        char * unitigs, size_t total_unitigs, size_t kmerSize
        ,unordered_map<size_t, Sequence> & sequence_map)
{
    cout << "Building hash..."<<endl;
    auto start = chrono::steady_clock::now();
    size_t num_unitig = 0;
    /*
     * Sequence Iterator
     */
    IBank* inputBank = Bank::open (unitigs);
    ProgressIterator<Sequence> it (*inputBank, "Iterating sequences");
    /*
     * Kmer models
     */
    Kmer<SPAN>::ModelCanonical kmerModel (kmerSize);
    Kmer<SPAN>::ModelCanonical::Iterator kmerIt (kmerModel);
    for (it.first(); !it.isDone(); it.next())
    {
        Sequence& seq = it.item();
        size_t l = seq.getDataSize(), pos = 0;
        kmerIt.setData(seq.getData());
        for (kmerIt.first(); !kmerIt.isDone();kmerIt.next())
        {
            kmer_map[kmerIt->forward()] = stored_info(pos, l-pos-1, num_unitig);
            kmer_map[kmerIt->revcomp()] = stored_info(l-pos-1, pos, total_unitigs/2 + num_unitig);
            pos++;
        }
        sequence_map[num_unitig] = seq;
        num_unitig++;
    }
    cout <<"Number of unitigs: "<<num_unitig<<endl;
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds (FULL) : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
}
/*
 * Traversion alternative
 */
void _traverseReadsHash(char * file_left, char * file_right
        ,const unordered_map<Kmer<SPAN>::Type,stored_info> kmer_map
        ,size_t kmerSize, DBG & g)
{
    cout << "Traversing reads!"<< endl;
    auto start = chrono::steady_clock::now();
    /*
     * Paired_end input banks
     */
    IBank* inputBank_left = Bank::open (file_left);
    IBank* inputBank_right = Bank::open(file_right);

    PairedIterator <Sequence> *itPair = new PairedIterator<Sequence>(inputBank_left->iterator(), inputBank_right->iterator());
    ProgressIterator <std::pair<Sequence, Sequence>> progress_iter(itPair, "paired-end", inputBank_left->estimateNbItems());
    /*
     * Kmer Models
     */
    Kmer<SPAN>::ModelCanonical kmerModel (kmerSize);
    Kmer<SPAN>::ModelCanonical::Iterator kmerItLeft (kmerModel), kmerItRight (kmerModel);
    cout << "Starting reads traversion"<<endl;
    for (progress_iter.first(); !progress_iter.isDone(); progress_iter.next()) {
        Sequence &s1 = itPair->item().first;
        Sequence &s2 = itPair->item().second;
        size_t l1 = s1.getDataSize(), l2 = s2.getDataSize();
        kmerItLeft.setData(s1.getData());
        kmerItRight.setData(s2.getData());
        /*
         * Traverse kmers
         */
        kmerItRight.first();
        for (kmerItLeft.first(); !kmerItLeft.isDone(); kmerItLeft.next())
        {
            l1--;l2--;
            bool forward = true;
            unordered_map<Kmer<SPAN>::Type, stored_info>::const_iterator map_iterator_left = kmer_map.find(kmerItLeft->forward())
                    , map_iterator_right;
            if (map_iterator_left == kmer_map.end())
            {
                forward = false;
                map_iterator_left = kmer_map.find(kmerItLeft->revcomp());
            }
            if (map_iterator_left != kmer_map.end())
            {
                map_iterator_right = (forward)?kmer_map.find(kmerItRight->revcomp()):kmer_map.find(kmerItRight->forward());
                if (map_iterator_right != kmer_map.end())
                {
                    size_t unitig_left = (*map_iterator_left).second._unitig
                            ,unitig_right = (*map_iterator_right).second._unitig;
                    size_t remaining_unitig_left = (*map_iterator_left).second._remaining
                            ,remaining_unitig_right = (*map_iterator_right).second._remaining;
                    g.addPair(unitig_left, unitig_right);
                    if (max(l1,l2) < min(remaining_unitig_left, remaining_unitig_right))
                        break;
                    for (size_t i = 0; i < min(remaining_unitig_left, remaining_unitig_right)-1; ++i)
                    {
                        kmerItLeft.next();
                        //cout << kmerModel.toString(kmerItLeft->forward())<<endl;
                        kmerItRight.next();
                        l1--;l2--;
                        if (kmerItLeft.isDone() || kmerItRight.isDone())
                            break;
                    }
                }
            }
            kmerItRight.next();
            if (kmerItRight.isDone())
                break;
        }
    }
}
/*
 * Traversion
 */
void _traverseReadsFR(char * file_left, char * file_right ,const FMIndex & fm, const RankOnes & rank,
                    const SelectOnes & select, size_t kmerSize, DBG & g)
{
    cout << "Traversing reads!"<<endl;
    auto start = chrono::steady_clock::now();
    /*
     * Paired_end input banks
     */
    IBank* inputBank_left = Bank::open (file_left);
    IBank* inputBank_right = Bank::open(file_right);

    PairedIterator <Sequence> *itPair = new PairedIterator<Sequence>(inputBank_left->iterator(), inputBank_right->iterator());
    ProgressIterator <std::pair<Sequence, Sequence>> progress_iter(itPair, "paired-end", inputBank_left->estimateNbItems());

    /*
     * Kmer Models
     */
    Kmer<SPAN>::ModelCanonical kmerModel (kmerSize);
    Kmer<SPAN>::ModelCanonical::Iterator kmerItLeft (kmerModel), kmerItRight (kmerModel);
    cout << "Starting reads traversion"<<endl;
    for (progress_iter.first(); !progress_iter.isDone(); progress_iter.next())
    {
        Sequence &s1 = itPair->item().first;
        Sequence &s2 = itPair->item().second;
        size_t l1 = s1.getDataSize(), l2 = s2.getDataSize();
        kmerItLeft.setData(s1.getData());
        kmerItRight.setData(s2.getData());
        /*
         * Traverse kmers
         */
        kmerItRight.first();
        for (kmerItLeft.first(); !kmerItLeft.isDone();kmerItLeft.next())
        {
            l1--;l2--;
            string query = "", query2 ="";
            bool dirLeft = true, dirRight = true;
            query = kmerModel.toString(kmerItLeft->forward());
            query2 = kmerModel.toString(kmerItRight->forward());
            query = kmerModel.toString(kmerItLeft->revcomp());
            query2 = kmerModel.toString(kmerItRight->revcomp());
            kmerItRight.next();
            if (kmerItRight.isDone())
                break;
        }
    }
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds (TFR) : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
}

void _traverseReadsL(char * file_left, char * file_right ,const FMIndex & fm, const RankOnes & rank,
                    const SelectOnes & select, size_t kmerSize, DBG & g)
{
    cout << "Traversing reads!"<<endl;
    auto start = chrono::steady_clock::now();
    /*
     * Paired_end input banks
     */
    IBank* inputBank_left = Bank::open (file_left);
    IBank* inputBank_right = Bank::open(file_right);

    PairedIterator <Sequence> *itPair = new PairedIterator<Sequence>(inputBank_left->iterator(), inputBank_right->iterator());
    ProgressIterator <std::pair<Sequence, Sequence>> progress_iter(itPair, "paired-end", inputBank_left->estimateNbItems());

    /*
     * Kmer Models
     */
    Kmer<SPAN>::ModelCanonical kmerModel (kmerSize);
    Kmer<SPAN>::ModelCanonical::Iterator kmerItLeft (kmerModel), kmerItRight (kmerModel);
    cout << "Starting reads traversion"<<endl;
    for (progress_iter.first(); !progress_iter.isDone(); progress_iter.next())
    {
        Sequence &s1 = itPair->item().first;
        Sequence &s2 = itPair->item().second;
        size_t l1 = s1.getDataSize(), l2 = s2.getDataSize();
        kmerItLeft.setData(s1.getData());
        kmerItRight.setData(s2.getData());
        /*
         * Traverse kmers
         */
        kmerItRight.first();
        for (kmerItLeft.first(); !kmerItLeft.isDone();kmerItLeft.next())
        {
            l1--;l2--;
            string query = "", query2 ="";
            bool dirLeft = true, dirRight = true;
            query = kmerModel.toString(kmerItLeft->forward());
            query2 = kmerModel.toString(kmerItRight->forward());
            auto locations1 = locate(fm, query.begin(), query.begin()+query.size());
            if (locations1.size() == 0)
            {
                query = kmerModel.toString(kmerItLeft->revcomp());
                locations1 = locate(fm, query.begin(), query.begin()+query.size());
                dirLeft = false;
            }
            if (locations1.size() != 0)
            {
                auto locations2 = locate(fm, query2.begin(), query2.begin() + query2.size());
                if (locations2.size() == 0) {
                    query2 = kmerModel.toString(kmerItRight->revcomp());
                    locations2 = locate(fm, query2.begin(), query2.begin() + query2.size());
                    dirRight = false;
                }
                if (locations2.size() != 0)
                {
                    size_t hit1 = rank(locations1[0]), hit2 = rank(locations2[0]);
                    size_t remainingUnitigLeft = (dirLeft) ? select(hit1 + 1) - locations1[0] : locations1[0] -((hit1)?
                                                                                                                select(hit1):hit1)
                    , remainingUnitigRight = (dirRight) ? select(hit2 + 1) - locations2[0] : locations2[0] -((hit2)?
                                                                                                             select(hit2):hit2);
                    /*cout << "Query1: " << query << " Query2: " << query2 << " Remaining Unitig1: "
                         << remainingUnitigLeft << " Remaining Unitig2: " << remainingUnitigRight << endl;*/
                    if (max(l1,l2) < min(remainingUnitigLeft, remainingUnitigRight))
                        break;
                    for (size_t i = 0; i < min(remainingUnitigLeft, remainingUnitigRight); ++i)
                    {
                        kmerItLeft.next();
                        //cout << kmerModel.toString(kmerItLeft->forward())<<endl;
                        kmerItRight.next();
                        l1--;l2--;
                        if (kmerItLeft.isDone() || kmerItRight.isDone())
                            break;
                    }
                }
            }
            kmerItRight.next();
            if (kmerItRight.isDone())
                break;
        }
    }
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds (L) : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
}

void _traverseReads(char * file_left, char * file_right ,const FMIndex & fm, const RankOnes & rank,
        const SelectOnes & select, size_t kmerSize, DBG & g)
{
    cout << "Traversing reads!"<<endl;
    auto start = chrono::steady_clock::now();
    /*
     * Paired_end input banks
     */
    IBank* inputBank_left = Bank::open (file_left);
    IBank* inputBank_right = Bank::open(file_right);

    PairedIterator <Sequence> *itPair = new PairedIterator<Sequence>(inputBank_left->iterator(), inputBank_right->iterator());
    ProgressIterator <std::pair<Sequence, Sequence>> progress_iter(itPair, "paired-end", inputBank_left->estimateNbItems());

    /*
     * Kmer Models
     */
    Kmer<SPAN>::ModelCanonical kmerModel (kmerSize);
    Kmer<SPAN>::ModelCanonical::Iterator kmerItLeft (kmerModel), kmerItRight (kmerModel);
    cout << "Starting reads traversion"<<endl;
    for (progress_iter.first(); !progress_iter.isDone(); progress_iter.next())
    {
        Sequence &s1 = itPair->item().first;
        Sequence &s2 = itPair->item().second;
        size_t l1 = s1.getDataSize(), l2 = s2.getDataSize();
        kmerItLeft.setData(s1.getData());
        kmerItRight.setData(s2.getData());
        /*
         * Traverse kmers
         */
        kmerItRight.first();
        for (kmerItLeft.first(); !kmerItLeft.isDone();kmerItLeft.next())
        {
            l1--;l2--;
            string query = "", query2 ="";
            bool dirLeft = true, dirRight = true;
            query = kmerModel.toString(kmerItLeft->forward());
            query2 = kmerModel.toString(kmerItRight->forward());
            auto locations1 = locate(fm, query.begin(), query.begin()+query.size());
            if (locations1.size() == 0)
            {
                query = kmerModel.toString(kmerItLeft->revcomp());
                locations1 = locate(fm, query.begin(), query.begin()+query.size());
                dirLeft = false;
            }
            if (locations1.size() != 0)
            {
                auto locations2 = locate(fm, query2.begin(), query2.begin() + query2.size());
                if (locations2.size() == 0) {
                    query2 = kmerModel.toString(kmerItRight->revcomp());
                    locations2 = locate(fm, query2.begin(), query2.begin() + query2.size());
                    dirRight = false;
                }
                if (locations2.size() != 0)
                {
                    size_t hit1 = rank(locations1[0]), hit2 = rank(locations2[0]);
                    size_t remainingUnitigLeft = (dirLeft) ? select(hit1 + 1) - locations1[0] - (kmerSize-1)
                            : locations1[0] -((hit1)? select(hit1):hit1)
                    , remainingUnitigRight = (dirRight) ? select(hit2 + 1) - locations2[0] - (kmerSize-1)
                            : locations2[0] -((hit2)? select(hit2):hit2);
                    /*cout << "Query1: " << query << " Query2: " << query2 << " Remaining Unitig1: "
                         << remainingUnitigLeft << " Remaining Unitig2: " << remainingUnitigRight << " "<<rank(locations1[0])<<" "<<rank(locations2[0])<<endl;*/
                    if (locations2.size() == 1) {
                        size_t node1 = rank(locations1[0]), node2 = rank(locations2[0]);
                        g.addPair(node1, node2);
                    }
                    //Adjust this part
                    if (max(l1,l2) < min(remainingUnitigLeft, remainingUnitigRight))
                        break;
                    for (size_t i = 0; i < min(remainingUnitigLeft, remainingUnitigRight)-1; ++i)
                    {
                        kmerItLeft.next();
                        //cout << kmerModel.toString(kmerItLeft->forward())<<endl;
                        kmerItRight.next();
                        l1--;l2--;
                        if (kmerItLeft.isDone() || kmerItRight.isDone())
                            break;
                    }
                }
            }
            kmerItRight.next();
            if (kmerItRight.isDone())
                break;
        }
    }
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds (FULL) : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
}

/*
 * Reporting
 */
void _write_unitigs(vector<vector<size_t>> unitigs
        , char * write_path, const unordered_map<size_t, Sequence> & sequence_map
        , size_t num_unitigs_fw)
{
    cout << "Writing unitigs"<<endl;
    size_t num_unitigs = 0;
    ofstream outputFile(write_path);
    for (auto unitig:unitigs)
    {
        outputFile << ">"<<num_unitigs++<<endl;
        for (auto u:unitig)
        {
            bool reverse = false;
            if (u >= (num_unitigs_fw/2))
            {
                u -= (num_unitigs_fw/2);
                reverse = true;
            }
            Sequence seq = sequence_map.at(u);
            cout << "Reverse:"<<endl;
            cout << seq.toString()<<endl;
            //outputFile << ((reverse)?seq.getRevcomp():seq.toString())<<endl;
            outputFile << seq.toString()<<endl;
        }
    }
    cout << "End!"<< endl;
    outputFile.close();
}


void _build_process_cliques(DBG & g, const unordered_map<size_t, Sequence> & sequence_map, char * write_path)
{
    cout << "Processing cliques from DBG... This will take a while!" <<endl;
    DBG apdbg;
    g.build_process_cliques(apdbg);
    vector<vector<size_t>> unitigs = apdbg.export_unitigs();
    _write_unitigs(unitigs, write_path, sequence_map, g.vertices());
}

DBG _buildGraph(char * file)
{
    cout << "Build DBG graph from file"<<endl;
    DBG graph(file);
    return graph;
}

int main (int argc, char* argv[])
{
    cout << "Params: unitigFile, dolarsFile (placement), leftRead, rightRead, kmerSize, graphFile, unitigsFasta, unitigsfile"<<endl;
    char * unitigs = argv[1], * dolars = argv[2], * file1 = argv[3]
            , *file2 = argv[4] , * graphFile = argv[6], * unitigsFa = argv[7], * unitigs_file = argv[8];
    size_t kmerSize = atoi(argv[5]);
    FMIndex fmIndex;
    unordered_map<Kmer<SPAN>::Type,stored_info> kmer_map;
    bit_vector bDolars;
    RankOnes bDolarRank;
    SelectOnes bDolarSelect;
    /*
     * Trial
     */
    unordered_map<size_t, Sequence> sequence_map;

    DBG g = _buildGraph(graphFile);
    if (VERSION == 1)
    {
        /*
         * FM_Index version - V.0.0.1
         */
        _buildBitMap(bDolars, bDolarRank, bDolarSelect, dolars);
        _buildFmIndex(fmIndex, unitigs);
        //_traverseReadsFR(file1, file2, fmIndex,bDolarRank, bDolarSelect, kmerSize, g);
        //_traverseReadsL(file1, file2, fmIndex,bDolarRank, bDolarSelect, kmerSize, g);
        _traverseReads(file1, file2, fmIndex,bDolarRank, bDolarSelect, kmerSize, g);
    } else {
        /*
         * Hash version - V.0.0.2
         */
        _buildHash(kmer_map, unitigsFa, g.vertices(), kmerSize, sequence_map);
        _traverseReadsHash(file1, file2, kmer_map, kmerSize, g);
        _build_process_cliques(g, sequence_map, unitigs_file);
    }
    //g.print();
}
