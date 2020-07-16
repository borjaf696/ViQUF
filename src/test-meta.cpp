//
// Created by borja on 29/10/19.
//
// We include what we need for the test
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
#define MIN_LENGTH 500

using namespace std;
using namespace sdsl;

typedef csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> FMIndex;
typedef rank_support_v<1> RankOnes;
typedef select_support_mcl<1> SelectOnes;
struct stored_info
{
    stored_info(){}
    stored_info(size_t left,size_t right, size_t unitig, size_t u_length)
            :_pos(left), _remaining(right), _unitig(unitig), __u_length(u_length){}
    stored_info& operator=(const stored_info & info)
    {
        _pos = info._pos;
        _remaining = info._remaining;
        _unitig = info._unitig;
        __u_length = info.__u_length;
        return *this;
    }

    void print()
    {
        cout << " "<<_unitig<<"-"<<_pos<<"-"<<_remaining<<"-"<<__u_length<<endl;
    }
    size_t _pos, _remaining, _unitig = 0, __u_length = 0;
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
        ,vector<string> & sequence_map, DBG & g)
{
    cout << "Building hash..."<<endl;
    cout << "File: "<<unitigs<<endl;
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
        g.set_length(num_unitig,l);
        g.set_length(num_unitig + total_unitigs/2, l);
        kmerIt.setData(seq.getData());
        for (kmerIt.first(); !kmerIt.isDone();kmerIt.next())
        {
            if (kmer_map.find(kmerIt->forward()) != kmer_map.end() || kmer_map.find(kmerIt->revcomp()) != kmer_map.end())
            {
                if (pos != 0 && pos != l - kmerSize)
                {
                    kmer_map[kmerIt->forward()] = stored_info(pos, l - pos - 1, num_unitig, l);
                    kmer_map[kmerIt->revcomp()] = stored_info(l - pos - 1, pos, total_unitigs / 2 + num_unitig,l);
                    pos++;
                }
            }else {
                kmer_map[kmerIt->forward()] = stored_info(pos, l - pos - 1, num_unitig,l);
                kmer_map[kmerIt->revcomp()] = stored_info(l - pos - 1, pos, total_unitigs / 2 + num_unitig,l);
                pos++;
            }
        }
        sequence_map.push_back(seq.toString());
        num_unitig++;
    }
    /*// Showing unitigs
    cout << "Unitigs: "<<endl;
    for (auto u: sequence_map)
    {
        cout << "   Number: "<<u.first<<" ";
        cout << "   Sequence: "<<u.second<<endl;
    }
    // Displaying kmers
    cout <<"Kmers: "<<endl;
    for (auto u: kmer_map)
    {
        cout <<" Kmer: "<<kmerModel.toString(u.first)<<endl;
        u.second.print();
    }*/
    cout <<"Number of unitigs: "<<num_unitig<<endl;
    cout <<"Number of k-mers: "<<kmer_map.size()<<endl;
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds (FULL) : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
}
/*
 * Add frequencies
 */
void _add_frequencies(const unordered_map<Kmer<SPAN>::Type,stored_info> & kmer_map,
                char * append_file, size_t total_unitigs, size_t kmerSize
        ,vector<string> & sequence_map, DBG & g)
{
    cout << "Traversing merged file: "<<append_file<< endl;
    auto start = chrono::steady_clock::now();
    /*
     * Sequence Iterator
     */
    IBank* inputBank = Bank::open (append_file);
    ProgressIterator<Sequence> it (*inputBank, "Iterating sequences");
    /*
     * Kmer models
     */
    Kmer<SPAN>::ModelCanonical kmerModel (kmerSize);
    Kmer<SPAN>::ModelCanonical::Iterator kmerIt (kmerModel);
    for (it.first(); !it.isDone(); it.next())
    {
        /*
         * Edges frequency:
         *  * f_s_l -> forward strain left? (cur_unitig_left)
         *  * f_s_r -> forward strain right? (cur_unitig_right)
         */
        size_t cur_unitig_left = NO_NEIGH, quantity = 1;
        bool f_s_l = true;
        Sequence& seq = it.item();
        size_t l = seq.getDataSize(), pos = 0;
        //cout <<"Sequence: "<<seq.toString()<<endl;
        kmerIt.setData(seq.getData());
        for (kmerIt.first(); !kmerIt.isDone();kmerIt.next()) {
            bool forward = true;
            unordered_map<Kmer<SPAN>::Type, stored_info>::const_iterator map_iterator_left = kmer_map.find(kmerIt->forward());
            if (map_iterator_left == kmer_map.end())
            {
                forward = false;
                map_iterator_left = kmer_map.find(kmerIt->revcomp());
            }
            if (map_iterator_left != kmer_map.end()) {
                size_t unitig_left = (*map_iterator_left).second._unitig;
                //cout <<"Find: "<<kmerModel.toString(kmerIt->forward())<<" "<<unitig_left<<endl;
                size_t rev_comp_unitig_left = (unitig_left >= total_unitigs/2)?unitig_left-total_unitigs/2:unitig_left+total_unitigs/2;
                /*if (cur_unitig_left == 924 || cur_unitig_left == (924 - total_unitigs/2))
                    cout << cur_unitig_left<<" "<<unitig_left<<endl;*/
                if (cur_unitig_left != NO_NEIGH) {
                    /*cout << "Kmer Nuevo: "<<kmerModel.toString(kmerIt->forward())<<" "<<unitig_left<<endl;
                    cout << "Current: "<<cur_unitig_left<<" Unitig: "<<unitig_left<<endl;*/
                    if (unitig_left != (cur_unitig_left + total_unitigs / 2) && unitig_left != cur_unitig_left) {
                        size_t cur_rev_comp = (cur_unitig_left >= total_unitigs / 2) ?
                                              cur_unitig_left - total_unitigs / 2 : cur_unitig_left + total_unitigs / 2;
                        if (f_s_l) {
                            if (forward) {
                                g.add_read(cur_unitig_left, unitig_left, quantity, true);
                                g.add_read(rev_comp_unitig_left, cur_rev_comp, quantity, false);
                            } else {
                                g.add_read(cur_unitig_left, rev_comp_unitig_left, quantity, true);
                                g.add_read(unitig_left, cur_rev_comp, quantity, false);
                            }
                        } else {
                            if (forward) {
                                g.add_read(cur_rev_comp, unitig_left, quantity, true);
                                g.add_read(rev_comp_unitig_left, cur_unitig_left, quantity, false);
                            } else {
                                g.add_read(cur_rev_comp, rev_comp_unitig_left, quantity, true);
                                g.add_read(unitig_left, cur_unitig_left, quantity, false);
                            }
                        }
                        cur_unitig_left = unitig_left;
                        f_s_l = forward;
                        quantity = 1;
                    } else
                        quantity++;
                } else {
                    cur_unitig_left = unitig_left;
                    //cout << "Kmer Inicial: "<<kmerModel.toString(kmerIt->forward())<<endl;
                    f_s_l = forward;
                    quantity = 1;
                }
            } else {
                if (quantity != 0 && cur_unitig_left != NO_NEIGH)
                {
                    g.add_read(cur_unitig_left, quantity);
                    size_t cur_rev_comp = (cur_unitig_left >= total_unitigs / 2) ?
                                          cur_unitig_left - total_unitigs / 2 : cur_unitig_left + total_unitigs / 2;
                    g.add_read(cur_rev_comp, quantity);
                }
                cur_unitig_left = NO_NEIGH;
                quantity = 1;
                f_s_l = true;
            }
        }
        if (quantity != 0 && cur_unitig_left != NO_NEIGH) {
            g.add_read(cur_unitig_left, quantity);
            size_t cur_rev_comp = (cur_unitig_left >= total_unitigs / 2) ?
                                  cur_unitig_left - total_unitigs / 2 : cur_unitig_left + total_unitigs / 2;
            g.add_read(cur_rev_comp, quantity);
        }
    }
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds (FULL) : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
    cout << "Polishing the graph. Wait..."<<endl;
    g.subsane(sequence_map);
}
/*
 * Traversion alternative
 */
void _traverseReadsHash(char * file_left, char * file_right
        ,const unordered_map<Kmer<SPAN>::Type,stored_info> kmer_map
        ,size_t kmerSize, DBG & g, size_t total_unitigs
        ,const vector<string> & sequence_map)
{
    /*
     * Debug options
     */
    std::ofstream pe_information;
    if (Parameters::get().debug)
    {
        pe_information.open("debug/pe_information.txt", std::ofstream::out);
    }
    cout << "Traversing paired-end reads: "<< endl;
    cout << "Left: "<<file_left<<endl;
    cout << "Right: "<<file_right<<endl;
    size_t rate = ceil(float(kmerSize)*float(FIRSTLASTLENGTH)/100);
    cout << "Accepted distances: "<<rate<<endl;
    auto start = chrono::steady_clock::now();
    /*
     * Paired_end input banks
     */
    IBank* inputBank_left = Bank::open (file_left);
    IBank* inputBank_right = Bank::open(file_right);
    /*
     * Number of reads traversed
     */
    size_t number_reads = 0;
    PairedIterator <Sequence> *itPair = new PairedIterator<Sequence>(inputBank_left->iterator(), inputBank_right->iterator());
    ProgressIterator <std::pair<Sequence, Sequence>> progress_iter(itPair, "paired-end", inputBank_left->estimateNbItems());
    /*
     * (TRIAL) Bloom filter for paired-end
     */
    bit_vector pairs_added = bit_vector(pow(2,30),0);
    cout << "Paired-end traversion"<<endl;
    for (progress_iter.first(); !progress_iter.isDone(); progress_iter.next()) {
        Sequence &s1 = itPair->item().first;
        Sequence &s2 = itPair->item().second;
        size_t l1 = s1.getDataSize(), l2 = s2.getDataSize();
        /*
         * Kmer Models
         */
        Kmer<SPAN>::ModelCanonical kmerModel (kmerSize);
        Kmer<SPAN>::ModelCanonical::Iterator kmerItLeft (kmerModel), kmerItRight (kmerModel);
        kmerItLeft.setData(s1.getData());
        kmerItRight.setData(s2.getData());
        /*
         * Traverse kmers
         */
        kmerItRight.first();
        /*
         * (TRIAL) Unordered_set composed unitigs as keys.
         */
        unordered_set<uint32_t> pairs_set_1;
        for (kmerItLeft.first(); !kmerItLeft.isDone(); kmerItLeft.next())
        {
            bool show = false;
            bool forward = true, forward_right = true;
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
                if (map_iterator_right == kmer_map.end())
                {
                    map_iterator_right = (forward)?kmer_map.find(kmerItRight->forward()):kmer_map.find(kmerItRight->revcomp());
                    forward_right = false;
                }
                if (map_iterator_right != kmer_map.end()) {
                    size_t pos_left = (*map_iterator_left).second._pos, u_length = (*map_iterator_left).second.__u_length;
                    size_t unitig_left = (*map_iterator_left).second._unitig
                    , unitig_right = (*map_iterator_right).second._unitig;
                    size_t remaining_unitig_left = (*map_iterator_left).second._remaining
                    , remaining_unitig_right = (*map_iterator_right).second._remaining;
                    size_t rev_comp_unitig_left = (unitig_left >= total_unitigs / 2) ? unitig_left - total_unitigs / 2 :
                                                  unitig_left + total_unitigs / 2,
                            rev_comp_unitig_right = (unitig_right >= total_unitigs / 2) ? unitig_right -
                                                                                          total_unitigs / 2 :
                                                    unitig_right + total_unitigs / 2;
                    if (pos_left >= rate || pos_left >= (u_length - rate))
                    {
                        uint32_t key, rc_key;
                        if (forward) {
                            if (forward_right) {
                                key = (unitig_left << 15)| unitig_right
                                        , rc_key = (rev_comp_unitig_right << 15) | rev_comp_unitig_left;
                                if (pairs_added[key] == 1 & pairs_set_1.find(key) == pairs_set_1.end())
                                {
                                    g.addPair(unitig_left, unitig_right, forward, forward_right);
                                    g.addPair(rev_comp_unitig_right, rev_comp_unitig_left, forward, forward_right);
                                }
                            } else {
                                key = (unitig_left << 15)| rev_comp_unitig_right
                                        , rc_key = (unitig_right << 15) | rev_comp_unitig_left;
                                if (pairs_added[key] == 1 & pairs_set_1.find(key) == pairs_set_1.end())
                                {
                                    g.addPair(unitig_left, rev_comp_unitig_right, forward, forward_right);
                                    g.addPair(unitig_right, rev_comp_unitig_left, forward, forward_right);
                                }
                            }
                        } else {
                            if (forward_right) {
                                key = (rev_comp_unitig_left << 15)| unitig_right
                                        , rc_key = (rev_comp_unitig_right << 15) | unitig_left;
                                if (pairs_added[key] == 1 & pairs_set_1.find(key) == pairs_set_1.end())
                                {
                                    g.addPair(rev_comp_unitig_left, unitig_right, forward, forward_right);
                                    g.addPair(rev_comp_unitig_right, unitig_left, forward, forward_right);
                                }
                            } else {
                                key = (rev_comp_unitig_left << 15)| rev_comp_unitig_right
                                        , rc_key = (unitig_right << 15) | unitig_left;
                                if (pairs_added[key] == 1 & pairs_set_1.find(key) == pairs_set_1.end())
                                {
                                    g.addPair(rev_comp_unitig_left, rev_comp_unitig_right, forward, forward_right);
                                    g.addPair(unitig_right, unitig_left, forward, forward_right);
                                }
                            }
                        }
                        if (pairs_added[key] == 0)
                        {
                            pairs_added[key] = 1;pairs_added[rc_key] = 1;
                            pairs_set_1.emplace(key);pairs_set_1.emplace(rc_key);
                        }
                    }
                }
            }else{
                /*cout << "Kmer no localizado: "<<kmerModel.toString(kmerItLeft->forward())
                <<" Read: "<<number_reads<<" Pos: "<<(s1.getDataSize()-l1)<<endl;*/
            }
            l1--;l2--;
            kmerItRight.next();
            if (kmerItRight.isDone())
                break;
        }
        /*
         * Number of reads
         */
        number_reads++;
    }
    cout << "Number of reads traversed: "<<number_reads<<endl;
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds (Paired-end) : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
    cout << "Post processing pairs. It might take a while..."<<endl;
    start = chrono::steady_clock::now();
    g.post_process_pairs(sequence_map);
    end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds (Paired-end) : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
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
        , char * write_path, const vector<string> & sequence_map
        , size_t num_unitigs_fw,size_t kmer_size, bool force_write = false)
{
    cout << "Writing unitigs"<<endl;
    size_t num_unitigs = 0;
    ofstream outputFile(write_path);
    for (auto unitig:unitigs)
    {
        //cout << "Unitig: "<<endl;
        string full_unitig = "";
        for (size_t i = 0; i < unitig.size();++i)
        {
            //cout << unitig[i]<<" ";
            size_t u = unitig[i];
            bool reverse = false;
            if (u >= (num_unitigs_fw/2))
            {
                u -= (num_unitigs_fw/2);
                reverse = true;
            }
            string seq = sequence_map.at(u);
            seq = ((reverse)?Sequence(seq.c_str()).getRevcomp():seq);
            if (i == 0 && i != unitig.size() -1)
                seq = seq;
            else if (i != unitig.size()-1)
                seq = seq.substr(kmer_size-1, seq.size()-kmer_size+1);
            else if (i == unitig.size()-1 && i != 0)
                seq = seq.substr(kmer_size-1, seq.size());
            /*char * cstr = new char [seq.length()+1];
            std::strcpy (cstr, seq.c_str());*/
            full_unitig += seq;
            //cout <<" Seq "<<seq<<" "<<full_unitig<<" ";
        }
        //cout << endl;
        if (full_unitig.length() > MIN_LENGTH || force_write)
        {
            outputFile << ">"<<num_unitigs++<<endl;
            outputFile << full_unitig;
            outputFile << endl;
        }
    }
    cout << "End!"<< endl;
    outputFile.close();
}


void _build_process_cliques(DBG & g, const vector<string> & sequence_map,
        char * write_path, size_t kmer_size, size_t num_unitigs)
{
    cout << "Processing cliques from DBG... This will take a while!" <<endl;
    DBG apdbg;
    g.build_process_cliques(apdbg, sequence_map, num_unitigs);
    vector<vector<size_t>> unitigs = apdbg.export_unitigs(sequence_map, false);
    vector<vector<size_t>> unitigs_2 = apdbg.export_unitigs_basic(sequence_map, false);
    /*cout << "Unitigs!"<<endl;
    for (auto u: unitigs) {
        for (auto u2: u)
            cout << " " << u2;
        cout << endl;
    }*/
    _write_unitigs(unitigs, write_path, sequence_map, g.vertices(), kmer_size);
}

DBG _buildGraph(char * file)
{
    cout << "Build DBG graph from file"<<endl;
    DBG graph(file);
    return graph;
}

int main (int argc, char* argv[])
{
    cout << "Params: unitigFile, dolarsFile (placement),tmp_file, kmerSize, graphFile, unitigsFasta, unitigsfile "<<endl;
    char * unitigs = argv[1], * dolars = argv[2]
            , * graphFile = argv[5], * unitigsFa = argv[6], * unitigs_file = argv[7], * append_file = argv[8];
    string file_1 =  string(argv[3]) + "/0.fasta", file_2 = string(argv[3]) + "/1.fasta";
    char * file1 = file_1.c_str(), * file2 = file_2.c_str();
    size_t kmerSize = atoi(argv[4]);
    /*
     * Naive!!!
     */
    for (size_t i = 0; i < argc; i++)
    {
        if (strcmp(argv[i],"--debug") == 0)
        {
            cout << "############## Debug mode enabled #################"<<endl;
            Parameters::get().debug = true;
        }
    }
    Parameters::get().kmerSize = kmerSize;
    FMIndex fmIndex;
    unordered_map<Kmer<SPAN>::Type,stored_info> kmer_map;
    bit_vector bDolars;
    RankOnes bDolarRank;
    SelectOnes bDolarSelect;
    /*
     * Trial
     */
    vector<string> sequence_map;
    //unordered_map<size_t, string> sequence_map;
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
        cout << "Creating hash from: "<<unitigsFa<<endl;
        _buildHash(kmer_map, unitigsFa, g.vertices(), kmerSize, sequence_map, g);
        _add_frequencies(kmer_map, append_file, g.vertices(), kmerSize, sequence_map, g);
        _traverseReadsHash(file1, file2, kmer_map, kmerSize, g, g.vertices(), sequence_map);
        cout << "Sanity check!"<<endl;
        vector<vector<size_t>> unitigs = g.export_unitigs(sequence_map);
        _write_unitigs(unitigs, "tmp/testing_unitigs.fa", sequence_map, g.vertices(),kmerSize, false);
        if (Parameters::get().debug)
            g.stats();
        /*cout << "Grafo original"<<endl;
        g.print();*/
        _build_process_cliques(g, sequence_map, unitigs_file, kmerSize,g.vertices());
    }
    //g.print();
}
