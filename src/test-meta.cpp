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

#define STATIC_RESERVE 1024*1024*1024
#define HAND_THRESHOLD 20 //Previo 30

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
        sequence_map.push_back(seq.toString());
        if (!g.getNodeState(num_unitig))
        {
            num_unitig++;
            continue;
        }
        for (kmerIt.first(); !kmerIt.isDone();kmerIt.next())
        {
            if (pos == l-kmerSize)
            {
                kmer_map[kmerIt->forward()] = stored_info(pos, l - pos - 1, num_unitig, l);
                kmer_map[kmerIt->revcomp()] = stored_info(l - pos - 1, pos, total_unitigs / 2 + num_unitig,l);
                pos++;
            } else if (pos == 0) {
                if (kmer_map.find(kmerIt->forward()) == kmer_map.end())
                {
                    kmer_map[kmerIt->forward()] = stored_info(pos, l - pos - 1, num_unitig, l);
                    kmer_map[kmerIt->revcomp()] = stored_info(l - pos - 1, pos, total_unitigs / 2 + num_unitig,l);
                    pos++;
                }
            } else {
                kmer_map[kmerIt->forward()] = stored_info(pos, l - pos - 1, num_unitig,l);
                kmer_map[kmerIt->revcomp()] = stored_info(l - pos - 1, pos, total_unitigs / 2 + num_unitig,l);
                pos++;
            }
        }
        num_unitig++;
    }
    cout <<"Number of unitigs: "<<num_unitig<<"\nNumber of active nodes: "<<g.vertices(true)<<"\nNumber of active edges: "<<g.edges()<<endl;
    cout <<"Number of k-mers: "<<kmer_map.size()<<endl;
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds (FULL) : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
}
/*
 * Add frequencies
 */
bit_vector _add_frequencies(const unordered_map<Kmer<SPAN>::Type,stored_info> & kmer_map,
                char * append_file, size_t total_unitigs, size_t kmerSize
        ,vector<string> & sequence_map, DBG & g)
{
    /*
     * Debug options
     */
    std::ofstream reads_transitions;
    if (Parameters::get().debug)
    {
        reads_transitions.open("debug/reads_transitions.txt", std::ofstream::out);
    }
    cout << "Traversing merged file: "<<append_file<< endl;
    auto start = chrono::steady_clock::now();
    /*
     * Sequence Iterator
     */
    IBank* inputBank = Bank::open (append_file);
    size_t number_of_reads = (*inputBank).getSize(), totalSize, maxSize;
    inputBank->estimate (number_of_reads,totalSize, maxSize);
    ProgressIterator<Sequence> it (*inputBank, "Iterating sequences");
    number_of_reads = ceil(1.1*number_of_reads);
    cout << "Number of reads: "<<number_of_reads<<endl;
    cout << "Number of k-mers under processment: "<<kmer_map.size()<<endl;
    /*
     * MOD! - bitvector de 1's originalmente 0's
     */
    bit_vector reads_with_transition(number_of_reads,0);
    size_t n_read = 0, cuenta = 0;
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
        //cout <<"Sequence: "<<seq.toString()<<" "<<n_read<<endl;
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
                size_t rev_comp_unitig_left = Common::rev_comp_index(total_unitigs,unitig_left);
                if (cur_unitig_left != NO_NEIGH) {
                    size_t difference = abs((int)unitig_left-(int)cur_unitig_left);
                    if (difference != 0 && difference != (total_unitigs/2)) {
                        if (!reads_with_transition[n_read]) {
                            if (Parameters::get().debug)
                                reads_transitions << "Read included: "<<n_read<<endl;
                            reads_with_transition[n_read] = 1;
                        }
                        size_t cur_rev_comp = Common::rev_comp_index(total_unitigs, cur_unitig_left);
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
                    if (g.out_degree(unitig_left) == 1)
                        reads_with_transition[n_read] = 1;
                    cur_unitig_left = unitig_left;
                    f_s_l = forward;
                    quantity = 1;
                }
            } else {
                cur_unitig_left = NO_NEIGH;
                quantity = 1;
                f_s_l = true;
            }
        }
        n_read++;
    }
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in milliseconds (FULL) : "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
    cout << "Polishing the graph. Wait..."<<endl;
    if (Parameters::get().debug)
        g.print(INF, INF, "graphs/dbg.txt");
    g.subsane(sequence_map);
    g.polish();
    cout << "Correcting kmer map:"<<endl;
    vector<Kmer<SPAN>::Type> kmers_to_remove;
    for (auto k_v:kmer_map) {
        if (!g.getNodeState(k_v.second._unitig))
            kmers_to_remove.push_back(k_v.first);
    }
    for (auto k: kmers_to_remove)
        kmer_map.erase(k);
    cout <<"**After polishing**\nNumber of active nodes: "<<g.vertices(true)<<"\nNumber of active edges: "<<g.edges()<<"\nNumber of kmers: "<<kmer_map.size()<<endl;
    rank_support_v5<> rb(&reads_with_transition);
    cout << "Number ofreads with transitions: "<<rb(number_of_reads)<<" out of "<<((float)rb(n_read) / (float)n_read)*100<<"%"<<endl;
    return reads_with_transition;
}
/*
 * Traversion of paired-end reads
 */
void _traverseReadsHash(char * file_left, char * file_right
        ,const unordered_map<Kmer<SPAN>::Type,stored_info> kmer_map
        ,size_t kmerSize, DBG & g, size_t total_unitigs
        ,const vector<string> & sequence_map
        ,const bit_vector reads_with_transitions)
{
    /*
     * Debug options
     */
    std::ofstream pe_information, histogram_file;
    if (Parameters::get().debug)
    {
        pe_information.open("debug/pe_information.txt", std::ofstream::out);
        histogram_file.open("debug/histogram_pe.txt", std::ofstream::out);
    }
    cout << "Traversing paired-end reads: "<< endl;
    cout << "Left: "<<file_left<<endl;
    cout << "Right: "<<file_right<<endl;
    size_t rate = ceil(float(kmerSize)*float(FIRSTLASTLENGTH)/100);
    /*
     * Testing
     */
    rate = 1;
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
    size_t number_reads = 0, cnr = 0;
    PairedIterator <Sequence> *itPair = new PairedIterator<Sequence>(inputBank_left->iterator(), inputBank_right->iterator());
    ProgressIterator <std::pair<Sequence, Sequence>> progress_iter(itPair, "paired-end", inputBank_left->estimateNbItems());
    /*
     * (TRIAL) Bloom filter for paired-end
     */
    auto preCounters = new std::atomic<unsigned char>[STATIC_RESERVE];
    for (size_t i = 0; i < STATIC_RESERVE; ++i)
        preCounters[i] = 0;

    cout << "Paired-end traversion"<<endl;
    for (progress_iter.first(); !progress_iter.isDone(); progress_iter.next())
    {
        size_t number_of_kmer = 0;
        if (Parameters::get().t_data == "virus")
        {
            if (!reads_with_transitions[cnr++])
                continue;
        }
        if (Parameters::get().debug)
        {
            pe_information << "Read traversed: "<<(cnr-1)<<endl;
        }
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
        unordered_set<size_t> pairs_set_1;
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
                    if (!g.getNodeState(unitig_left) || !g.getNodeState(unitig_right))
                        cout << "Impossible!"<<endl;
                    if (pos_left >= rate || pos_left >= (u_length - rate))
                    {
                        uint32_t key, rc_key;
                        if (forward) {
                            if (forward_right) {
                                key = (unitig_left << 15)| unitig_right
                                        , rc_key = (rev_comp_unitig_right << 15) | rev_comp_unitig_left;
                                if (pairs_set_1.find(key) == pairs_set_1.end())
                                {
                                    if (preCounters[key] < 255)
                                    {
                                        preCounters[key]++;preCounters[rc_key]++;
                                        pairs_set_1.emplace(key);pairs_set_1.emplace(rc_key);
                                        if (preCounters[key] >= HAND_THRESHOLD) {
                                            g.addPair(unitig_left, unitig_right, forward, forward_right);
                                            g.addPair(rev_comp_unitig_right, rev_comp_unitig_left, forward,forward_right);
                                        }
                                    }
                                }
                            } else {
                                key = (unitig_left << 15)| rev_comp_unitig_right
                                        , rc_key = (unitig_right << 15) | rev_comp_unitig_left;
                                if (pairs_set_1.find(key) == pairs_set_1.end())
                                {
                                    if (preCounters[key] < 255) {
                                        preCounters[key]++;preCounters[rc_key]++;
                                        pairs_set_1.emplace(key);pairs_set_1.emplace(rc_key);
                                        if (preCounters[key] >= HAND_THRESHOLD) {
                                            g.addPair(unitig_left, rev_comp_unitig_right, forward, forward_right);
                                            g.addPair(unitig_right, rev_comp_unitig_left, forward, forward_right);
                                        }
                                    }
                                }
                            }
                        } else {
                            if (forward_right) {
                                key = (rev_comp_unitig_left << 15)| unitig_right
                                        , rc_key = (rev_comp_unitig_right << 15) | unitig_left;
                                if (pairs_set_1.find(key) == pairs_set_1.end())
                                {
                                    if (preCounters[key] < 255) {
                                        preCounters[key]++;preCounters[rc_key]++;
                                        pairs_set_1.emplace(key);pairs_set_1.emplace(rc_key);
                                        if (preCounters[key] >= HAND_THRESHOLD) {
                                            g.addPair(rev_comp_unitig_left, unitig_right, forward, forward_right);
                                            g.addPair(rev_comp_unitig_right, unitig_left, forward, forward_right);
                                        }
                                    }
                                }
                            } else {
                                key = (rev_comp_unitig_left << 15)| rev_comp_unitig_right
                                        , rc_key = (unitig_right << 15) | unitig_left;
                                    if (pairs_set_1.find(key) == pairs_set_1.end())
                                    {
                                        if (preCounters[key] < 255) {
                                            preCounters[key]++;preCounters[rc_key]++;
                                            pairs_set_1.emplace(key);pairs_set_1.emplace(rc_key);
                                            if (preCounters[key] >= HAND_THRESHOLD) {
                                                g.addPair(rev_comp_unitig_left, rev_comp_unitig_right, forward,
                                                          forward_right);
                                                g.addPair(unitig_right, unitig_left, forward, forward_right);
                                            }
                                        }
                                }
                            }
                        }
                        if (preCounters[key] == 0)
                        {
                            preCounters[key] = 1;preCounters[rc_key] = 1;
                            pairs_set_1.emplace(key);pairs_set_1.emplace(rc_key);
                        }
                    } 
                }
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
    if (Parameters::get().debug)
    {
        std::ofstream outfile("graphs/pair_end_info.txt", std::ofstream::binary);
        for (size_t i = 0; i < STATIC_RESERVE;++i)
        {
            if (preCounters[i] >= 1)
            {
                size_t node_left = (i >> 15), node_right = (i & 32767);
                //auto neighbors = g.getNeighbor(node_left);
                    outfile <<" Left: "<<node_left<<" - "<<node_right<<" Freq: "<<((int)preCounters[i])<<endl;
                /*for (auto n: neighbors)
                {
                    size_t key = (n << 15) | node_right;
                    if ((preCounters[key] > 0) & (preCounters[i] < HAND_THRESHOLD)){
                        int new_count = ((int)preCounters[i]) + ((int)preCounters[key]);
                        outfile << "Left: "<<node_left<<" - "<<node_right<<" Freq: "<<((int)preCounters[i])<<" + "<<((int)preCounters[key])<<" = "<<new_count<<endl;
                        if (new_count > HAND_THRESHOLD)
                            g.addPair(node_left, node_right, true, true);
                    } else {
                        outfile <<"Neigh: "<<n<<" Left: "<<node_left<<" - "<<node_right<<" Freq: "<<((int)preCounters[i])<<endl;
                    }
                }*/
                histogram_file<<node_left<<","<<node_right<<","<<((int)preCounters[i])<<endl;
            }
        }
        histogram_file.close();
        outfile.close();
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
    if (Parameters::get().debug) {
        g.print(INF, INF, "graphs/dbg_postsubsane.txt");
        g.export_to_gfa(sequence_map, "graphs/dbg_postsubsane.gfa");
        cout << "Graphs exported as 'txt' and 'gfa'"<<endl;
    }
}

/*
 * Traverse Paired-end reads
 */
void _traverseReads(char * file_left, char * file_right ,const FMIndex & fm, const RankOnes & rank,
        const SelectOnes & select, size_t kmerSize, DBG & g)
{
    cout << "Paired-end traversal "<<file_left<<" "<<file_right<<endl;
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
static vector<size_t> DEFAULT_VECTOR;
vector<bool> _write_unitigs(vector<vector<size_t>> unitigs
        , string write_path, const vector<string> & sequence_map
        , size_t num_unitigs_fw,size_t kmer_size, bool force_write = false, vector<size_t> & flows = DEFAULT_VECTOR)
{
    vector<bool> wrote = vector<bool>(unitigs.size(), false);
    cout << "Writing unitigs: "<<unitigs.size()<<endl;
    bool write_flows = (flows.size() > 0);
    size_t num_unitigs = 0;
    ofstream outputFile(write_path);
    for (size_t j = 0; j < unitigs.size(); ++j)
    {
        vector<size_t> unitig = unitigs[j];
        string full_unitig = "";
        for (size_t i = 0; i < unitig.size();++i)
        {
            //cout << unitig[i]<<" ";
            size_t u = unitig[i];
            string seq = Common::return_unitig(sequence_map, u);
            if (i != 0)
                seq = seq.substr(kmer_size-1, seq.size() - kmer_size + 1);
            /*char * cstr = new char [seq.length()+1];
            std::strcpy (cstr, seq.c_str());*/
            full_unitig += seq;
            //cout <<" Seq "<<seq<<" "<<full_unitig<<" ";
        }
        //cout << endl;
        if (full_unitig.length() > MIN_LENGTH || force_write)
        {
            wrote[j] = true;
            if (write_flows)
                outputFile << ">"<<num_unitigs++<<"-"<<full_unitig.length()<<"-"<<flows[j]<<endl;
            else
                outputFile << ">"<<num_unitigs++<<endl;
            outputFile << full_unitig;
            outputFile << endl;
        }
    }
    cout << "End!"<< endl;
    outputFile.close();
    return wrote;
}

void _write_freqs(vector<size_t> unitig_flow, string write_path, vector<bool> index)
{
    cout << "Writing flow: "<<unitig_flow.size()<<endl;
    ofstream outputFile(write_path);
    outputFile << "Freqs"<<endl;
    for (size_t i = 0; i < unitig_flow.size(); ++i)
    {
        size_t flow = unitig_flow[i];
        if (index[i])
            outputFile << flow<<endl;
    }
    outputFile.close();
}
void _build_process_cliques(DBG & g, const vector<string> & sequence_map,
        char * write_path, size_t kmer_size, size_t num_unitigs)
{
    cout << "Processing cliques from DBG... This will take a while!" <<endl;
    DBG apdbg;
    g.build_process_cliques(apdbg, sequence_map, num_unitigs);
    if (Parameters::get().debug) {
        cout << "Exporting ADPBG" << endl;
        apdbg.print(1, INF);
        apdbg.export_to_gfa(sequence_map);
    }
    /*
     * Max flow
     */
    /*cout << "MCP via max flow"<<endl;
    float pre_flow = apdbg.to_max_flow_solution();
    priority_queue<pair<size_t,vector<OwnNode_t>>> unitigs_nf_with_freqs = apdbg.get_min_cost_flow_paths(pre_flow);
    vector<vector<OwnNode_t>> unitigs_nf;
    vector<size_t> flows;
    while(!unitigs_nf_with_freqs.empty()) {
        pair<size_t,vector<OwnNode_t>> u = unitigs_nf_with_freqs.top();
        unitigs_nf_with_freqs.pop();
        unitigs_nf.push_back(u.second);
        flows.push_back(u.first);
    }
    string str_obj_nf("tmp/unitigs-viaDBG-nf.fasta");
    char * write_path_nf = &str_obj_nf[0];
    if (unitigs_nf.size() > 0) {
        vector<bool> index_unitigs = _write_unitigs(unitigs_nf, write_path_nf, sequence_map, g.vertices(), kmer_size, false, flows);
        _write_freqs(flows, "tmp/flows.csv", index_unitigs);
    }*/
    /*
     * Testing - standard approach
     */
    cout << "MCP standard approach"<<endl;
    priority_queue<pair<size_t,vector<OwnNode_t>>> unitigs_nf_with_freqs_2 = apdbg.solve_std_mcp(sequence_map);
    cout << "End std_mcp"<<endl;
    vector<vector<OwnNode_t>> unitigs_nf_2;
    vector<size_t> flows_2;
    while(!unitigs_nf_with_freqs_2.empty()) {
        pair<size_t,vector<OwnNode_t>> u = unitigs_nf_with_freqs_2.top();
        unitigs_nf_with_freqs_2.pop();
        unitigs_nf_2.push_back(u.second);
        flows_2.push_back(u.first);
        cout << "Flow: "<<u.first<<endl;
        for (auto u:u.second)
            cout << " " << u <<" -> ";
        cout << endl;
    }
    string str_obj_nf_2("tmp/unitigs-viaDBG-nf-std.fasta");
    char * write_path_nf_2 = &str_obj_nf_2[0];
    if (unitigs_nf_2.size() > 0) {
        vector<bool> index_unitigs = _write_unitigs(unitigs_nf_2, write_path_nf_2, sequence_map, g.vertices(), kmer_size, false, flows_2);
        _write_freqs(flows_2, "tmp/flows_2.csv", index_unitigs);
    }
    //vector<vector<size_t>> unitigs = apdbg.export_unitigs(sequence_map, false);
    //vector<vector<size_t>> unitigs_2 = apdbg.export_unitigs_basic(sequence_map, false);
    /*cout << "Unitigs!"<<endl;
    for (auto u: unitigs) {
        for (auto u2: u)
            cout << " " << u2;
        cout << endl;
    }*/
    //_write_unitigs(unitigs, write_path, sequence_map, g.vertices(), kmer_size);
    string str_obj("tmp/unitigs-viaDBG-greedy.fasta");
    char * write_path_extra = &str_obj[0];
    //_write_unitigs(unitigs_2, write_path_extra,  sequence_map, g.vertices(), kmer_size);
}

DBG _buildGraph(char * file)
{
    cout << "Build DBG graph from file"<<endl;
    DBG graph(file);
    return graph;
}

int main (int argc, char* argv[])
{
    cout << "Params: unitigFile, dolarsFile (placement),tmp_file, kmerSize, graphFile, unitigsFasta, unitigsfile, appendfile, paired_end (optional) "<<argc<<endl;
    char * unitigs = argv[1], * dolars = argv[2]
            , * graphFile = argv[5], * unitigsFa = argv[6], * unitigs_file = argv[7], * append_file = argv[8];
    string paired_end = "";
    if (argc == 12)
    {
        paired_end = string(argv[9]);
    }
    //Ajustar para pear
    string file_1 = string(""), file_2 = string("");
    if (paired_end == ""){
        file_1 =  string(argv[3]) + "/0.fasta";
        file_2 = string(argv[3]) + "/1.fasta";
    } else {
        auto files = OwnSystem::get_directories(paired_end);
        file_1 = files[0];
        file_2 = files[1];
    }
    cout << "File Left: "<<file_1<<" File Right: "<<file_2<<endl;
    char * file1 = &file_1[0], * file2 = &file_2[0];
    // Parameters read
    Parameters::get().check_cmd_line(argc, argv);
    Parameters::get().print_info();
    size_t kmerSize = Parameters::get().kmerSize;
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
    /*
     * Polishing the graph
     */
    g.polish();
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
        bit_vector reads_with_transitions = _add_frequencies(kmer_map, append_file, g.vertices(), kmerSize, sequence_map, g);
        _traverseReadsHash(file1, file2, kmer_map, kmerSize, g, g.vertices(), sequence_map, reads_with_transitions);
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
