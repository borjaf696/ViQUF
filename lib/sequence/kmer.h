#ifndef KMER_H
#define KMER_H
    #include <unordered_map>
    #include <memory>
    #include <assert.h>
    #include "sequence_container.h"

    using namespace std;
    /*
     * Single_end Kmer
     */
    class Kmer
    {
    public:
        bool exist = true;
        Kmer():exist(false){}
        Kmer(const DnaSequence&,size_t,size_t);
        Kmer(DnaSequence seq):_seq(seq){}
        Kmer(const Kmer & kmer){
            _seq = kmer.getSeq();
            _standard = kmer._standard;
        }
        Kmer(const Kmer *kmer){
            _seq = kmer->getSeq();
            _standard = kmer->_standard;
        }

        Kmer(const string &string_own)
        {
            _seq = DnaSequence(string_own);
        }

        Kmer rc(){
            return Kmer(_seq.complement());
        }

        void appendRight(DnaSequence::NuclType);
        void appendLeft(DnaSequence::NuclType);
        void appendRightReplace(DnaSequence::NuclType);
        void appendLeftReplace(DnaSequence::NuclType);

        const DnaSequence getSeq() const {
            return _seq;
        }

        const DnaSequence* getSeq_ref() const{
            return &_seq;
        }

        void standard(){
            if (!_standard)
            {
                Kmer rc = this->rc();
                if ((*this) > rc) {
                    (*this) = rc;
                    _standard = true;
                }
                _standard = true;
            }
        }

        DnaSequence::NuclType at(size_t index) const;

        bool operator>(const Kmer&) const;
        bool operator<(const Kmer&) const;
        bool operator==(const Kmer&) const;
        bool operator!=(const Kmer&) const;
        Kmer& operator=(const Kmer&);

        std::size_t hash() const{
            return _seq.hash();
        }

        std::string str() const
        {
            return _seq.str();
        }

        DnaSequence substr(size_t start, size_t end) const
        {
            return _seq.substr(start,end);
        }

        vector<DnaSequence> firstLastSubstr(size_t length, bool canonical) const
        {
            return _seq.firstLastSubstr(length, canonical);
        }

        pair<Kmer, Kmer> preffixsuffix()
        {
            pair<DnaSequence,DnaSequence> pre_suf =_seq.preffixsuffix();
            return pair<Kmer,Kmer>(Kmer(pre_suf.first),Kmer(pre_suf.second));
        };
        /*
         * Odd but need :(
         */
        size_t length() const
        {
            return _seq.length();
        }
    private:
        DnaSequence _seq;
        bool _standard = false;
    };

    /*
     * Pair_End Kmers
     */
    class Pair_Kmer
    {
    public:
        bool exist = true;
        Pair_Kmer():exist(false){}
        Pair_Kmer(const DnaSequence&,const DnaSequence&,size_t,size_t);
        Pair_Kmer(DnaSequence seq_left, DnaSequence seq_right):
                _seq_left(seq_left),_seq_right(seq_right){}
        Pair_Kmer(const Pair_Kmer & kmer){
            pair<DnaSequence,DnaSequence> seqs = kmer.getSeq();
            _seq_left = seqs.first;
            _seq_right = seqs.second;
            _standard = kmer._standard;
        }
        Pair_Kmer(const Pair_Kmer *kmer){
            pair<DnaSequence,DnaSequence> seqs = kmer->getSeq();
            _seq_left = seqs.first;
            _seq_right = seqs.second;
            _standard = kmer->_standard;
        }

        Pair_Kmer(const string &string_own_left, const string &string_own_right)
        {
            _seq_left = DnaSequence(string_own_left);
            _seq_right = DnaSequence(string_own_right);
        }

        Pair_Kmer rc(){
            return Pair_Kmer(_seq_left.complement(),_seq_right.complement());
        }

        void standard()
        {
            if (!_standard)
            {
                if (_seq_left > _seq_left.complement())
                    _seq_left = _seq_left.complement();
                if (_seq_right > _seq_right.complement())
                    _seq_right = _seq_right.complement();
                _standard = true;
            }
        }

        void appendRight(DnaSequence::NuclType, DnaSequence::NuclType);
        void appendLeft(DnaSequence::NuclType, DnaSequence::NuclType);

        pair<DnaSequence,DnaSequence> getSeq() const
        {
            return pair<DnaSequence,DnaSequence>(_seq_left,_seq_right);
        }

        pair<const DnaSequence*,const DnaSequence*> getSeq_ref() const
        {
            return pair<const DnaSequence*,const DnaSequence*>(&_seq_left,&_seq_right);
        }

        pair<Kmer,Kmer> getKmers()
        {
            return pair<Kmer,Kmer>{Kmer(_seq_left),Kmer(_seq_right)};
        }

        pair<DnaSequence::NuclType,DnaSequence::NuclType> at(size_t index) const;

        bool operator>(const Pair_Kmer&) const;
        bool operator<(const Pair_Kmer&) const;
        bool operator==(const Pair_Kmer&) const;
        bool operator!=(const Pair_Kmer&) const;
        Pair_Kmer& operator=(const Pair_Kmer&);

        std::size_t hash() const
        {
            size_t seed = 0;
            boost::hash_combine(seed,_seq_left.hash());
            boost::hash_combine(seed,_seq_right.hash());
            return seed;
        }

        pair<string,string> str() const
        {
            return pair<string,string>(_seq_left.str(),_seq_right.str());
        }

        pair<DnaSequence,DnaSequence> substr(size_t start, size_t end) const
        {
            return pair<DnaSequence,DnaSequence>(_seq_left.substr(start,end),_seq_right.substr(start,end));
        }
    private:
        DnaSequence _seq_left,_seq_right;
        bool _standard = false;
    };
    namespace std{
        template <>
        struct hash<Kmer>{
            size_t operator()(const Kmer& kmer) const{
                return kmer.hash();
            }
        };
        template<>
        struct hash<Pair_Kmer>{
            size_t operator()(const Pair_Kmer& pair_kmer) const{
                return pair_kmer.hash();
            }
        };
        template<>
        struct hash<DnaSequence>{
            size_t operator()(const DnaSequence& seq) const {
                return seq.hash();
            }
        };
    };
    /*
     * Complementos Kmer
     */
    template <bool> struct KmerInfo;
    template<> struct KmerInfo<false>{
        KmerInfo()
        {}
        KmerInfo(Kmer kmer1, size_t pos)
                :kmer(kmer1),kmer_pos(pos)
        {}
        KmerInfo(const KmerInfo & k_info)
                :kmer(k_info.kmer),kmer_pos(k_info.kmer_pos)
        {
        }
        KmerInfo& operator=(const KmerInfo& k_info)
        {
            kmer = k_info.kmer;
            kmer_pos = k_info.kmer_pos;
            return *this;
        }
        size_t hash() const
        {
            size_t seed = 0;
            boost::hash_combine(seed,kmer.hash());
            boost::hash_combine(seed, kmer_pos);
            return seed;
        }
        bool operator==(const KmerInfo & k_info) const
        {
            return ((kmer == k_info.kmer) && (kmer_pos == k_info.kmer_pos));
        }
        bool operator!=(const KmerInfo & k_info) const
        {
            return !(*this == k_info);
        }
        string str() const {
            return " KmerInfo<false>: "+kmer.str()+";"+to_string(kmer_pos);
        }
        Kmer kmer;
        size_t kmer_pos;
    };

    template<> struct KmerInfo<true>{
        KmerInfo()
        {}
        KmerInfo(Pair_Kmer kmer, size_t pos)
                :pair_kmer(kmer),kmer_pos(pos)
        {}
        size_t hash() const
        {
            size_t seed = 0;
            boost::hash_combine(seed, pair_kmer.hash());
            boost::hash_combine(seed, dist);
            return seed;
        }
        KmerInfo(const KmerInfo & k_info)
                :pair_kmer(k_info.pair_kmer),kmer_pos(k_info.kmer_pos),dist(k_info.dist)
        {
        }
        KmerInfo& operator=(const KmerInfo& k_info)
        {
            pair_kmer = k_info.pair_kmer;
            kmer_pos = k_info.kmer_pos;
            dist = k_info.dist;
            return *this;
        }
        bool operator==(const KmerInfo & k_info) const
        {
            return (pair_kmer == k_info.pair_kmer) && (kmer_pos == k_info.kmer_pos) && (dist == k_info.dist);
        }
        bool operator!=(const KmerInfo & k_info) const
        {
            return !((pair_kmer == k_info.pair_kmer) && (kmer_pos == k_info.kmer_pos) && (dist == k_info.dist));
        }
        string str()
        {
            pair<string,string> str_pair = pair_kmer.str();
            return "KmerInfo<true>: "+str_pair.first+"|"+str_pair.second+";"+to_string(kmer_pos);
        }
        Pair_Kmer pair_kmer;
        size_t kmer_pos, dist = 0;
    };

    namespace std
    {
        template <bool P>
        struct hash<KmerInfo<P>>
        {
            size_t operator()(const KmerInfo<P> & k_info) const
            {
                return k_info.hash();
            }
        };
    }
    template<bool P> class KmerIt;
    /*
     * Iterator Single_end
     */
    template<>
    class KmerIt<false>{
    public:
        typedef forward_iterator_tag iterator_category;
        KmerIt(const DnaSequence*,size_t, size_t);

        bool operator==(const KmerIt&) const;
        bool operator!=(const KmerIt&) const;

        //Return the kmer which is being read right now
        KmerInfo<false> operator*() const;
        KmerIt& operator++();
    protected:
        const DnaSequence* _own_seq;
        bool _paired;
        size_t _pos, _kmer_length;
        Kmer _kmer;
    };

    /*
     * Iterator Pair_end
     */
    template<>
    class KmerIt<true>{
    public:
        typedef forward_iterator_tag iterator_category;
        KmerIt(const DnaSequence*, const DnaSequence*,size_t,size_t);

        bool operator==(const KmerIt&) const;
        bool operator!=(const KmerIt&) const;

        //Return the kmer which is being read right now
        KmerInfo<true> operator*() const;
        KmerIt& operator++();
    protected:
        const DnaSequence * _seq_left, * _seq_right;
        bool _paired;
        size_t _pos,  _kmer_length;
        Pair_Kmer _pair_kmer;
    };

    template<bool P> class IterKmers;
    /*
     * Single_end
     */
    template<>
    class IterKmers<false>{
    public:
        IterKmers(const DnaSequence& seq, size_t kmer_size):_seq(seq),_kmer_length(kmer_size){}
        KmerIt<false> begin();
        KmerIt<false> end();

    private:
        const DnaSequence& _seq;
        size_t _kmer_length;
    };

    /*
     * Pair_end
     */
    template<>
    class IterKmers<true>{
    public:
        IterKmers(const DnaSequence &seq_left, const DnaSequence &seq_right, size_t kmer_size):
                _seq_left(seq_left),_seq_right(seq_right), _kmer_length(kmer_size)
        {
            //assert(_seq_right.length() >= _seq_left.length());
        }
        KmerIt<true> begin();
        KmerIt<true> end();

    private:
        const DnaSequence& _seq_left, _seq_right;
        size_t _kmer_length;
    };
#endif