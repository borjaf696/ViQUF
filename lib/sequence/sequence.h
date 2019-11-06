#include <boost/functional/hash.hpp>
#include <cassert>
#include <string>
#include <unistd.h>
#include <vector>
#include <iostream>
#include <algorithm>

class DnaSequence
{
public:
    typedef size_t NuclType;

private:
    /*
     * I should take a look on that to do it automatic
     */
    static const size_t NUCL_BITS = 2;
    static const size_t operator_and = 3;
    static const size_t NUCL_IN_CHUNK = (sizeof(NuclType) * 8) / NUCL_BITS;
    static const size_t offset = NUCL_BITS;

    struct SharedBuffer
    {
        SharedBuffer(): useCount(0), length(0) {}
        size_t useCount;
        size_t length;
        std::vector<size_t> chunks;
    };

public:
    static char nfi(size_t nuc)
    {
        static char table[] = {'A', 'C', 'G', 'T'};
        return table[nuc];
    }

    static size_t getNumPointers()
    {
        return num_pointers;
    }

    DnaSequence():
            _complement(false)
    {
        _data = new SharedBuffer();
        num_pointers++;
        ++_data->useCount;
    }

    ~DnaSequence()
    {
        if (_data != nullptr)
        {
            --_data->useCount;
            if (_data->useCount == 0)
            {
                _data->chunks.clear();
                _data->chunks.shrink_to_fit();
                num_pointers--;
                delete _data;
            }
        }
    }

    explicit DnaSequence(const std::string& string):
            _complement(false)
    {
        size_t i_rep = 0, num_n = 0;
        _data = new SharedBuffer;
        num_pointers++;
        ++_data->useCount;

        if (string.empty()) return;

        _data->length = string.length();
        _data->chunks.assign((_data->length - 1) / NUCL_IN_CHUNK + 1, 0);
        for (size_t i = 0; i < string.length(); ++i) {
            size_t chunkId = i_rep / NUCL_IN_CHUNK;
            size_t nt;
            if (string[i] == 'N') {
                if (i == 0)
                {
                    nt = rand()%4;
                    _data->chunks[chunkId] |= nt << (i_rep % NUCL_IN_CHUNK) * offset;
                }else{
                    num_n++;
                    continue;
                }
            }else
            {
                nt = dnaToId(string[i]);
                _data->chunks[chunkId] |= nt << (i_rep % NUCL_IN_CHUNK) * offset;
            }
            i_rep++;
        }
        _data->length = _data->length-num_n;
    }

    //Only accepts 1 nuc insertion -> Recordar que estan al reves en _data
    void append_nuc_right (DnaSequence::NuclType dnaSymbol) const{
        if ((_data->length==0) && (_data->chunks.size() == 0))
        {
            _data->chunks.push_back(dnaSymbol);
            _data->length++;
            return;
        }
        if ((_data->length / NUCL_IN_CHUNK) == _data->chunks.size())
            _data->chunks.push_back(dnaSymbol);
        else
            _data->chunks[_data->chunks.size()-1] |= (dnaSymbol
                    << ((_data->length % NUCL_IN_CHUNK)*offset));
        _data->length++;
    }

    void append_seq_right (DnaSequence sequence) const {
        for (size_t i = 0; i < sequence.length(); ++i)
            append_nuc_right(sequence.atRaw(i));
    }

    void append_with_replace_right(DnaSequence::NuclType symbol) const {
        for (int i = 0; i < std::max(0,(int)_data->chunks.size()-1); ++i)
            _data->chunks[i] = (_data->chunks[i] >> offset) |
                               (_data->chunks[i+1] << (NUCL_IN_CHUNK-1)*offset);
        _data->chunks[_data->chunks.size()-1] = (_data->chunks[_data->chunks.size()-1]>> offset) |
                                                (symbol << ((_data->length-1) % NUCL_IN_CHUNK)*offset);
    }

    void append_nuc_left (DnaSequence::NuclType dnaSymbol) const{
        _data->length++;
        if ((_data->length/NUCL_IN_CHUNK) == _data->chunks.size())_data->chunks.push_back(0);
        for (uint i = _data->chunks.size()-1; i > 0; --i) {
            _data->chunks[i] = (_data->chunks[i] << offset) | (
                    (_data->chunks[i-1] >> (NUCL_IN_CHUNK-1)*offset));
        }
        _data->chunks[0] = (_data->chunks[0] << offset) | dnaSymbol;
    }

    void append_seq_left(DnaSequence sequence) const {
        for (size_t i = 0; i < sequence.length(); ++i)
            append_nuc_left(sequence.atRaw(i));
    }

    void append_with_replace_left(DnaSequence::NuclType symbol) const{
        for (int i = _data->chunks.size()-1; i > 0; --i){
            _data->chunks[i] = ((_data->chunks[i]<<offset) |
                                (_data->chunks[i-1] >> (NUCL_IN_CHUNK-1)*offset));
        }
        _data->chunks[0] = (_data->chunks[0] << offset) | symbol ;
    }

    void set(DnaSequence::NuclType dnaSymbol, size_t index) const
    {
        size_t aux_v = ~(operator_and << ((index % NUCL_IN_CHUNK)*offset));
        _data->chunks[index / NUCL_IN_CHUNK] &= aux_v;
        _data->chunks[index / NUCL_IN_CHUNK] |= (dnaSymbol << ((index % NUCL_IN_CHUNK)*offset));
    }

    DnaSequence(const DnaSequence& other):
            _complement(other._complement)
    {
        _data = new SharedBuffer();
        num_pointers++;
        _data->length = other.length();
        _data->chunks = other._data->chunks;
        ++_data->useCount;
        /*Original*/
        //		++_data->useCount;
    }

    DnaSequence(DnaSequence&& other):
            _data(other._data),
            _complement(other._complement)
    {
        other._data = nullptr;
    }

    DnaSequence& operator=(const DnaSequence& other)
    {
        --_data->useCount;
        if (_data->useCount == 0){
            --num_pointers;
            _data->chunks.clear();
            _data->chunks.shrink_to_fit();
            delete _data;
        }

        _complement = other._complement;
        _data = other._data;
        ++_data->useCount;
        return *this;
    }

    DnaSequence& operator=(DnaSequence&& other)
    {
        --_data->useCount;
        if (_data->useCount == 0) {
            --num_pointers;
            _data->chunks.clear();
            _data->chunks.shrink_to_fit();
            delete _data;
        }

        _complement = other._complement;
        _data = other._data;
        other._data = nullptr;
        return *this;
    }

    std::vector<size_t> getChunk()const {
        return _data->chunks;
    }

    NuclType operator[](int index){
        if (_complement)
            index = _data->length-index+1;
        size_t id = (_data->chunks[index / NUCL_IN_CHUNK] >> (index % NUCL_IN_CHUNK)*offset)&operator_and;
        if (id > 3)
            return 4;
        return !_complement?id : ~id & 3;
    }

    //TODO: Revisar las constantes esas
    const DnaSequence& operator*() const{
        return *this;
    }

    /*
     * Seria necesario extenderlo para tratar la orientacion de las lecturas
     */
    bool operator==(const DnaSequence& other) const {
        if (other.length() != this->length())
            return false;
        return (this->_data->chunks == other._data->chunks);
    }

    bool operator!=(const DnaSequence& other) const {
        if (other.length() != this->length())
            return true;
        return !(this->_data->chunks == other._data->chunks);
    }

    bool higher(const DnaSequence & other) const
    {
        if (other.length() > length())
            return false;
        if (other.length() < length())
            return true;
        for (uint i = 0; i < _data->length;++i)
        {
            if (atRaw(i) > other.atRaw(i))
                return true;
            if (atRaw(i) < other.atRaw(i))
                return false;
        }
        return false;
    }

    bool lower(const DnaSequence & other) const
    {
        if (other.length() > length())
            return true;
        if (other.length() < length())
            return false;
        for (uint i = 0; i < _data->length;++i)
        {
            if (atRaw(i) > other.atRaw(i))
                return false;
            if (atRaw(i) < other.atRaw(i))
                return true;
        }
        return false;
    }

    bool operator<(const DnaSequence & other) const
    {
        return (lower(other));
    }

    bool operator>(const DnaSequence & other) const
    {
        return (higher(other));
    }

    size_t length() const {return _data->length;}

    size_t hash() const
    {
        size_t seed = 0;
        for (uint i = 0; i < _data->chunks.size();++i)
            boost::hash_combine(seed, _data->chunks[i]);
        return seed;
    }

    char at(size_t index) const
    {
        if (_complement)
        {
            index = _data->length - index - 1;
        }
        if (index >= length())
            index = length()-1;
        size_t id = (_data->chunks[index / NUCL_IN_CHUNK] >>
                                                          (index % NUCL_IN_CHUNK) * offset ) & operator_and;
        if (id > 3)
            return 'N';
        return idToDna(!_complement ? id : ~id & 3);
    }

    NuclType atRaw(size_t index) const
    {
        if (_complement)
        {
            index = _data->length - index - 1;
        }
        size_t id = (_data->chunks[index / NUCL_IN_CHUNK] >>
                                                          (index % NUCL_IN_CHUNK) * offset ) & operator_and;
        if (id > 3)
            return 4;
        return !_complement ? id : ~id & 3;
    }

    DnaSequence complement() const
    {
        DnaSequence complSequence(*this);
        complSequence._complement = !(this->_complement);
        complSequence = complSequence.substr(0, complSequence._data->length);
        complSequence._complement = !(complSequence._complement);
        return complSequence;
    }

    DnaSequence substr(size_t start, size_t length) const;
    std::pair<DnaSequence,DnaSequence> preffixsuffix() const;
    std::vector<DnaSequence> firstLastSubstr(size_t length, bool) const;
    std::string str() const;

private:
    static std::vector<size_t> _dnaTable;

    static size_t dnaToId(char c)
    {
        return _dnaTable[(size_t)c];
    }

    static char idToDna(size_t id)
    {
        static char table[] = {'A', 'C', 'G', 'T'};
        return table[id];
    }

    struct TableFiller
    {
        TableFiller()
        {
            static bool tableFilled = false;
            if (!tableFilled)
            {
                tableFilled = true;
                _dnaTable.assign(256, -1);	//256 chars
                _dnaTable[(size_t)'A'] = 0;
                _dnaTable[(size_t)'a'] = 0;
                _dnaTable[(size_t)'C'] = 1;
                _dnaTable[(size_t)'c'] = 1;
                _dnaTable[(size_t)'G'] = 2;
                _dnaTable[(size_t)'g'] = 2;
                _dnaTable[(size_t)'T'] = 3;
                _dnaTable[(size_t)'t'] = 3;
                _dnaTable[(size_t)'N'] = 4;
                _dnaTable[(size_t)'n'] = 4;
            }
        }
    };
    static TableFiller _filler;
    static size_t num_pointers;

    SharedBuffer* _data;
    bool _complement, show = false;
};

inline std::string DnaSequence::str() const
{
    std::string result;
    result.reserve(this->length());
    for (size_t i = 0; i < this->length(); ++i)
    {
        result.push_back(this->at(i));
    }
    return result;
}

inline DnaSequence DnaSequence::substr(size_t start, size_t length) const
{
    DnaSequence newSequence;
    if (length == 0)
        return newSequence;
    if (start >= _data->length)
        return newSequence;
    if (start + length > _data->length)
    {
        length = _data->length - start;
    }

    newSequence._data->length = length;
    newSequence._data->chunks.assign((length - 1) / NUCL_IN_CHUNK + 1, 0);

    for (size_t i = 0; i < length; ++i)
    {
        size_t nucId = this->atRaw(start + i);
        size_t newChunkId = i / NUCL_IN_CHUNK;
        newSequence._data->chunks[newChunkId] |= nucId << (i % NUCL_IN_CHUNK) * offset;
    }
    newSequence._complement = _complement;

    return newSequence;
}
/*
 * Preffix-Suffix DnaSequence
 * */
inline std::pair<DnaSequence, DnaSequence> DnaSequence::preffixsuffix() const
{
    DnaSequence preffix((*this)), suffix((*this));
    suffix.append_with_replace_right(0);
    suffix._data->length--;
    preffix.set(0, preffix.length()-1);
    preffix._data->length--;
    return std::pair<DnaSequence, DnaSequence>{preffix, suffix};
}
/*
 * AdHoc method to speed it up
 */
inline std::vector<DnaSequence> DnaSequence::firstLastSubstr(size_t length, bool canonical) const
{
    std::vector<DnaSequence> firstLastSubStrs;
    if (length >= _data->length)
        return firstLastSubStrs;
    DnaSequence firstSubSeq((*this)), lastSubSeq((*this));
    lastSubSeq.append_with_replace_right(0);
    firstSubSeq.set(0,length);
    firstSubSeq._data->length = length;
    lastSubSeq._data->length = length;
    /*
     * Forward
     */
    firstLastSubStrs.push_back(firstSubSeq);
    firstLastSubStrs.push_back(lastSubSeq);
    /*
     * Complementary
     */
    if (canonical)
    {
        DnaSequence comp(*this->complement());
        DnaSequence rc_firstSubSeq(comp), rc_lastSubSeq(comp);
        rc_lastSubSeq.append_with_replace_right(0);
        rc_firstSubSeq.set(0,length);
        rc_firstSubSeq._data->length = length;
        rc_lastSubSeq._data->length = length;
        firstSubSeq._complement = _complement;
        lastSubSeq._complement = _complement;
        /*
         * Complementary
         */
        firstLastSubStrs.push_back(rc_firstSubSeq);
        firstLastSubStrs.push_back(rc_lastSubSeq);
    };
    return firstLastSubStrs;
}