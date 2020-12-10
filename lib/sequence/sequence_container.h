//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#ifndef SEQUENCE_CONTAINER_H
#define SEQUENCE_CONTAINER_H
	#pragma once
	#include <vector>
	#include <unordered_map>
	#include <boost/filesystem.hpp>
	#include <map>
	#include <string>
	#include <limits>
	#include "sequence.h"

	using namespace std;
	namespace fs = boost::filesystem;

	struct FastaRecord
	{
		class Id
		{
		public:
			Id(): _id(std::numeric_limits<uint32_t>::max()) {}
			Id(uint32_t id): _id(id) {}

			bool operator==(const Id& other) const
				{return _id == other._id;}
			bool operator!=(const Id& other) const
				{return !(*this == other);}
			bool operator>(const Id &other) const
				{return _id > other._id;}
			bool operator<(const Id &other) const
				{return _id < other._id;}

			Id rc() const		//reverse complement
				{return Id(_id + 1 - (_id % 2) * 2);}
			Id pr() const   //pair_read
				{return (_id%4<2)?(_id%4==0)?Id(_id+3):Id(_id+1):Id(_id-2);}

			bool strand() const		//true = positive(forward), false = negative
				{return !(_id % 2);}
			size_t hash() const
			{
				size_t x = _id;
				size_t z = (x += 0x9E3779B97F4A7C15ULL);
				z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
				z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
				return z ^ (z >> 31);
			}
			uint32_t getId() const {return _id;}
			int signedId() const
				{return (_id % 2) ? -((int)_id + 1) / 2 : (int)_id / 2 + 1;}

			friend std::ostream& operator << (std::ostream& stream, const Id& id)
			{
				stream << std::to_string(id._id);
				return stream;
			}

			friend std::istream& operator >> (std::istream& stream, Id& id)
			{
				std::string buffer;
				stream >> buffer;
				id._id = std::stoi(buffer);
				return stream;
			}

		private:
			uint32_t _id;
		};
		static const Id ID_NONE;
		typedef std::tuple<Id, Id> IdPair;

		FastaRecord(): id(ID_NONE) {}
		FastaRecord(const DnaSequence& sequence, const std::string& description,
					Id id):
			id(id), sequence(sequence), description(description)
		{
		}
		FastaRecord(const DnaSequence& sequence, const std::string& description,
					Id id, uint8_t sample):
				id(id), sequence(sequence), description(description), sample(sample)
		{
		}

		FastaRecord(const FastaRecord& other):
			id(other.id), sequence(other.sequence),
			description(other.description) {}

		FastaRecord(FastaRecord&& other):
			id (other.id)
		{
			*this = std::move(other);
		}

		FastaRecord& operator=(const FastaRecord& other)
		{
			id = other.id;
			sequence = other.sequence;
			description = other.description;
			return *this;
		}

		FastaRecord& operator=(FastaRecord&& other)
		{
			id = other.id;
			sequence = std::move(other.sequence);
			description = std::move(other.description);
			return *this;
		}

		std::string getDescription(){
			return description;
		}

		Id getComplementaryId() const {
			return id.rc();
		}

		Id getPairId() const {
			return id.pr();
		}

		Id getId() const {
			return id.getId();
		}

		uint8_t getSample()
		{
			return sample;
		}

		Id id;
		DnaSequence sequence;
		std::string description;
		uint8_t sample;
	};

	namespace std
	{
		template <>
		//Esto se hace para que FastaRecord::Id pueda ser usado como claves en mapas y conjuntos
		struct hash<FastaRecord::Id>
		{
			size_t operator() (const FastaRecord::Id& h) const throw()
			{
				 return h.hash();
			}
		};

		template <>
		struct hash<FastaRecord::IdPair>
		{
			 size_t operator()(const FastaRecord::IdPair& k) const
			 {
				size_t lhs = std::get<0>(k).hash();
				size_t rhs = std::get<1>(k).hash();
				lhs ^= rhs + 0x9ddfea08eb382d69ULL + (lhs << 6) + (lhs >> 2);
				return lhs;
			 }
		};
	}

	class SequenceContainer
	{
	public:
		class ParseException : public std::runtime_error
		{
		public:
			ParseException(const std::string & what):
				std::runtime_error(what)
			{}
		};

		//typedef std::unordered_map<FastaRecord::Id,
		//						   FastaRecord> SequenceIndex;
		typedef std::unordered_map<FastaRecord::Id, FastaRecord> SequenceIndex;
		SequenceContainer() {}
		~SequenceContainer(){
			_seqIndex.clear();
		}

		void clear()
		{
			SequenceIndex sc_tmp;
			_seqIndex.clear();
			_seqIndex.swap(sc_tmp);
		}

		void load(const std::string& path, bool);

		void loadFromFile(const std::string& filename, bool,size_t sample = 0, size_t = 1 );
		static void writeFasta(const std::vector<FastaRecord>& records,
							   const std::string& fileName);
		const FastaRecord&  addSequence(const DnaSequence& sequence,
										const std::string& description);

		const FastaRecord& addPairedSequences(std::pair<DnaSequence,DnaSequence> pair_sequences
				,std::pair<std::string,std::string> pair_description);

		void setRead(FastaRecord::Id, DnaSequence);

		const SequenceIndex& getIndex() const
		{
			return _seqIndex;
		}
		const DnaSequence& getSeq(FastaRecord::Id readId) const
		{
			return _seqIndex.at(readId).sequence;
		}
		int32_t seqLen(FastaRecord::Id readId) const
		{
			return _seqIndex.at(readId).sequence.length();
		}
		std::string seqName(FastaRecord::Id readId) const
		{
			return _seqIndex.at(readId).description;
		}
		int computeNxStat(float fraction) const;

		void writeSequenceContainer(const std::string&) const;

		size_t getAvLength()
		{
			return (_totalBases==0)?0:_totalBases/_totalReads;
		}

		size_t getTotalBases()
		{
			return _totalBases;
		}

		size_t size()
		{
			return _seqIndex.size();
		}

		//Operators
		const SequenceContainer& operator=(const SequenceContainer & sc)
		{
			_seqIndex = sc._seqIndex;
			g_nextSeqId = sc.g_nextSeqId;
			return *this;
		}

		//ShowInfo
		void ShowInfo()
		{
			for (auto fr:_seqIndex)
			{
				std::cout << fr.first << " "<<fr.second.sequence.str()<<"\n";
			}
		}

		void write_left_chains();
		void write_right_chains();

	private:
		size_t readFasta(std::vector<FastaRecord>& record
				,const std::string&,bool,size_t, size_t);
		size_t readFastq(std::vector<FastaRecord>& record
				,const std::string&,bool,size_t, size_t);
		bool isFasta(const std::string& fileName);
		void 	validateSequence(std::string& sequence);
		void 	validateHeader(std::string& header);
		size_t _totalBases = 0, _totalReads = 0;
		SequenceIndex _seqIndex;
		size_t g_nextSeqId = 0, g_nextRightSeqId = 2;
	};
#endif
