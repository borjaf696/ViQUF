//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iostream>
#include <random>
#include <algorithm>
#include <zlib.h>
#include <queue>

#include "sequence_container.h"

//TODO: Anhadir la parte del procesamiento de las Ns


const FastaRecord::Id FastaRecord::ID_NONE = 
			Id(std::numeric_limits<uint32_t>::max());


bool SequenceContainer::isFasta(const std::string& fileName)
{
	std::string withoutGz = fileName;
	if (fileName.substr(fileName.size() - 3) == ".gz")
	{
		withoutGz = fileName.substr(0, fileName.size() - 3);
	}

	size_t dotPos = withoutGz.rfind(".");
	if (dotPos == std::string::npos)
	{
		throw ParseException("Can't identify input file type");
	}
	std::string suffix = withoutGz.substr(dotPos + 1);
	if (suffix == "fasta" || suffix == "fa")
	{
		return true;
	}
	else if (suffix == "fastq" || suffix == "fq")
	{
		return false;
	}
	throw ParseException("Can't identify input file type");
}

void SequenceContainer::load(const std::string &path, bool is_paired)
{
    fs::path doc_path(path);
    if (fs::is_directory(doc_path))
    {
        std::queue<fs::path> dir_queue;
        dir_queue.emplace(doc_path);
        size_t num_files_checked = 0;
        std::cout << "Is directory: ";
        (is_paired)?std::cout << " Pair_end\n":std::cout << " Single_End\n";
        while(!dir_queue.empty())
        {
			doc_path = dir_queue.front();
			dir_queue.pop();
            for (auto &f : fs::directory_iterator(doc_path)) {
                if (fs::is_directory(f)) {
					dir_queue.emplace(f);
                    continue;
                }
                num_files_checked++;
                std::ostringstream oss;
                oss << f;
                std::string converted_path = oss.str().substr(1, oss.str().size() - 2);
                std::cout << "File: " << converted_path << "\n";
                if (!is_paired)
                    loadFromFile(converted_path, is_paired);
                else {
                    std::string side_read = converted_path.substr(converted_path.rfind('.') - 1, 1);
                    std::cout << "Read: " << side_read << "\n";
                    if (side_read != "1" && side_read != "2")
                        throw ParseException("Fail in pair_end reads");
                    else if (side_read == "1")
                        loadFromFile(converted_path, is_paired);
                    else
                        loadFromFile(converted_path, is_paired, 0);
                }
            }
        }
        if (num_files_checked % 2 && is_paired)
        {
            std::cout << "Wrong number of Pair_End files\n";
            exit(1);
        }
    }
	if (fs::is_regular_file(doc_path))
    {
		if (is_paired){
			std::cout << "-.-' It is only a file dude!\n";
			exit(1);
		}
        std::cout << path<<"->Sequence Container\n";
        loadFromFile(path, false);
    }
	if (!fs::exists(doc_path))
		std::cout << "Path does not exist\n";
}

void SequenceContainer::loadFromFile(const std::string& fileName, bool is_paired,size_t sample, size_t side_read)
{
	std::vector<FastaRecord> records;
	if (this->isFasta(fileName))
		this->readFasta(records, fileName, is_paired, side_read, sample);
	else
		this->readFastq(records, fileName, is_paired, side_read, sample);
    if (side_read)
	    records.reserve(records.size() * (is_paired)?4:2);
	std::vector<FastaRecord> complements;
	for (auto &record : records)
	{
		std::string header = "-" + record.description;
		record.description = "+" + record.description;
		DnaSequence revComplement = record.sequence.complement();
		_totalBases += 2*record.sequence.length();
		_totalReads += 2;
		complements.push_back(FastaRecord((is_paired)?revComplement:DnaSequence(""), header, record.id.rc()));
	}
	for (auto& rec : complements)
	{
		records.push_back(std::move(rec));
	}
	
	//shuffling input reads
	std::vector<size_t> indicesPerm(records.size());
	for (size_t i = 0; i < indicesPerm.size(); ++i) indicesPerm[i] = i;
	std::random_shuffle(indicesPerm.begin(), indicesPerm.end());
	//

	//_seqIndex.reserve(_seqIndex.size() + records.size());
	for (size_t i : indicesPerm)
	{
		//_seqIndex[records[i].id] = std::move(records[i]);
		_seqIndex.emplace(std::pair<FastaRecord::Id,FastaRecord>
							(records[i].id, records[i]));
	}
}


int SequenceContainer::computeNxStat(float fraction) const
{
	std::vector<int32_t> readLengths;
	int64_t totalLengh = 0;
	for (auto& read : _seqIndex) 
	{
		readLengths.push_back(read.second.sequence.length());
		totalLengh += read.second.sequence.length();
	}
	std::sort(readLengths.begin(), readLengths.end(),
			  [](int32_t a, int32_t b) {return a > b;});

	int32_t nx = 0;
    int64_t cummulativeLen = 0;
	for (auto l : readLengths)
	{
        cummulativeLen += l;
        if (cummulativeLen > fraction * totalLengh)
		{
            nx = l;
            break;
		}
	}
	return nx;
}


const FastaRecord& 
	SequenceContainer::addSequence(const DnaSequence& sequence, 
								   const std::string& description)
{
	FastaRecord::Id newId = FastaRecord::Id(g_nextSeqId);
	FastaRecord fwdRecord(sequence, "+" + description, newId);
	_seqIndex[fwdRecord.id] = std::move(fwdRecord);

	FastaRecord revRecord(sequence.complement(), "-" + description, 
						  newId.rc());
	_seqIndex[revRecord.id] = std::move(revRecord);

	g_nextSeqId += 2;
	return _seqIndex[newId];
}

const FastaRecord&
SequenceContainer::addPairedSequences(std::pair<DnaSequence, DnaSequence> pair_sequences,
                                      std::pair<std::string, std::string> pair_description)
{
    /*
     * First_read + reverse_complement
     */
    FastaRecord::Id newId = FastaRecord::Id(g_nextSeqId);
    FastaRecord fwdRecord(pair_sequences.first, "+" + pair_description.first, newId);
    _seqIndex[fwdRecord.id] = std::move(fwdRecord);

    FastaRecord revRecord(pair_sequences.first.complement(), "-" + pair_description.first,
                          newId.rc());
    _seqIndex[revRecord.id] = std::move(revRecord);
    /*
     * Second_read + reverse_complement
     */
    g_nextSeqId += 2;
    FastaRecord::Id paired_id = FastaRecord::Id(g_nextSeqId);
    FastaRecord paired_fwdRecord(pair_sequences.second, "+" + pair_description.second, paired_id);
    _seqIndex[paired_fwdRecord.id] = std::move(paired_fwdRecord);

    FastaRecord paired_revRecord(pair_sequences.second.complement(), "-" + pair_description.second,
                          paired_id.rc());
    _seqIndex[paired_revRecord.id] = std::move(paired_revRecord);
    /*
     * Update seqId
     */
    g_nextSeqId += 2;
    return _seqIndex[newId];
}

size_t SequenceContainer::readFasta(std::vector<FastaRecord>& record
        ,const std::string& fileName, bool is_paired, size_t side_read, size_t sample)
{
	size_t BUF_SIZE = 32 * 1024 * 1024;
	char* rawBuffer = new char[BUF_SIZE];
	auto* fd = gzopen(fileName.c_str(), "rb");
	if (!fd)
	{
		throw ParseException("Can't open reads file");
	}

	record.clear();
	int lineNo = 1;
	std::string header; 
	std::string sequence;
	try
	{
		while(!gzeof(fd))
		{
			//get a new line
			std::string nextLine;
			for (;;)
			{
				char* read = gzgets(fd, rawBuffer, BUF_SIZE);
				if (!read) break;
				nextLine += read;
				if (nextLine.empty()) break;
				if (nextLine.back() == '\n')
				{
					nextLine.pop_back();
					break;
				}
			}

			if (nextLine.empty()) continue;
			if (nextLine.back() == '\r') nextLine.pop_back();

			if (nextLine[0] == '>')
			{
				if (!header.empty())
				{
					if (sequence.empty()) throw ParseException("empty sequence");

					record.push_back(FastaRecord(DnaSequence(sequence), header,
										FastaRecord::Id((side_read)?g_nextSeqId:g_nextRightSeqId),sample));
                    (side_read)?g_nextSeqId += (is_paired)?4:2:g_nextRightSeqId+=4;
					sequence.clear();
					header.clear();
				}
				this->validateHeader(nextLine);
				header = nextLine;
			}
			else
			{
				this->validateSequence(nextLine);
				std::copy(nextLine.begin(), nextLine.end(), 
						  std::back_inserter(sequence));
			}
			++lineNo;
		}
		
		if (sequence.empty()) throw ParseException("empty sequence");
		if (header.empty())
		{
			throw ParseException("Fasta fromat error");
		}
		record.push_back(FastaRecord(DnaSequence(sequence), header,
							FastaRecord::Id((side_read)?g_nextSeqId:g_nextRightSeqId)));
        (side_read)?g_nextSeqId += (is_paired)?4:2:g_nextRightSeqId+=4;
	}
	catch (ParseException& e)
	{
		std::stringstream ss;
		ss << "parse error in " << fileName << " on line " << lineNo << ": " << e.what();
		gzclose(fd);
		throw ParseException(ss.str());
	}

	delete[] rawBuffer;
	gzclose(fd);
	return record.size();
}

size_t SequenceContainer::readFastq(std::vector<FastaRecord>& record
        ,const std::string& fileName, bool is_paired, size_t side_read, size_t sample)
{
	size_t BUF_SIZE = 32 * 1024 * 1024;
	char* rawBuffer = new char[BUF_SIZE];
	auto* fd = gzopen(fileName.c_str(), "rb");
	if (!fd)
	{
		throw ParseException("Can't open reads file");
	}
	record.clear();
	int lineNo = 1;
	int stateCounter = 0;
	std::string header; 
	try
	{
		while (!gzeof(fd))
		{
			//get a new line
			std::string nextLine;
			for (;;)
			{
				char* read = gzgets(fd, rawBuffer, BUF_SIZE);
				if (!read) {
					break;
				}
				nextLine +=read;
				if (nextLine.empty()) break;
				if (nextLine.back() == '\n')
				{
					nextLine.pop_back();
					break;
				}
			}

			if (nextLine.empty()) 
			{
				stateCounter = (stateCounter + 1) % 4;
				continue;
			}
			if (nextLine.back() == '\r') nextLine.pop_back();

			if (stateCounter == 0)
			{
				if (nextLine[0] != '@') throw ParseException("Fastq format error");
				header = nextLine;
				this->validateHeader(header);
			}
			else if (stateCounter == 1)
			{
				this->validateSequence(nextLine);
				record.push_back(FastaRecord(DnaSequence(nextLine), header,
									FastaRecord::Id((side_read)?g_nextSeqId:g_nextRightSeqId), sample));
                (side_read)?g_nextSeqId += (is_paired)?4:2:g_nextRightSeqId+=4;
			}
			else if (stateCounter == 2)
			{
				if (nextLine[0] != '+') throw ParseException("Fastq fromat error");
			}
			stateCounter = (stateCounter + 1) % 4;
			++lineNo;
		}
	}
	catch (ParseException& e)
	{
		std::stringstream ss;
		ss << "parse error in " << fileName << " on line " << lineNo << ": " << e.what();
		gzclose(fd);
		throw ParseException(ss.str());
	}
	gzclose(fd);
	delete[] rawBuffer;
	return record.size();
}


void SequenceContainer::validateHeader(std::string& header)
{
	size_t delim = header.find(' ');
	if (delim == std::string::npos)
	{
		delim = header.length() - 1;
	}
	else
	{
		--delim;
	}

	header = header.substr(1, delim);
	if (header.empty()) throw ParseException("empty header");
}

void SequenceContainer::validateSequence(std::string& sequence)
{
	static const std::string VALID_CHARS = "ACGTN";
	for (size_t i = 0; i < sequence.length(); ++i)
	{
		if (VALID_CHARS.find(toupper(sequence[i])) == std::string::npos) 
		{
			sequence[i] = VALID_CHARS[rand() % 4];
		}
	}
}

void SequenceContainer::writeFasta(const std::vector<FastaRecord>& records, 
								   const std::string& filename)
{
	static const size_t FASTA_SLICE = 80;

	FILE* fout = fopen(filename.c_str(), "w");
	if (!fout) throw std::runtime_error("Can't open " + filename);
	
	for (auto& rec : records)
	{
		std::string contigSeq;
		for (size_t c = 0; c < rec.sequence.length(); c += FASTA_SLICE)
		{
			contigSeq += rec.sequence.substr(c, FASTA_SLICE).str() + "\n";
		}
		std::string header = ">" + rec.description + "\n";
		fwrite(header.data(), sizeof(header.data()[0]), 
			   header.size(), fout);
		fwrite(contigSeq.data(), sizeof(contigSeq.data()[0]), 
			   contigSeq.size(), fout);
	}
	fclose(fout);
}

void SequenceContainer::setRead(FastaRecord::Id readId, DnaSequence newSeq)
{
    _totalBases+=(_seqIndex[readId]).sequence.length()-newSeq.length();
	_seqIndex[readId] = FastaRecord(newSeq,(_seqIndex[readId]).getDescription()
			,FastaRecord::Id(readId));
}

void SequenceContainer::writeSequenceContainer(const std::string & filename) const
{
	FILE* fout = fopen(filename.c_str(),"w");
	if (!fout) throw std::runtime_error("Can't open " + filename);
	size_t cont = 0;
	for (auto& id_rec: _seqIndex)
	{
		if (id_rec.first.getId() % 2 == 0)
		{
			std::string corrected_read = id_rec.second.sequence.str()+"\n";
			std::string header = ">" + id_rec.second.description+"\n";
			fwrite(header.data(), sizeof(header.data()[0]), header.size(),fout);
			fwrite(corrected_read.data(), sizeof(corrected_read.data()[0]),corrected_read.size(),fout);
		}
		cont++;
	}
    fclose(fout);
}

void SequenceContainer::write_left_chains()
{
	FILE * fout = fopen("left_reads.fa","w");
	if (!fout) throw std::runtime_error("Cant't open left_reads.fa");
	for (auto & id_rec:_seqIndex)
	{
		if (id_rec.first.getId() % 4 == 0)
		{
			std::string corrected_read = id_rec.second.sequence.str()+"\n";
			std::string header = ">"+std::to_string(id_rec.first.getId())+"\n";
			fwrite(header.data(), sizeof(header.data()[0]), header.size(), fout);
			fwrite(corrected_read.data(), sizeof(corrected_read.data()[0]), corrected_read.size(), fout);
		}
	}
	fclose(fout);
}

void SequenceContainer::write_right_chains()
{
	FILE * fout = fopen("right_reads.fa","w");
	if (!fout) throw std::runtime_error("Cant't open right_reads.fa");
	for (auto & id_rec:_seqIndex)
	{
		if (id_rec.first.getId() % 4 == 0)
		{
			std::string corrected_read = getSeq(id_rec.second.getPairId()).str()+"\n";
			std::string header = ">"+std::to_string(id_rec.first.getId())+"\n";
			fwrite(header.data(), sizeof(header.data()[0]), header.size(), fout);
			fwrite(corrected_read.data(), sizeof(corrected_read.data()[0]), corrected_read.size(), fout);
		}
	}
	fclose(fout);
}