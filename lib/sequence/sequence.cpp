#include "sequence.h"

std::vector<size_t> DnaSequence::_dnaTable;
DnaSequence::TableFiller DnaSequence::_filler;
size_t DnaSequence::num_pointers = 0;
