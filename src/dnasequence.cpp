#include "dnasequence.h"

namespace SyntenyBuilder
{
	inline std::string ConstructComplementarityTable()
	{
		std::string ret(1 << 8, ' ');
		ret['a'] = 't';
		ret['t'] = 'a';
		ret['g'] = 'c';
		ret['c'] = 'g';
		ret['A'] = 'T';
		ret['T'] = 'A';
		ret['G'] = 'C';
		ret['C'] = 'G';
		return ret;
	}
	
	const std::string DNASequence::alphabet("agct");	
	const std::string DNASequence::complementary_(ConstructComplementarityTable());
	DNASequence::PositiveReadingStrategy<IndexIterator> DNASequence::positiveReading_;		
	DNASequence::NegativeReadingStrategy<IndexIterator> DNASequence::negativeReading_;		
	DNASequence::PositiveReadingStrategy<IndexConstIterator> DNASequence::positiveConstReading_;
	DNASequence::NegativeReadingStrategy<IndexConstIterator> DNASequence::negativeConstReading_;
}