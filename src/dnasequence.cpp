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
	
	const size_t DNASequence::HASH_BASE = 57;
	const char DNASequence::DELETED_CHARACTER = -1;
	const std::string DNASequence::alphabet("agct");	
	const std::string DNASequence::complementary_(ConstructComplementarityTable());
}