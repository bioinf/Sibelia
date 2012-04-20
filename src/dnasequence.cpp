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
	
	const char DNASequence::INVALID_HASH = -1;
	const char DNASequence::DELETED_CHARACTER = -1;
	const size_t DNASequence::HASH_BASE = 57;
	const size_t DNASequence::MOD = (size_t(1) << (sizeof(size_t) * 8 - 2));
	const std::string DNASequence::alphabet("agct");	
	const std::string DNASequence::complementary_(ConstructComplementarityTable());
}