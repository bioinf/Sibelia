#include "syntenyutility.h"

namespace SyntenyBuilder
{
	namespace 
	{
		std::string ConstructComplementarityTable()
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
	}

	const char DELETED_CHARACTER = -1;
	const std::string alphabet_("acgt");	
	const std::string complement_(ConstructComplementarityTable());
}