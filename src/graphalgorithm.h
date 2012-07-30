#ifndef _GRAPH_ALGORITHM_H
#define _GRAPH_ALGORITHM_H

#include "kmerindex.h"

namespace SyntenyBuilder
{
	class GraphAlgorithm
	{
	public:		
		static void SimplifyGraph(DNASequence & sequence, size_t k, size_t minBranchSize);		
		static void SerializeGraph(const DNASequence & seq, size_t k, std::ostream & out);
		static void ListNonBranchingPaths(const DNASequence & sequence, size_t k, std::ostream & out, std::ostream & indexOut);
	};
}

#endif
