#ifndef _GRAPH_ALGORITHM_H
#define _GRAPH_ALGORITHM_H

#include "kmerindex.h"

namespace SyntenyBuilder
{
	class GraphAlgorithm
	{
	public:		
		static void SimplifyGraph(KMerIndex & index, DNASequence & sequence, size_t minBranchSize);		
		static void SerializeGraph(const KMerIndex & index, const DNASequence & seq, std::ostream & out);
		static void ListNonBranchingPaths(const KMerIndex & index, const DNASequence & sequence, std::ostream & out, std::ostream & indexOut);
	};
}

#endif
