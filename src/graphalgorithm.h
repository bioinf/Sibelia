#ifndef _GRAPH_ALGORITHM_H
#define _GRAPH_ALGORITHM_H

#include "kmerindex.h"

namespace SyntenyBuilder
{
	class GraphAlgorithm
	{
	public:		
		static int SimplifyGraph(KMerIndex & index, DNASequence & sequence, size_t minCycleSize);		
		static void SerializeGraph(KMerIndex & index, const DNASequence & seq, std::ostream & out);
		static void ListNonBranchingPaths(const KMerIndex & index, const DNASequence & sequence, std::ostream & out, std::ostream & indexOut);
	};

}

#endif
