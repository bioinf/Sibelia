#ifndef _GRAPH_ALGORITHM_H
#define _GRAPH_ALGORITHM_H

#include "debruijngraph.h"

namespace SyntenyBuilder
{
	class GraphAlgorithm
	{
	public:
		typedef DeBruijnGraph Graph;
		static void Simplify(Graph & g, int minCycleSize);		
		static void ListNonBranchingPaths(Graph & g, std::ostream & out);
	};

	std::ostream& operator << (std::ostream & out, DeBruijnGraph & g);	
}

#endif
