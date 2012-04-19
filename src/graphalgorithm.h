#ifndef _GRAPH_ALGORITHM_H
#define _GRAPH_ALGORITHM_H

#include "debruijngraph.h"

namespace SyntenyBuilder
{
	class GraphAlgorithm
	{
	public:
		typedef DeBruijnGraph Graph;
		static int Simplify(Graph & g, int minCycleSize);		
		static void ListNonBranchingPaths(Graph & g, std::ostream & out, std::ostream & indexOut);
	};

	std::ostream& operator << (std::ostream & out, DeBruijnGraph & g);	
}

#endif
