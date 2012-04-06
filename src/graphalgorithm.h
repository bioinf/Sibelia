#ifndef _GRAPH_ALGORITHM_H
#define _GRAPH_ALGORITHM_H

#include "debruijngraph.h"

namespace SyntenyBuilder
{
	class GraphAlgorithm
	{
	public:
		static void Simplify(DeBruijnGraph & g, int minCycleSize);
		static void ListNonBranchingPaths(const DeBruijnGraph & g);
		static void DebugOutput(const DeBruijnGraph & g, std::ostream & out);
	};
}

#endif _GRAPH_ALGORITHM_H