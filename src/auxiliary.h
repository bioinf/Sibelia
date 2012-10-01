#ifndef _AUXILIARY_H_
#define _AUXILIARY_H_

#include "fasta.h"
#include "graphalgorithm.h"

namespace SyntenyFinder
{
	inline bool CompareBlocksById(const SyntenyFinder::GraphAlgorithm::BlockInstance & a, const SyntenyFinder::GraphAlgorithm::BlockInstance & b)
	{
		return a.GetBlockId() < b.GetBlockId();
	}

	void GenerateReport(const ChrList & chrList, std::vector<GraphAlgorithm::BlockInstance> & block, std::ostream & out);
}

#endif