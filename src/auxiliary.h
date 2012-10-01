#ifndef _AUXILIARY_H_
#define _AUXILIARY_H_

#include "fasta.h"
#include "graphalgorithm.h"

namespace SyntenyBuilder
{
	inline bool CompareBlocksById(const SyntenyBuilder::GraphAlgorithm::BlockInstance & a, const SyntenyBuilder::GraphAlgorithm::BlockInstance & b)
	{
		return a.GetBlockId() < b.GetBlockId();
	}

	void GenerateReport(const ChrList & chrList, std::vector<GraphAlgorithm::BlockInstance> & block, std::ostream & out);
}

#endif