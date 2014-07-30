//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "blockinstance.h"
#include "indexedsequence.h"

//#include <sparsehash/sparse_hash_map>

#ifndef _DE_BRUIJN_GRAPH_H_
#define _DE_BRUIJN_GRAPH_H_

namespace SyntenyFinder
{	
	class DeBruijnIndex
	{
	public:
		static const uint32_t NO_BIFURCATION;
		uint32_t GetBifurcationId(FastaRecord::Iterator it) const;
		void GetBifurcationPositions(uint32_t bifId, std::vector<FastaRecord::Iterator> & pos) const;
		DeBruijnIndex(DeBruijnIndex & g, size_t k);
		DeBruijnIndex(const std::vector<FastaRecord> & chrList, size_t k, const std::string & tempDir);
	private:

	};

	void SimplifyGraph(DeBruijnIndex & g, size_t minBranchSize);	
	std::vector<BlockInstance> GenerateSyntenyBlocks(const DeBruijnIndex & g);
}

#endif
