//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "debruijnindex.h"

#ifndef _BLOCK_BUILDER_H_
#define _BLOCK_BUILDER_H_

namespace SyntenyFinder
{	
	class BlockBuilder
	{
	public:
		BlockBuilder(const std::vector<FastaRecord> * originalChr);
		void ConstructIndex(size_t k);
		void Simplify(size_t minBranchSize, size_t maxIterations);
		void GenerateBlocks(std::vector<BlockInstance> & ret) const;
	private:
		DISALLOW_COPY_AND_ASSIGN(BlockBuilder);
		size_t lastK_;
		DeBruijnIndex * index_;		
		const std::vector<FastaRecord> * originalChr_;
	};
}

#endif