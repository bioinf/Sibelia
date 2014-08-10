//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "blockbuilder.h"

namespace SyntenyFinder
{
	BlockBuilder::BlockBuilder(const std::vector<FastaRecord> * originalChr):
		index_(0), originalChr_(originalChr)
	{
	}

	void BlockBuilder::ConstructIndex(size_t k)
	{
	}

	void BlockBuilder::Simplify(size_t minBranchSize, size_t maxIterations)
	{
	}

	void BlockBuilder::GenerateBlocks(std::vector<BlockInstance> & ret) const
	{
	}
}