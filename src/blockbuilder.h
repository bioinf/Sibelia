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
		enum State
		{
			start,
			run,
			end
		};

		typedef boost::function<void(size_t, State)> ProgressCallBack;

		BlockBuilder(const std::vector<FastaRecord> * originalChr, const std::string & tempDir);
		void ConstructIndex(size_t k);		
		void GenerateBlocks(std::vector<BlockInstance> & ret, size_t minBlockSize) const;
		void Simplify(size_t minBranchSize, size_t maxIterations, ProgressCallBack f = ProgressCallBack());
	private:
		DISALLOW_COPY_AND_ASSIGN(BlockBuilder);
		size_t lastK_;
		std::auto_ptr<DeBruijnIndex> index_;		
		const std::vector<FastaRecord> * originalChr_;
		std::string tempDir_;
	};
}

#endif