//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _POSTPROCESSOR_H_
#define _POSTPROCESSOR_H_

#include <seqan/align.h>
#undef min
#undef max
#include "fasta.h"
#include "blockinstance.h"

namespace SyntenyFinder
{
	class Postprocessor
	{
	public:		
		Postprocessor(const std::vector<FASTARecord> & chr, size_t minBlockSize);
		void GlueStripes(std::vector<BlockInstance> & block);
		bool ImproveBlockBoundaries(std::vector<BlockInstance> & block);
		void MatchRepeats(std::vector<BlockInstance> & block, const std::set<size_t> & referenceSequenceId);		
	private:
		DISALLOW_COPY_AND_ASSIGN(Postprocessor);
		const std::vector<FASTARecord> * chr_;
		std::set<size_t> referenceSequenceId_;
		std::vector<std::vector<BlockInstance> > history_;
		size_t minBlockSize_;
		size_t correctionRange_;
		void ImprovePairwiseBlockBoundaries(std::vector<BlockInstance> & block);
		const BlockInstance* PreviousBlock(const BlockInstance & block, const std::vector<BlockInstance> & blockList);
		const BlockInstance* NextBlock(const BlockInstance & block, const std::vector<BlockInstance> & blockList);
		std::pair<size_t, size_t> DetermineLeftProbableBoundaries(const std::vector<BlockInstance> & blockList, size_t block);
		std::pair<size_t, size_t> DetermineRightProbableBoundaries(const std::vector<BlockInstance> & blockList, size_t block);
		std::pair<BlockInstance, BlockInstance> CorrectBlocksBoundaries(const std::vector<BlockInstance> & blockList, size_t referenceBlock, size_t assemblyBlock);
		BlockInstance UpdateBlockBoundaries(const BlockInstance & block, std::pair<size_t, size_t> leftBoundaries, std::pair<size_t, size_t> rightBoundaries, std::pair<size_t, size_t> startAlignmentCoords, std::pair<size_t, size_t> endAlignmentCoords);
		void LocalAlignment(const std::string & sequence1, const std::string & sequence2, std::pair<size_t, size_t> & coord1, std::pair<size_t, size_t> & coord2);
		void GetBoundariesSequence(const BlockInstance & block, std::pair<size_t, size_t> leftBoundaries, std::pair<size_t, size_t> rightBoundaries, std::string & start, std::string & end);
	};

}

#endif