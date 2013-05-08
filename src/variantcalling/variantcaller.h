//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _VARIANT_CALLER_
#define _VARIANT_CALLER_

#include "variant.h"
#include "rearrangement.h"

namespace SyntenyFinder
{
	class VariantCaller
	{
	public:		
		VariantCaller(const std::vector<FASTARecord> & chr, const std::set<size_t> & referenceSequenceId, const std::vector<std::vector<BlockInstance> > & history, size_t trimK, size_t minBlockSize);
		void CallVariants(std::vector<Variant> & variantList);
		void GetHistory(std::vector<std::vector<BlockInstance> > & history);
		void CallRearrangements(std::vector<Reversal> & reversal, std::vector<Translocation> & translocation);
	private:
		DISALLOW_COPY_AND_ASSIGN(VariantCaller);
		const std::vector<FASTARecord> * chr_;
		std::set<size_t> referenceSequenceId_;
		std::vector<std::vector<BlockInstance> > history_;
		size_t trimK_;
		size_t minBlockSize_;
		std::vector<saidx_t> indexOut_;
		std::vector<saidx_t> referenceSuffixArray_;
		std::vector<saidx_t> assemblySuffixArray_;
		std::string assemblySequence_;
		std::string referenceSequence_;
		bool SearchInAssembly(const std::string & pattern);
		bool SearchInReference(const std::string & pattern);
		const BlockInstance* PreviousBlock(const BlockInstance & block, const std::vector<BlockInstance> & blockList);
		const BlockInstance* NextBlock(const BlockInstance & block, const std::vector<BlockInstance> & blockList);
		std::pair<size_t, size_t> DetermineLeftProbableBoundaries(std::vector<BlockInstance> & blockList, size_t block);
		std::pair<size_t, size_t> DetermineRightProbableBoundaries(std::vector<BlockInstance> & blockList, size_t block);
		void CorrectBlocksBoundaries(std::vector<BlockInstance> & blockList, size_t referenceBlock, size_t assemblyBlock);
		bool ConfirmVariant(StrandIterator refenceStart, StrandIterator referenceEnd, StrandIterator assemblyStart, StrandIterator assemblyEnd);
		void UpdateBlockBoundaries(BlockInstance & block, std::pair<size_t, size_t> leftBoundaries, std::pair<size_t, size_t> rightBoundaries, std::pair<size_t, size_t> startAlignmentCoords, std::pair<size_t, size_t> endAlignmentCoords);
		void LocalAlignment(const std::string & sequence1, const std::string & sequence2, std::pair<size_t, size_t> & coord1, std::pair<size_t, size_t> & coord2);
		void GetBoundariesSequence(const BlockInstance & block, std::pair<size_t, size_t> leftBoundaries, std::pair<size_t, size_t> rightBoundaries, std::string & start, std::string & end);
		void AlignSyntenyBlocks(const BlockInstance & reference, const BlockInstance & assembly, std::vector<Variant> & variantList);
		void AlignBulgeBranches(size_t blockId, StrandIterator referenceBegin, StrandIterator referenceEnd, StrandIterator assemblyBegin, StrandIterator assemblyEnd, const FASTARecord & referenceSequence, const FASTARecord & assemblySequence, std::vector<Variant> & variantList);
	};
}

#endif 