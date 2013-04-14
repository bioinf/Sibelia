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
		VariantCaller(const std::vector<FASTARecord> & chr, size_t refSeqId, const std::vector<BlockInstance> & blockList,
			const std::vector<BlockInstance> & smallBlockList, size_t trimK, size_t minBlockSize);
		void CallVariants(std::vector<Variant> & variantList) const;
		void GetBlockList(std::vector<BlockInstance> & blockList) const;
		void CallRearrangements(std::vector<Reversal> & reversal, std::vector<Translocation> & translocation) const;
	private:
		DISALLOW_COPY_AND_ASSIGN(VariantCaller);
		const std::vector<FASTARecord> * chr_;
		size_t refSeqId_;
		mutable std::vector<BlockInstance> blockList_;
		size_t trimK_;
		size_t minBlockSize_;
		mutable std::vector<saidx_t> indexOut_;
		std::vector<saidx_t> referenceSuffixArray_;
		std::vector<saidx_t> assemblySuffixArray_;
		std::string assemblySequence_;
		bool SearchInAssembly(const std::string & pattern) const;
		bool SearchInReference(const std::string & pattern) const;		
		bool ConfirmVariant(StrandIterator refenceStart, StrandIterator referenceEnd, StrandIterator assemblyStart, StrandIterator assemblyEnd) const;
		void AlignSyntenyBlocks(const BlockInstance & reference, const BlockInstance & assembly, std::vector<Variant> & variantList) const;
		void AlignBulgeBranches(size_t blockId, StrandIterator referenceBegin, StrandIterator referenceEnd,
			StrandIterator assemblyBegin, StrandIterator assemblyEnd, const FASTARecord & assemblySequence, std::vector<Variant> & variantList) const;
	};
}

#endif 