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
		VariantCaller(const FASTARecord & refSeq, const std::vector<BlockInstance> & blockList, size_t trimK);
		void CallVariants(std::vector<Variant> & variantList) const;
		void CallRearrangements(std::vector<Reversal> & reversal, std::vector<Translocation> & translocation) const;
	private:
		DISALLOW_COPY_AND_ASSIGN(VariantCaller);
		const FASTARecord & refSeq_;
		mutable std::vector<BlockInstance> blockList_;
		size_t trimK_;
		mutable std::vector<saidx_t> indexOut_;
		std::vector<saidx_t> suffixArray_;
		bool ConfirmVariant(StrandIterator refenceStart, StrandIterator referenceEnd, StrandIterator assemblyStart, StrandIterator assemblyEnd) const;
		void AlignSyntenyBlocks(const BlockInstance & reference, const BlockInstance & assembly, std::vector<Variant> & variantList) const;
		void AlignBulgeBranches(size_t blockId, StrandIterator referenceBegin, StrandIterator referenceEnd,
			StrandIterator assemblyBegin, StrandIterator assemblyEnd, std::vector<Variant> & variantList) const;
	};
}

#endif 