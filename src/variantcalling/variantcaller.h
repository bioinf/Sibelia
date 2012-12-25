//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _VARIANT_CALLER_
#define _VARIANT_CALLER_

#include "variant.h"

namespace SyntenyFinder
{
	class VariantCaller
	{
	public:
		void CallVariants(std::vector<Variant> & variantList) const;
		VariantCaller(size_t refSeqId, const std::vector<BlockInstance> & blockList, size_t trimK);
	private:
		size_t refSeqId_;
		mutable std::vector<BlockInstance> blockList_;
		size_t trimK_;
		
		void AlignSyntenyBlocks(const BlockInstance & reference, const BlockInstance & assembly, std::vector<Variant> & variantList) const;
		void AlignBulgeBranches(size_t blockId, StrandIterator referenceBegin, StrandIterator referenceEnd, StrandIterator assemblyBegin, StrandIterator assemblyEnd, std::vector<Variant> & variantList) const;
	};
}

#endif 