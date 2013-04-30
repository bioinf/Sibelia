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
		void CallVariants(std::vector<Variant> & variantList) const;
		void GetHistory(std::vector<std::vector<BlockInstance> > & history) const;
		void CallRearrangements(std::vector<Reversal> & reversal, std::vector<Translocation> & translocation) const;
	private:
		DISALLOW_COPY_AND_ASSIGN(VariantCaller);
		const std::vector<FASTARecord> * chr_;
		std::set<size_t> referenceSequenceId_;
		mutable std::vector<std::vector<BlockInstance> > history_;
		size_t trimK_;
		size_t minBlockSize_;
		mutable std::vector<saidx_t> indexOut_;
		std::vector<saidx_t> referenceSuffixArray_;
		std::vector<saidx_t> assemblySuffixArray_;
		std::string assemblySequence_;
		std::string referenceSequence_;
		bool SearchInAssembly(const std::string & pattern) const;
		bool SearchInReference(const std::string & pattern) const;		
		bool ConfirmVariant(StrandIterator refenceStart, StrandIterator referenceEnd, StrandIterator assemblyStart, StrandIterator assemblyEnd) const;
		void AlignSyntenyBlocks(const BlockInstance & reference, const BlockInstance & assembly, std::vector<Variant> & variantList) const;
		void AlignBulgeBranches(size_t blockId, StrandIterator referenceBegin, StrandIterator referenceEnd,
			StrandIterator assemblyBegin, StrandIterator assemblyEnd, const FASTARecord & assemblySequence, std::vector<Variant> & variantList) const;
	};
}

#endif 