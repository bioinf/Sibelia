//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _VARIANT_H_
#define _VARIANT_H_

#include "../blockfinder.h"

namespace SyntenyFinder
{
	class Variant
	{
	public:
		size_t GetBlockId() const;
		size_t GetReferencePos() const;
		const FASTARecord* GetAssemblySequence() const;
		const FASTARecord* GetReferenceSequence() const;
        const std::string& GetReferenceAllele() const;
        const std::string& GetAlternativeAllele() const;		
		const std::string& GetReferenceContext() const;
		const std::string& GetAlternativeContext() const;
		bool Equal(const Variant & toCompare) const;
		Variant(size_t referencePos, size_t blockId, const std::string & referenceAllele, const std::string & alternativeAllele,
			const FASTARecord * referenceSequence, const FASTARecord * assemblySequence, const std::string & referenceContext, const std::string alternativeContext);
		bool operator < (const Variant & toCompare) const;		
		static const size_t UNKNOWN_POS;
		static const size_t UNKNOWN_BLOCK;		
	private:		
		size_t referencePos_;
		size_t blockId_;
		std::string referenceAllele_;
		std::string alternativeAllele_;
		std::string referenceContext_;
		std::string alternativeContext_;
		const FASTARecord * referenceSequence_;
		const FASTARecord * assemblySequence_;
		std::string alignment_;
		friend std::ostream& operator << (std::ostream & out, const Variant & variant);
	};
}

#endif
