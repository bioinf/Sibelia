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
		const FASTARecord & GetSequence() const;
		const std::string & GetAlignment() const;		
        const std::string & GetReferenceAllele() const;
        const std::string & GetAlternativeAllele() const;		
		bool Equal(const Variant & toCompare) const;
		Variant(size_t refPos, size_t blockId, bool collinear, const std::string & refAllele, const std::string & altAllele, const FASTARecord & sequence, const std::string & alignment = "");
		bool operator < (const Variant & toCompare) const;
	private:
		static const size_t UNKNOWN_BLOCK;
		size_t refPos_;
		size_t blockId_;
		bool collinear_;
		std::string refAllele_;
		std::string altAllele_;
		const FASTARecord * sequence_;
		std::string alignment_;				
		friend std::ostream& operator << (std::ostream & out, const Variant & variant);
	};
}

#endif
