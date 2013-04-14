//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "variant.h"

namespace SyntenyFinder
{
	const size_t Variant::UNKNOWN_BLOCK = -1;

	Variant::Variant(size_t referencePos, size_t blockId, const std::string & referenceAllele, const std::string & alternativeAllele,
		const FASTARecord & sequence, const std::string & referenceContext, const std::string alternativeContext):
		referencePos_(referencePos), blockId_(blockId), referenceAllele_(referenceAllele), alternativeAllele_(alternativeAllele),
			sequence_(&sequence), referenceContext_(referenceContext), alternativeContext_(alternativeContext)
	{

	}

	bool Variant::operator < (const Variant & toCompare) const
	{
		return referencePos_ < toCompare.referencePos_;
	}

	size_t Variant::GetBlockId() const
	{
		return blockId_;
	}

	size_t Variant::GetReferencePos() const
	{
		return referencePos_;
	}

    const std::string & Variant::GetReferenceAllele() const
    {
        return referenceAllele_;
    }

    const std::string & Variant::GetAlternativeAllele() const
    {
        return alternativeAllele_;
    }

	bool Variant::Equal(const Variant & toCompare) const
	{
		return toCompare.referencePos_ == referencePos_ && referenceAllele_ == toCompare.referenceAllele_ && alternativeAllele_ == toCompare.alternativeAllele_;
	}

	const FASTARecord & Variant::GetSequence() const
	{
		return *sequence_;
	}

	const std::string & Variant::GetReferenceContext() const
	{
		return referenceContext_;
	}
	
	const std::string & Variant::GetAlternativeContext() const
	{
		return alternativeContext_;
	}

	std::ostream& operator << (std::ostream & out, const Variant & variant)
	{
		out << variant.referencePos_ << '\t' << (variant.referenceAllele_.empty() ? "." : variant.referenceAllele_) << '\t';
		out << (variant.alternativeAllele_.empty() ? "." : variant.alternativeAllele_) << '\t';
		if(variant.blockId_ == Variant::UNKNOWN_BLOCK)
		{
			out << "None";
		}
		else
		{
			out << variant.blockId_;
		}
		
		out << '\t' << variant.GetSequence().GetDescription();
		out << '\t' << (variant.referenceContext_.empty() ? "." : variant.referenceContext_);
		return out << '\t' << (variant.alternativeContext_.empty() ? "." : variant.alternativeContext_);
	}
}
