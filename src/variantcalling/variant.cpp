//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "variant.h"

namespace SyntenyFinder
{
	const size_t Variant::UNKNOWN_BLOCK = -1;

	Variant::Variant(size_t refPos, size_t blockId, bool collinear, const std::string & refAllele, const std::string & altAllele, const FASTARecord & sequence, const std::string & alignment):
		refPos_(refPos), blockId_(blockId), collinear_(collinear), refAllele_(refAllele), altAllele_(altAllele), sequence_(&sequence), alignment_(alignment)
	{

	}

	bool Variant::operator < (const Variant & toCompare) const
	{
		return refPos_ < toCompare.refPos_;
	}

	size_t Variant::GetBlockId() const
	{
		return blockId_;
	}

	size_t Variant::GetReferencePos() const
	{
		return refPos_;
	}

	const std::string& Variant::GetAlignment() const
	{
		return alignment_;
	}

    const std::string & Variant::GetReferenceAllele() const
    {
        return refAllele_;
    }

    const std::string & Variant::GetAlternativeAllele() const
    {
        return altAllele_;
    }

	bool Variant::Equal(const Variant & toCompare) const
	{
		return toCompare.refPos_ == refPos_ && refAllele_ == toCompare.refAllele_ && altAllele_ == toCompare.altAllele_;
	}

	const FASTARecord & Variant::GetSequence() const
	{
		return *sequence_;
	}

	std::ostream& operator << (std::ostream & out, const Variant & variant)
	{
		out << variant.refPos_ << '\t' << (variant.refAllele_.empty() ? "." : variant.refAllele_) << '\t';
		return out << (variant.altAllele_.empty() ? "." : variant.altAllele_) << '\t' << variant.blockId_ << '\t' << variant.GetSequence().GetDescription();
	}
}
