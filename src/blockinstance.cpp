//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "fasta.h"
#include "blockinstance.h"

namespace SyntenyFinder
{
	typedef boost::function<size_t(const BlockInstance&)> SizeF;
	SizeF getId = boost::bind(&BlockInstance::GetBlockId, _1);
	SizeF getChrId = boost::bind(&BlockInstance::GetChrId, _1);
	SizeF getStart = boost::bind(&BlockInstance::GetStart, _1);
	const BlockComparer compareById = boost::bind(CompareBlocks<SizeF>, _1, _2, getId);
	const BlockComparer compareByChrId = boost::bind(CompareBlocks<SizeF>, _1, _2, getChrId);
	const BlockComparer compareByStart = boost::bind(CompareBlocks<SizeF>, _1, _2, getStart);

	int BlockInstance::GetSignedBlockId() const
	{
		return id_;
	}

	FastaRecord::Direction BlockInstance::GetDirection() const
	{
		return id_ > 0 ? FastaRecord::positive : FastaRecord::negative;
	}

	int BlockInstance::GetSign() const
	{
		return GetSignedBlockId() > 0 ? +1 : -1;
	}

	int BlockInstance::GetBlockId() const
	{
		return Abs(id_);
	}

	const FastaRecord& BlockInstance::GetChrInstance() const
	{
		return *chr_;
	}

	size_t BlockInstance::GetChrId() const
	{
		return GetChrInstance().GetId();
	}

	size_t BlockInstance::GetStart() const
	{
		return start_;
	}

	size_t BlockInstance::GetEnd() const
	{
		return end_;
	}

	size_t BlockInstance::GetConventionalStart() const
	{
		if(GetDirection() == FastaRecord::positive)
		{
			return start_ + 1;
		}

		return end_;
	}

	size_t BlockInstance::GetConventionalEnd() const
	{
		if(GetDirection() == FastaRecord::positive)
		{
			return end_;
		}

		return start_ + 1;
	}

	std::pair<size_t, size_t> BlockInstance::CalculateOverlap(const BlockInstance & instance) const
	{
		if(GetChrId() == instance.GetChrId())
		{			
			size_t overlap = 0;
			if(GetStart() >= instance.GetStart() && GetStart() <= instance.GetEnd())
			{
				return std::pair<size_t, size_t>(GetStart(), std::min(GetEnd(), instance.GetEnd()));
			}

			if(instance.GetStart() >= GetStart() && instance.GetStart() <= GetEnd())
			{
				return std::pair<size_t, size_t>(instance.GetStart(), std::min(GetEnd(), instance.GetEnd()));
			}
		}

		return std::pair<size_t, size_t>(0, 0);
	}

	bool BlockInstance::operator == (const BlockInstance & toCompare) const
	{
		return start_ == toCompare.start_ && end_ == toCompare.end_ && chr_->GetId() == toCompare.chr_->GetId() && id_ == toCompare.id_;
	}

	bool BlockInstance::operator != (const BlockInstance & toCompare) const
	{
		return !(*this == toCompare);
	}

	void BlockInstance::Reverse()
	{
		id_ = -id_;
	}

	size_t BlockInstance::GetLength() const
	{
		return end_ - start_;
	}

	bool BlockInstance::operator < (const BlockInstance & toCompare) const
	{
		return std::make_pair(GetChrInstance().GetId(), GetStart()) < std::make_pair(toCompare.GetChrInstance().GetId(), toCompare.GetStart());
	}
}