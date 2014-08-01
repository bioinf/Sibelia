//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "debruijnindex.h"

namespace SyntenyFinder
{
	const size_t DeBruijnIndex::NO_BIFURCATION = -1;
	DeBruijnIndex::DeBruijnIndex(const std::vector<FastaRecord> & chrList, size_t k, const std::string & tempDir)
	{
	}

	DeBruijnIndex::Location::Location(size_t pos, size_t chrId, FastaRecord::Direction dir):
		pos_(static_cast<uint32_t>(pos)),
		chrId_(dir == FastaRecord::positive ? static_cast<int32_t>(chrId) : -static_cast<int32_t>(chrId))
	{
	}

	size_t DeBruijnIndex::Location::GetPosition() const
	{
		return pos_;
	}
	
	size_t DeBruijnIndex::Location::GetChromosomeId() const
	{
		return abs(chrId_);
	}
	
	FastaRecord::Direction DeBruijnIndex::Location::GetDirection() const
	{
		return chrId_ > 0 ? FastaRecord::positive : FastaRecord::negative;
	}

	DeBruijnIndex::Edge::Edge(size_t pos, char mark): pos_(static_cast<uint32_t>(pos)), mark_(mark)
	{
	}

	char DeBruijnIndex::Edge::GetMark() const
	{
		return mark_;
	}

	size_t DeBruijnIndex::Edge::GetPositon() const
	{
		return pos_;
	}
}