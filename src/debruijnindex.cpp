//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "debruijnindex.h"

namespace SyntenyFinder
{
	const size_t DeBruijnIndex::NO_BIFURCATION = -1;

	DeBruijnIndex::DeBruijnIndex(size_t chrNumber, size_t bifNumber):
		bifurcationPosition_(bifNumber)
	{
		positionEdge_[0].resize(chrNumber);
		positionEdge_[1].resize(chrNumber);
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

	DeBruijnIndex::EdgeData::EdgeData(size_t pos, char mark): pos_(static_cast<uint32_t>(pos)), mark_(mark)
	{
	}

	char DeBruijnIndex::EdgeData::GetMark() const
	{
		return mark_;
	}

	DeBruijnIndex::Edge::Edge(size_t bifId, EdgeData data, Location location):
		bifId_(bifId), EdgeData(data), Location(location)
	{
	}
	
}