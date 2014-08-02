//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "debruijnindex.h"

namespace SyntenyFinder
{
	const size_t DeBruijnIndex::NO_BIFURCATION = -1;

	DeBruijnIndex::Location::Location(size_t chrId, size_t pos, FastaRecord::Direction dir):
		chrId_(dir == FastaRecord::positive ? static_cast<int32_t>(chrId) : -static_cast<int32_t>(chrId)),
		pos_(static_cast<uint32_t>(pos))
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

	size_t DeBruijnIndex::EdgeData::GetBifurcationId() const
	{
		return bifId_;
	}

	DeBruijnIndex::EdgeData::EdgeData(size_t bifId, char mark): bifId_(static_cast<uint32_t>(bifId)), mark_(mark)
	{
	}

	char DeBruijnIndex::EdgeData::GetMark() const
	{
		return mark_;
	}

	bool DeBruijnIndex::Edge::Valid() const
	{
		return GetBifurcationId() != NO_BIFURCATION;
	}

	DeBruijnIndex::Edge::Edge(EdgeData data, Location location):
		EdgeData(data), Location(location)
	{
	}

	DeBruijnIndex::DeBruijnIndex(size_t chrNumber, size_t bifNumber):
		bifurcationPosition_(bifNumber)
	{
		positionEdge_[0].resize(chrNumber);
		positionEdge_[1].resize(chrNumber);
	}

	void DeBruijnIndex::AddEdge(size_t chrId, size_t pos, FastaRecord::Direction dir, size_t bifId, char mark)
	{
		size_t strand = dir == FastaRecord::positive ? 0 : 1;
		positionEdge_[strand][chrId].insert(std::make_pair(static_cast<uint32_t>(pos), EdgeData(bifId, mark)));
		bifurcationPosition_[bifId].push_back(Location(chrId, pos, dir));
	}
}
