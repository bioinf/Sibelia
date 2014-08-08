//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "debruijnindex.h"

namespace SyntenyFinder
{
	const size_t DeBruijnIndex::NO_BIFURCATION = -1;

	size_t DeBruijnIndex::GetStrand(FastaRecord::Direction dir)
	{
		return dir == FastaRecord::positive ? 0 : 1;
	}

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

	DeBruijnIndex::DeBruijnIndex(const std::vector<FastaRecord> * chr, size_t bifNumber):
		chr_(chr), bifurcationPosition_(bifNumber)
	{
		positionEdge_[0].resize(chr->size());
		positionEdge_[1].resize(chr->size());
	}

	void DeBruijnIndex::AddEdge(size_t chrId, size_t pos, FastaRecord::Direction dir, size_t bifId, char mark)
	{		
		positionEdge_[GetStrand(dir)][chrId].insert(std::make_pair(static_cast<uint32_t>(pos), EdgeData(bifId, mark)));
		bifurcationPosition_[bifId].push_back(Location(chrId, pos, dir));
	}

	void  DeBruijnIndex::RemoveEdge(Edge edge, FastaRecord::Direction dir)
	{
	
	}
	
	DeBruijnIndex::Edge DeBruijnIndex::GetEdgeAtPosition(FastaRecord::Iterator it) const
	{
		EdgeData data = positionEdge_[GetStrand(it.GetDirection())][it.GetSequence()->GetId()].find(static_cast<uint32_t>(it.GetPosition()))->second;
		Location location(it.GetSequence()->GetId(), it.GetPosition(), it.GetDirection());
		return Edge(data, location);
	}
	
	void DeBruijnIndex::GetEdgesOfVertex(size_t bifId, std::vector<Edge> & e) const
	{
		e.clear();
		for(size_t i = 0; i < bifurcationPosition_[bifId].size(); i++)
		{
			Location location = bifurcationPosition_[bifId][i];
			FastaRecord::Iterator it = (*chr_)[location.GetChromosomeId()].Begin(location.GetDirection());
			e.push_back(GetEdgeAtPosition(it + location.GetPosition()));
		}
	}
}
