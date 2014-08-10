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

	bool DeBruijnIndex::Location::EquivalentLocation(const Location & a, const Location & b)
	{
		return a.GetChromosomeId() == b.GetChromosomeId() && a.GetDirection() == b.GetDirection() && a.GetPosition() == b.GetPosition();
	}
	
	FastaRecord::Direction DeBruijnIndex::Location::GetDirection() const
	{
		return chrId_ > 0 ? FastaRecord::positive : FastaRecord::negative;
	}

	size_t DeBruijnIndex::EdgeData::GetBifurcationId() const
	{
		return bifId_;
	}

	size_t DeBruijnIndex::EdgeData::GetProjection() const
	{
		return projection_;
	}

	DeBruijnIndex::EdgeData::EdgeData(size_t pos):
		pos_(static_cast<uint32_t>(pos))
	{
	}

	DeBruijnIndex::EdgeData::EdgeData(size_t pos, size_t bifId, char mark, size_t projection):
		pos_(static_cast<uint32_t>(pos)), bifId_(static_cast<uint32_t>(bifId)), mark_(mark),
		projection_(static_cast<uint32_t>(projection))
	{
	}

	size_t DeBruijnIndex::EdgeDataKey::operator () (const EdgeData & data) const
	{
		return data.GetVirtualPosition();
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

	size_t DeBruijnIndex::EdgeData::GetVirtualPosition() const
	{
		return pos_;
	}

	void DeBruijnIndex::AddEdge(size_t chrId, size_t pos, FastaRecord::Direction dir, size_t bifId, char mark, size_t projection) 
	{		
		positionEdge_[GetStrand(dir)][chrId].insert(EdgeData(pos, bifId, mark, projection));
		bifurcationPosition_[bifId].push_back(Location(chrId, pos, dir));
	}

	void DeBruijnIndex::RemoveEdge(Edge edge, FastaRecord::Direction dir)
	{
		LocationVector & v = bifurcationPosition_[edge.GetBifurcationId()];
		positionEdge_[GetStrand(edge.GetDirection())][edge.GetChromosomeId()].erase(static_cast<uint32_t>(edge.GetPosition())) ;
		Location location(edge.GetChromosomeId(), edge.GetPosition(), edge.GetDirection());
		LocationVector::iterator it = std::find_if(v.begin(), v.end(), boost::bind(Location::EquivalentLocation, location, _1));
		v.erase(it);
	}
	
	DeBruijnIndex::Edge DeBruijnIndex::GetEdgeAtPosition(size_t chrId, size_t pos, FastaRecord::Direction dir) const
	{
		EdgeData lookUp(pos);
		EdgeData data = *positionEdge_[dir][chrId].find(lookUp);
		Location location(chrId, pos, dir);
		return Edge(data, location);
	}
	
	void DeBruijnIndex::GetEdgesOfVertex(size_t bifId, std::vector<Edge> & e) const
	{
		e.clear();
		for(size_t i = 0; i < bifurcationPosition_[bifId].size(); i++)
		{
			Location location = bifurcationPosition_[bifId][i];
			e.push_back(GetEdgeAtPosition(location.GetChromosomeId(), location.GetPosition(), location.GetDirection()));
		}
	}

	bool DeBruijnIndex::EdgeDataEquivalence::operator () (const EdgeData & a, const EdgeData & b) const
	{
		return a.GetVirtualPosition() == b.GetVirtualPosition();
	}
}
