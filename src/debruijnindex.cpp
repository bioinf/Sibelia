//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "debruijnindex.h"

namespace SyntenyFinder
{
	const uint32_t DeBruijnIndex::EdgeData::NO_POSITION = -1;
	const uint32_t DeBruijnIndex::EdgeData::NO_BIFURCATION = -1;

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

	DeBruijnIndex::EdgeData::EdgeData(): pos_(NO_POSITION), bifId_(NO_BIFURCATION)
	{
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

	bool DeBruijnIndex::EdgeData::Valid() const
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
		for(size_t strand = 0; strand < 2; strand++)
		{			
			positionEdge_[strand].resize(chrNumber);
			for(size_t chr = 0; chr < chrNumber; chr++)
			{
				positionEdge_[strand][chr].set_deleted_key(EdgeData());
			}
		}		
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
		Edge ret;
		EdgeData lookUp(pos);
		PositionEdgeMap::const_iterator it = positionEdge_[dir][chrId].find(lookUp);
		if(it != positionEdge_[dir][chrId].end())
		{
			ret = Edge(*it, Location(chrId, pos, dir));
		}
		
		return ret;
	}

	size_t DeBruijnIndex::CountEdges(size_t bifId) const
	{
		return bifurcationPosition_[bifId].size();
	}
	
	size_t DeBruijnIndex::GetEdgesOfVertex(size_t bifId, std::vector<Edge> & e) const
	{
		e.clear();
		for(size_t i = 0; i < bifurcationPosition_[bifId].size(); i++)
		{
			Location location = bifurcationPosition_[bifId][i];
			e.push_back(GetEdgeAtPosition(location.GetChromosomeId(), location.GetPosition(), location.GetDirection()));
		}

		return e.size();
	}

	bool DeBruijnIndex::EdgeDataEquivalence::operator () (const EdgeData & a, const EdgeData & b) const
	{
		return a.GetVirtualPosition() == b.GetVirtualPosition();
	}

	size_t DeBruijnIndex::GetBifurcationsNumber() const
	{
		return bifurcationPosition_.size();
	}
}
