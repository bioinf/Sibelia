//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "blockinstance.h"
#include "indexedsequence.h"

#ifndef _DE_BRUIJN_GRAPH_H_
#define _DE_BRUIJN_GRAPH_H_

namespace SyntenyFinder
{	
	class DeBruijnIndex
	{	
	public:
		static const size_t NO_BIFURCATION;
		static const size_t MAX_POSITION;
		static const size_t MAX_BIFURCATION_ID;		
		static const size_t MAX_SEQUENCE_NUMBER;		

		class Edge;		

		Edge GetEdgeAtPosition(FastaRecord::Iterator it) const;
		DeBruijnIndex(size_t chrNumber, size_t bifNumber);
		void RemoveEdge(Edge edge, FastaRecord::Direction dir);
		void GetEdgesFromVertex(size_t bifId, std::vector<Edge> & pos) const;		
		void AddEdge(size_t chrId, size_t pos, FastaRecord::Direction dir, size_t bifId, char mark);		

	private:
		class EdgeData
		{
		public:
			EdgeData() {}
			EdgeData(size_t pos, char mark);
			char GetMark() const;
			size_t GetPositon() const;
		private:
			char mark_;
			uint32_t pos_;
		};

		class Location
		{
		public:
			Location() {}
			Location(size_t pos, size_t chrId, FastaRecord::Direction dir);
			size_t GetPosition() const;
			size_t GetChromosomeId() const;
			FastaRecord::Direction GetDirection() const;
		private:
			int32_t chrId_;
			uint32_t pos_;
		};

		typedef std::vector<Location> LocationVector;
		typedef boost::unordered_map<uint32_t, EdgeData> PositionEdgeMap;
		std::vector<PositionEdgeMap> positionEdge_[2];
		std::vector<LocationVector> bifurcationPosition_;

	public:

		class Edge: public EdgeData, public Location
		{
		public:
			Edge() {};
			size_t GetBifurcationId() const;
		private:			
			friend class DeBruijnIndex;
			Edge(size_t bifId, EdgeData data, Location location);
			size_t bifId_;
		};
		
	};

	void SimplifyGraph(DeBruijnIndex & g, size_t minBranchSize);	
	std::vector<BlockInstance> GenerateSyntenyBlocks(const DeBruijnIndex & g);
}

#endif
