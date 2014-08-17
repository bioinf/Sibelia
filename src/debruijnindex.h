//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "fasta.h"
#include "blockinstance.h"
#include "indexedsequence.h"

#ifndef _DE_BRUIJN_GRAPH_H_
#define _DE_BRUIJN_GRAPH_H_

namespace SyntenyFinder
{	
	class DeBruijnIndex
	{	
	public:		
		static const size_t MAX_POSITION;
		static const size_t MAX_BIFURCATION_ID;		
		static const size_t MAX_SEQUENCE_NUMBER;		

		class Edge;
		
		DeBruijnIndex(size_t chrNumber, size_t bifNumber);
		size_t GetVirtualChrSize(size_t chrId) const;
		void RemoveEdge(Edge edge, FastaRecord::Direction dir);
		void GetEdgesOfVertex(size_t bifId, std::vector<Edge> & e) const;
		Edge GetEdgeAtPosition(size_t chrId, size_t pos, FastaRecord::Direction dir) const;
		void AddEdge(size_t chrId, size_t pos, FastaRecord::Direction dir, size_t bifId, char mark, size_t projection);		
	private:
		DISALLOW_COPY_AND_ASSIGN(DeBruijnIndex);
		class EdgeData
		{
		public:
			EdgeData();
			bool Valid() const;
			EdgeData(size_t pos);			
			EdgeData(size_t pos, size_t bifId, char mark, size_t projection);
			char GetMark() const;			
			size_t GetProjection() const;
			size_t GetBifurcationId() const;
			size_t GetVirtualPosition() const;
		private:
			uint32_t pos_;			
			uint32_t bifId_;
			char mark_;
			uint32_t projection_;
			static const uint32_t NO_BIFURCATION;
		};

		class Location
		{
		public:
			Location() {}
			Location(size_t pos, size_t chrId, FastaRecord::Direction dir);
			size_t GetPosition() const;
			size_t GetChromosomeId() const;
			FastaRecord::Direction GetDirection() const;
			static bool EquivalentLocation(const Location & a, const Location & b);
		private:
			int32_t chrId_;
			uint32_t pos_;
		};

		class EdgeDataKey
		{
		public:
			size_t operator () (const EdgeData & data) const;
		};

		class EdgeDataEquivalence
		{
		public:
			bool operator () (const EdgeData & a, const EdgeData & b) const;
		};

		static size_t GetStrand(FastaRecord::Direction dir);

		typedef std::vector<Location> LocationVector;
		typedef boost::unordered_set<EdgeData, EdgeDataKey, EdgeDataEquivalence> PositionEdgeMap;
		const std::vector<FastaRecord> * chr_;
		std::vector<size_t> virtualChrSize_;
		std::vector<PositionEdgeMap> positionEdge_[2];
		std::vector<LocationVector> bifurcationPosition_;		
	public:
		class Edge: public EdgeData, public Location
		{
		public:
			Edge() {};			
		private:			
			friend class DeBruijnIndex;
			Edge(EdgeData data, Location location);
		};				
	};
}

#endif
