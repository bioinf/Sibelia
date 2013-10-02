//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _BLOCK_FINDER_H_
#define _BLOCK_FINDER_H_

#include "blockinstance.h"
#include "indexedsequence.h"

namespace SyntenyFinder
{
	#define NEW_ENUMERATION

	struct VisitData
	{
		size_t kmerId;
		size_t distance;
		VisitData() {}
		VisitData(size_t kmerId, size_t distance): kmerId(kmerId), distance(distance) {}
	};

	typedef char Bool;	
	typedef std::vector<BifurcationStorage::IteratorProxy> IteratorProxyVector;
	
	class BlockFinder
	{
	public:
		enum State
		{
			start,
			run,
			end
		};

		static const char SEPARATION_CHAR;
		typedef boost::function<void(size_t, State)> ProgressCallBack;
		BlockFinder(const std::vector<FASTARecord> & chrList);
		BlockFinder(const std::vector<FASTARecord> & chrList, const std::string & tempDir);
		void SerializeGraph(size_t k, std::ostream & out);
		void SerializeCondensedGraph(size_t k, std::ostream & out, ProgressCallBack f = ProgressCallBack());
		void GenerateSyntenyBlocks(size_t k, size_t trimK, size_t minSize, std::vector<BlockInstance> & block, bool sharedOnly = false, ProgressCallBack f = ProgressCallBack());
		void PerformGraphSimplifications(size_t k, size_t minBranchSize, size_t maxIterations, ProgressCallBack f = ProgressCallBack());
	private:
		DISALLOW_COPY_AND_ASSIGN(BlockFinder);
		typedef std::vector<Pos> PosVector;
		typedef std::pair<size_t, size_t> ChrPos;
		std::string tempDir_;
		IndexedSequence * iseq_;		
		std::vector<std::string> rawSeq_;	
		std::vector<size_t> originalSize_;
		std::vector<PosVector> originalPos_;		
		const std::vector<FASTARecord> * originalChrList_;
		static const char POS_FREE;
		static const char POS_OCCUPIED;

		struct Edge
		{
		public:			
			Edge() {}
			Edge(size_t chr, DNASequence::Direction direction, size_t startVertex, size_t endVertex, size_t actualPosition, size_t actualLength, size_t originalPosition, size_t originalLength, char firstChar);
			size_t GetChr() const;
			DNASequence::Direction GetDirection() const;
			size_t GetStartVertex() const;
			size_t GetEndVertex() const;
			size_t GetActualPosition() const;
			size_t GetActualLength() const;
			size_t GetOriginalPosition() const;
			size_t GetOriginalLength() const;
			char GetFirstChar() const;
			static bool PositiveEdge(const Edge & e);
			bool operator < (const Edge & e);
		private:
			size_t chr;
			DNASequence::Direction direction;
			size_t startVertex;
			size_t endVertex;
			size_t actualPosition;
			size_t actualLength;
			size_t originalPosition;
			size_t originalLength;
			char firstChar;
		};		

		struct EdgeGroupComparer
		{
		public:
			bool operator () (const std::pair<size_t, size_t> & range1, const std::pair<size_t, size_t> & range2) const
			{
				size_t size1 = range1.second - range1.first;
				size_t size2 = range2.second - range2.first;
				return size1 > size2;
			}

			EdgeGroupComparer(const std::vector<Edge> * edge): edge_(edge) {}
		private:
			const std::vector<Edge> * edge_;			
		};		

		static bool EdgeEmpty(const Edge & a, size_t k);				
		static std::vector<size_t> EdgeToVector(const Edge & a);
		static bool CompareEdgesNaturally(const Edge & a, const Edge & b);
		static bool CompareEdgesByDirection(const Edge & a, const Edge & b);
		void PrintRaw(const DNASequence & s, std::ostream & out);
		void PrintPath(const DNASequence & s, StrandIterator e, size_t k, size_t distance, std::ostream & out);
		void SpellBulges(const DNASequence & sequence, size_t k, size_t bifStart, size_t bifEnd, const std::vector<StrandIterator> & startKMer, const std::vector<VisitData> & visitData);
		
		void Init(const std::vector<FASTARecord> & chrList);		
		size_t RemoveBulges(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, size_t minBranchSize, size_t bifId);		
		void ListEdges(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k, std::vector<Edge> & edge) const;
		bool TrimBlocks(std::vector<Edge> & block, size_t trimK, size_t minSize);
		size_t SimplifyGraph(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, size_t minBranchSize, size_t maxIterations, ProgressCallBack f = ProgressCallBack());
		void CollapseBulgeGreedily(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, IteratorProxyVector & startKMer, VisitData sourceData, VisitData targetData);
		void UpdateBifurcations(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, const IteratorProxyVector & startKMer, VisitData sourceData, VisitData targetData,
			const std::vector<std::pair<size_t, size_t> > & lookForward, const std::vector<std::pair<size_t, size_t> > & lookBack);
		typedef std::vector<Bool> Indicator;
		typedef std::vector<Edge>::iterator EdgeIterator;
		
		void ResolveOverlap(EdgeIterator start, EdgeIterator end, size_t minSize, std::vector<Indicator> & overlap, std::vector<Edge> & nowBlock) const;
		void Extend(std::vector<StrandIterator> current, std::vector<size_t> & start, std::vector<size_t> & end, std::vector<Indicator> & overlap, std::set<ChrPos> & localOverlap, IndexedSequence & iseq, bool forward);
	};
}

#endif
