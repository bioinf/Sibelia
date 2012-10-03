#ifndef _GRAPH_ALGORITHM_H
#define _GRAPH_ALGORITHM_H

#include "hashing.h"
#include "blockinstance.h"
#include "bifurcationstorage.h"

namespace SyntenyFinder
{
	struct VisitData
	{
		size_t kmerId;
		size_t distance;
		VisitData() {}
		VisitData(size_t kmerId, size_t distance): kmerId(kmerId), distance(distance) {}
	};

	

	typedef char Bool;	
	typedef DNASequence::StrandIterator StrandIterator;
	typedef std::pair<StrandIterator, StrandIterator> IteratorPair;
	typedef boost::unordered_multimap<size_t, VisitData> VertexVisitMap;
	typedef boost::unordered_multiset<StrandIterator, WindowHashFunction, KMerEqualTo> KMerMultiSet;	
	
	class GraphAlgorithm
	{
	public:
		typedef boost::function<void(size_t)> ProgressCallBack;
		static void SerializeGraph(const DNASequence & sequence, size_t k, std::ostream & out);
		static void SerializeCondensedGraph(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k, std::ostream & out);
		static void PrintRaw(const DNASequence & s, std::ostream & out);
		static void PrintPath(StrandIterator e, size_t k, size_t distance, std::ostream & out);
		static size_t EnumerateBifurcations(const DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, ProgressCallBack f = ProgressCallBack());
		static size_t SimplifyGraph(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, size_t minBranchSize, size_t maxIterations, ProgressCallBack f = ProgressCallBack());
		static void Test(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k);
		static void GenerateSyntenyBlocks(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k, std::vector<BlockInstance> & chrList);
	private:	
		struct Edge
		{
			size_t chr;
			DNASequence::Direction direction;
			size_t startVertex;
			size_t endVertex;
			size_t actualPosition;
			size_t actualLength;
			size_t originalPosition;
			size_t originalLength;
			char firstChar;
			Edge() {}
			Edge(size_t chr, DNASequence::Direction direction, size_t startVertex, size_t endVertex,
				size_t actualPosition, size_t actualLength, size_t originalPosition, size_t originalLength, char firstChar):
				chr(chr), direction(direction), startVertex(startVertex), endVertex(endVertex), 
					actualPosition(actualPosition), actualLength(actualLength), originalPosition(originalPosition), originalLength(originalLength), firstChar(firstChar) {}

			bool Coincide(const Edge & edge) const
			{
				return startVertex == edge.startVertex && endVertex == edge.endVertex && firstChar == edge.firstChar;
			}
		};

		static bool EdgeEmpty(const Edge & a, size_t k);
		static bool EdgeCompare(const Edge & a, const Edge & b);
		static std::vector<size_t> EdgeToVector(const Edge & a);
		static size_t RemoveBulges(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, size_t minBranchSize, size_t bifId);		
		static void ListEdges(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k, std::vector<Edge> & edge);
		
	};
}

#endif
