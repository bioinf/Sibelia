#ifndef _GRAPH_ALGORITHM_H
#define _GRAPH_ALGORITHM_H

#include "hashing.h"
#include "bifurcationstorage.h"

namespace SyntenyBuilder
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
		static void FindGraphBulges(const DNASequence & sequence, BifurcationStorage & bifStorage, size_t k);
		static void SimplifyGraph(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, size_t minBranchSize);		
		static void SerializeGraph(const DNASequence & sequence, size_t k, std::ostream & out);
		static void SerializeCondensedGraph(const DNASequence & sequence, const BifurcationStorage & bifStorage,
			size_t k, std::ostream & out);
		static void ListNonBranchingPaths(const DNASequence & sequence, const BifurcationStorage & bifStorage,
			size_t k, std::ostream & out, std::ostream & indexOut);
		static void PrintRaw(const DNASequence & s, std::ostream & out);
		static void PrintPath(StrandIterator e, size_t k, size_t distance, std::ostream & out);
		static size_t EnumerateBifurcations(const DNASequence & sequence, BifurcationStorage & bifStorage, size_t k);
		static void Test(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k);
	private:				
	//	static size_t FindBulges(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k, size_t bifId);
	//	static size_t RemoveWhirls(BifurcationStorage & bifStorage, DNASequence & sequence, size_t k, size_t minBranchSize);
		static size_t RemoveBulges(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, size_t minBranchSize, size_t bifId);		
	};
}

#endif
