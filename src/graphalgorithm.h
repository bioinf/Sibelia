#ifndef _GRAPH_ALGORITHM_H
#define _GRAPH_ALGORITHM_H

#include "kmerindex.h"
#include "bifurcationstorage.h"

namespace SyntenyBuilder
{
	typedef DNASequence::StrandIterator StrandIterator;
	typedef char Bool;	
	typedef std::pair<StrandIterator, StrandIterator> IteratorPair;
	typedef google::sparse_hash_set<StrandIterator, KMerIndex::KMerHashFunction, KMerIndex::KMerEqualTo> KMerSet;	

	struct VisitData
	{
		size_t kmerId;
		size_t distance;
		VisitData() {}
		VisitData(size_t kmerId, size_t distance): kmerId(kmerId), distance(distance) {}
	};
	
	typedef boost::unordered_multimap<size_t, VisitData> VertexVisitMap;

	class GraphAlgorithm
	{
	public:		
		static void SimplifyGraph(DNASequence & sequence, size_t k, size_t minBranchSize);		
		static void SerializeGraph(const DNASequence & seq, size_t k, std::ostream & out);
		static void ListNonBranchingPaths(const DNASequence & sequence, size_t k, std::ostream & out, std::ostream & indexOut);
		static void PrintRaw(DNASequence & s, std::ostream & out);
		static void PrintPath(StrandIterator e, size_t k, size_t distance, std::ostream & out);
	private:		
		static size_t EnumerateBifurcations(DNASequence & sequence, size_t k, BifurcationStorage & bifStorage);
		static size_t RemoveWhirls(BifurcationStorage & bifStorage, DNASequence & sequence, size_t k, size_t minBranchSize);
		static size_t RemoveBulges(BifurcationStorage & bifStorage, DNASequence & sequence, size_t k, size_t minBranchSize, size_t bifId);		
	};
}

#endif
