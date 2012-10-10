#ifndef _BLOCK_FINDER_H_
#define _BLOCK_FINDER_H_

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
	
	class BlockFinder
	{
	public:
		enum State
		{
			start,
			run,
			end
		};

		typedef boost::function<void(size_t, State)> ProgressCallBack;
		BlockFinder(const std::vector<FASTARecord> & chrList);
		void SerializeGraph(size_t k, std::ostream & out) const;
		void SerializeCondensedGraph(size_t k, std::ostream & out, ProgressCallBack f = ProgressCallBack()) const;
		void GenerateSyntenyBlocks(size_t k, std::vector<BlockInstance> & block, ProgressCallBack f = ProgressCallBack()) const;
		void PerformGraphSimplifications(size_t k, size_t minBranchSize, size_t maxIterations, ProgressCallBack f = ProgressCallBack());

		static void PrintRaw(const DNASequence & s, std::ostream & out);
		static void PrintPath(StrandIterator e, size_t k, size_t distance, std::ostream & out);
		static void Test(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k);
	private:	
		typedef std::vector<Pos> PosVector;
		std::vector<FASTARecord> chrList_;
		std::vector<PosVector> originalPos_;

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
			Edge(size_t chr, DNASequence::Direction direction, size_t startVertex, size_t endVertex, size_t actualPosition, size_t actualLength, size_t originalPosition, size_t originalLength, char firstChar);
			bool Coincide(const Edge & edge) const;
		};

		struct BifurcationInstance
		{
			Size bifId;
			Size chr;
			Size pos;
			BifurcationInstance() {}
			BifurcationInstance(Size bifId, Size chr, Size pos): bifId(bifId), chr(chr), pos(pos) {}
			bool operator < (const BifurcationInstance & toCompare)
			{
				return std::make_pair(chr, pos) < std::make_pair(toCompare.chr, toCompare.pos);
			}
		};

		static bool EdgeEmpty(const Edge & a, size_t k);
		static bool EdgeCompare(const Edge & a, const Edge & b);
		static std::vector<size_t> EdgeToVector(const Edge & a);	
				
		size_t RemoveBulges(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, size_t minBranchSize, size_t bifId);		
		void ListEdges(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k, std::vector<Edge> & edge) const;
		size_t EnumerateBifurcations(const DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, ProgressCallBack f = ProgressCallBack()) const;
		size_t EnumerateBifurcations(size_t k, std::vector<BifurcationInstance> & posBifurcation, std::vector<BifurcationInstance> & negBifurcation) const;
		void ConvertEdgesToBlocks(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k, std::vector<BlockInstance> & chrList) const;
		size_t SimplifyGraph(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, size_t minBranchSize, size_t maxIterations, ProgressCallBack f = ProgressCallBack());

	};
}

#endif
