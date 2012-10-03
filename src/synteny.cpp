#include "blockfinder.h"

namespace SyntenyFinder
{
	std::vector<size_t> BlockFinder::EdgeToVector(const Edge & a)
	{
		size_t feature[] = {a.startVertex, a.endVertex, a.firstChar, a.direction, a.chr};
		return std::vector<size_t>(feature, feature + sizeof(feature) / sizeof(feature[0]));
	}

	bool BlockFinder::EdgeEmpty(const Edge & a, size_t k)
	{
		return a.originalLength < k;
	}

	bool BlockFinder::EdgeCompare(const Edge & a, const Edge & b)
	{
		return EdgeToVector(a) < EdgeToVector(b);
	}

	void BlockFinder::ConvertEdgesToBlocks(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k, std::vector<BlockInstance> & block) const
	{
		int blockCount = 1;
		block.clear();
		std::vector<Edge> edge;
		BlockFinder::ListEdges(sequence, bifStorage, k, edge);
		std::vector<std::set<std::pair<size_t, size_t> > > visit(sequence.ChrNumber());
		edge.erase(std::remove_if(edge.begin(), edge.end(), boost::bind(EdgeEmpty, _1, k)), edge.end());
		std::sort(edge.begin(), edge.end(), EdgeCompare);

 		for(size_t now = 0; now < edge.size(); )
		{
			bool hit = false;
			size_t prev = now;
			for(; now < edge.size() && edge[prev].Coincide(edge[now]); now++)
			{
				std::pair<size_t, size_t> coord(edge[now].originalPosition, edge[now].originalLength);
				hit = hit || (visit[edge[now].chr].count(coord) != 0);
			}

			if(!hit && edge[prev].direction == DNASequence::positive && now - prev > 1)
			{
				for(; prev < now; prev++)
				{
					int strand = edge[prev].direction == DNASequence::positive ? +1 : -1;
					visit[edge[prev].chr].insert(std::make_pair(edge[prev].originalPosition, edge[prev].originalLength));
					block.push_back(BlockInstance(blockCount * strand, edge[prev].chr, edge[prev].originalPosition, edge[prev].originalPosition + edge[prev].originalLength));
				}

				blockCount++;
			}
		}

		std::sort(block.begin(), block.end(), CompareBlocksNaturally);
	}
}
