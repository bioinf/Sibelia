#include "graphalgorithm.h"

namespace SyntenyBuilder
{
	std::vector<size_t> GraphAlgorithm::EdgeToVector(const Edge & a)
	{
		size_t feature[] = {a.startVertex, a.endVertex, a.firstChar, a.direction, a.chr};
		return std::vector<size_t>(feature, feature + sizeof(feature) / sizeof(feature[0]));
	}

	bool GraphAlgorithm::EdgeEmpty(const Edge & a)
	{
		return a.originalLength == 0;
	}

	bool GraphAlgorithm::EdgeCompare(const Edge & a, const Edge & b)
	{
		return EdgeToVector(a) < EdgeToVector(b);
	}

	void GraphAlgorithm::GenerateSyntenyBlocks(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k, std::vector<std::vector<BlockInstance> > & block)
	{
		int blockCount = 1;
		std::vector<Edge> edge;
		block.resize(sequence.ChrNumber(), std::vector<BlockInstance>());
		GraphAlgorithm::ListEdges(sequence, bifStorage, k, edge);
		edge.erase(std::remove_if(edge.begin(), edge.end(), EdgeEmpty), edge.end());
		std::sort(edge.begin(), edge.end(), EdgeCompare);
		std::vector<std::set<std::pair<size_t, size_t> > > visit(sequence.ChrNumber());

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
					visit[edge[prev].chr].insert(std::make_pair(edge[prev].originalPosition, edge[prev].originalLength));
					block[edge[prev].chr].push_back(BlockInstance(blockCount, edge[prev].chr, edge[prev].originalPosition, edge[prev].originalPosition + edge[prev].originalLength));
				}

				blockCount++;
			}
		}

		for(size_t i = 0; i < block.size(); i++)
		{
			std::sort(block[i].begin(), block[i].end());
		}
	}
}
