//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "blockfinder.h"

namespace SyntenyFinder
{
	void BlockFinder::GenerateSyntenyBlocks(size_t k, size_t minSize, std::vector<BlockInstance> & block, bool sharedOnly, ProgressCallBack enumeration)
	{
		std::auto_ptr<DNASequence> sequence;
		std::auto_ptr<BifurcationStorage> bifStorage;
		ConstructIndex(sequence, bifStorage, k);

		int blockCount = 1;
		block.clear();
		std::vector<Edge> edge;		
		BlockFinder::ListEdges(*sequence, *bifStorage, k, edge);
		std::vector<std::set<size_t> > visit(sequence->ChrNumber());
		edge.erase(std::remove_if(edge.begin(), edge.end(), boost::bind(EdgeEmpty, _1, minSize)), edge.end());
		std::vector<std::pair<size_t, size_t> > group;
		GroupBy(edge, EdgeCompare, std::back_inserter(group));
		EdgeGroupComparer groupComparer(&edge);
		std::sort(group.begin(), group.end(), groupComparer);

		for(size_t g = 0; g < group.size(); g++)
		{
			size_t firstEdge = group[g].first;
			size_t lastEdge = group[g].second;
			std::sort(edge.begin() + firstEdge, edge.begin() + lastEdge, CompareEdgesByDirection);
			if(edge[firstEdge].direction != DNASequence::positive || lastEdge - firstEdge < 2)
			{
				continue;
			}
			
			bool hit = false;
			std::vector<Edge> nowBlock;
			std::vector<size_t> occur(rawSeq_.size(), 0);
			for(size_t nowEdge = firstEdge; nowEdge < lastEdge; nowEdge++)
			{
				occur[edge[nowEdge].chr]++;				
				hit = hit || (visit[edge[nowEdge].chr].count(edge[nowEdge].originalPosition) != 0);
				if(std::find_if(nowBlock.begin(), nowBlock.end(), boost::bind(&Edge::Overlap, boost::cref(edge[nowEdge]), _1)) == nowBlock.end())
				{
					nowBlock.push_back(edge[nowEdge]);
				}
			}

			if(nowBlock.size() > 1 && !hit && (!sharedOnly || std::count(occur.begin(), occur.end(), 1) == rawSeq_.size()))
			{
				for(size_t i = 0; i < nowBlock.size(); i++)
				{
					int strand = nowBlock[i].direction == DNASequence::positive ? +1 : -1;
					visit[nowBlock[i].chr].insert(nowBlock[i].originalPosition);
					block.push_back(BlockInstance(blockCount * strand, &(*originalChrList_)[nowBlock[i].chr], nowBlock[i].originalPosition, nowBlock[i].originalPosition + nowBlock[i].originalLength));
				}

				blockCount++;
			}
		}
		
		std::sort(block.begin(), block.end(), CompareBlocksNaturally);
	}
}
