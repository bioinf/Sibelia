//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "blockfinder.h"

namespace SyntenyFinder
{
	//To be refactored
	void BlockFinder::TrimBlocks(std::vector<Edge> & block, size_t trimK, size_t minSize)
	{		
		size_t pos = 0;
		std::vector<size_t> chrStart(block.size());
		std::vector<std::string> blockSeq(block.size());
		std::vector<std::vector<Pos> > chrPos(block.size());
		for(size_t i = 0; i < block.size(); i++)
		{
			std::string::const_iterator begin = (*originalChrList_)[block[i].GetChr()].GetSequence().begin();
			blockSeq[i].assign(begin + block[i].GetOriginalPosition(), begin + block[i].GetOriginalPosition() + block[i].GetOriginalLength());
			for(size_t j = 0; j < blockSeq[i].size(); j++)
			{
				chrPos[i].push_back(static_cast<Pos>(pos++));
			}

			chrStart[i] = chrPos[i][0];
		}
		
		const size_t oo = UINT_MAX;
		IteratorProxyVector startKMer;
		std::auto_ptr<DNASequence> sequence;
		std::auto_ptr<BifurcationStorage> bifStorage;
		{
			std::vector<std::vector<BifurcationInstance> > bifurcation(2);	
			size_t maxId = EnumerateBifurcationsSArrayInRAM(blockSeq, trimK, bifurcation[0], bifurcation[1]);
			bifStorage.reset(new BifurcationStorage(maxId));
			sequence.reset(new DNASequence(blockSeq, chrPos, true));
			ConstructBifStorage(*sequence, bifurcation, *bifStorage);
		}

		std::vector<Edge> ret;
		for(size_t chr = 0; chr < block.size(); chr++)
		{
			StrandIterator begin = sequence->Begin(block[chr].GetDirection(), chr);
			StrandIterator end = sequence->End(block[chr].GetDirection(), chr);
			size_t trimStart = oo;
			size_t trimEnd = oo;
			size_t minTrimStartSum = oo;
			size_t minTrimEndSum = oo;
			for(size_t step = 0; begin != end; ++begin, ++step)
			{												
				size_t bifId = bifStorage->GetBifurcation(begin);
				if(bifId != BifurcationStorage::NO_BIFURCATION)
				{
					startKMer.clear();
					bifStorage->ListPositions(bifId, std::back_inserter(startKMer));
					for(size_t pos = 0; pos < startKMer.size(); pos++)
					{
						StrandIterator kmer = *startKMer[pos];
						size_t bifChr = std::upper_bound(chrStart.begin(), chrStart.end(), kmer.GetOriginalPosition()) - chrStart.begin() - 1;
						if(chr != bifChr)
						{
							trimStart = std::min(trimStart, step);
							trimEnd = std::max(trimEnd, step);
						}
					}
				}
			}

			trimEnd += trimK - 1;
			if(trimEnd - trimStart + 1 >= minSize)
			{
				if(block[chr].GetDirection() == DNASequence::negative)
				{
					size_t oldTrimEnd = trimEnd;
					trimEnd = blockSeq[chr].size() - 1 - trimStart;
					trimStart = blockSeq[chr].size() - 1 - oldTrimEnd;
				}

				size_t start = block[chr].GetOriginalPosition() + trimStart;
				size_t end = block[chr].GetOriginalPosition() + trimEnd + 1;
				ret.push_back(Edge(block[chr].GetChr(), block[chr].GetDirection(), block[chr].GetStartVertex(), block[chr].GetEndVertex(),
					block[chr].GetOriginalPosition(), block[chr].GetOriginalLength(), start, end - start, block[chr].GetFirstChar()));
			}			
		}

		block.swap(ret);
	}

	void BlockFinder::GenerateSyntenyBlocks(size_t k, size_t trimK, size_t minSize, std::vector<BlockInstance> & block, bool sharedOnly, ProgressCallBack enumeration)
	{
		std::vector<Edge> edge;
		std::vector<std::set<size_t> > visit;
		{
			std::auto_ptr<DNASequence> sequence;
			std::auto_ptr<BifurcationStorage> bifStorage;
			ConstructIndex(sequence, bifStorage, k);
			BlockFinder::ListEdges(*sequence, *bifStorage, k, edge);
			visit.resize(sequence->ChrNumber());
		}
		
		block.clear();
		int blockCount = 1;
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
			if(edge[firstEdge].GetDirection() != DNASequence::positive || lastEdge - firstEdge < 2)
			{
				continue;
			}
			
			bool hit = false;
			std::vector<Edge> nowBlock;
			std::vector<size_t> occur(rawSeq_.size(), 0);
			for(size_t nowEdge = firstEdge; nowEdge < lastEdge; nowEdge++)
			{
				occur[edge[nowEdge].GetChr()]++;				
				hit = hit || (visit[edge[nowEdge].GetChr()].count(edge[nowEdge].GetOriginalPosition()) != 0);
				if(std::find_if(nowBlock.begin(), nowBlock.end(), boost::bind(&Edge::Overlap, boost::cref(edge[nowEdge]), _1)) == nowBlock.end())
				{
					nowBlock.push_back(edge[nowEdge]);
				}
			}

			if(!hit)
			{
				TrimBlocks(nowBlock, trimK, minSize);
				if(nowBlock.size() > 1 && (!sharedOnly || std::count(occur.begin(), occur.end(), 1) == rawSeq_.size()))
				{
					for(size_t i = 0; i < nowBlock.size(); i++)
					{
						int strand = nowBlock[i].GetDirection() == DNASequence::positive ? +1 : -1;
						visit[nowBlock[i].GetChr()].insert(nowBlock[i].GetOriginalPosition());
						block.push_back(BlockInstance(blockCount * strand, &(*originalChrList_)[nowBlock[i].GetChr()], nowBlock[i].GetOriginalPosition(), nowBlock[i].GetOriginalPosition() + nowBlock[i].GetOriginalLength()));
					}

					blockCount++;
				}
			}
		}
		
		std::sort(block.begin(), block.end(), CompareBlocksNaturally);
	}
}
