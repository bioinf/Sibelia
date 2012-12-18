//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "blockfinder.h"

namespace SyntenyFinder
{
	namespace
	{
		struct Stripe
		{
			int firstBlock;
			int secondBlock;
			Stripe() {}
			Stripe(int firstBlock, int secondBlock): firstBlock(firstBlock), secondBlock(secondBlock) {}
			bool operator < (const Stripe & toCompare) const
			{
				return firstBlock < toCompare.firstBlock;
			}
		};
	}

	void BlockFinder::GlueStripes(std::vector<BlockInstance> & block)
	{
		std::vector<std::vector<BlockInstance> > perm(originalChrList_->size());
		for(size_t i = 0; i < block.size(); i++)
		{
			perm[block[i].GetChrId()].push_back(block[i]);
		}

		for(size_t i = 0; i < perm.size(); i++)
		{
			std::sort(perm[i].begin(), perm[i].end(), compareByStart);
		}

		int sentinel = INT_MAX >> 1;
		bool glue = false;
		do
		{
			std::vector<Stripe> stripe;
			for(size_t chr = 0; chr < perm.size(); chr++)
			{
				for(size_t i = 0; i < perm[chr].size(); i++)
				{
					int bid = perm[chr][i].GetSignedBlockId();
					if(bid > 0)
					{
						int nextBid = i < perm[chr].size() - 1 ? perm[chr][i + 1].GetSignedBlockId() : sentinel;
						stripe.push_back(Stripe(bid, nextBid));
					}
					else
					{
						int prevBid = i > 0 ? perm[chr][i - 1].GetSignedBlockId() : -sentinel;
						stripe.push_back(Stripe(-bid, -prevBid));
					}
				}
			}

			size_t now = 0;
			size_t next = 0;
			std::sort(stripe.begin(), stripe.end());
			for(; now < stripe.size(); now = next)
			{
				glue = true;
				for(; next < stripe.size() && stripe[next].firstBlock == stripe[now].firstBlock; next++)
				{
					if(stripe[next].secondBlock != stripe[now].secondBlock || stripe[next].secondBlock == sentinel || Abs(stripe[next].secondBlock) == stripe[next].firstBlock)
					{
						glue = false;
					}
				}

				if(glue)
				{
					typedef std::vector<Stripe>::iterator It;
					std::pair<It, It> range = std::equal_range(stripe.begin(), stripe.end(), Stripe(Abs(stripe[now].secondBlock), 0));
					if(range.second - range.first != next - now)
					{
						glue = false;
					}
					else
					{
						break;
					}
				}
			}

			if(glue)
			{
				assert(next - now > 1);
				int glueBid = stripe[now].firstBlock;
				for(size_t chr = 0; chr < perm.size(); chr++)
				{
					for(size_t i = 0; i < perm[chr].size(); i++)
					{
						int bid = perm[chr][i].GetBlockId();
						if(bid == glueBid)
						{
							bid = perm[chr][i].GetSignedBlockId();
							if(bid > 0)
							{
								BlockInstance & a = perm[chr][i];
								BlockInstance & b = perm[chr][i + 1];
								a = BlockInstance(a.GetSignedBlockId(), &a.GetChrInstance(), a.GetStart(), b.GetEnd());
								perm[chr].erase(perm[chr].begin() + i + 1);
							}
							else
							{
								BlockInstance & a = perm[chr][--i];
								BlockInstance & b = perm[chr][i + 1];
								a = BlockInstance(b.GetSignedBlockId(), &a.GetChrInstance(), a.GetStart(), b.GetEnd());
								perm[chr].erase(perm[chr].begin() + i + 1);
							}
						}
					}
				}
			}
		}
		while(glue);

		block.clear();
		std::vector<int> oldId;
		for(size_t chr = 0; chr < perm.size(); chr++)
		{
			for(size_t i = 0; i < perm[chr].size(); i++)
			{
				block.push_back(perm[chr][i]);
				oldId.push_back(perm[chr][i].GetBlockId());
			}
		}

		std::sort(oldId.begin(), oldId.end());
		oldId.erase(std::unique(oldId.begin(), oldId.end()), oldId.end());
		for(std::vector<BlockInstance>::iterator it = block.begin(); it != block.end(); ++it)
		{
			int sign = it->GetSignedBlockId() > 0 ? +1 : -1;
			size_t newId = std::lower_bound(oldId.begin(), oldId.end(), it->GetBlockId()) - oldId.begin() + 1;
			*it = BlockInstance(static_cast<int>(newId) * sign, &it->GetChrInstance(), it->GetStart(), it->GetEnd());
		}
	}

	//To be refactored
	size_t chop = 0;
	bool BlockFinder::TrimBlocks(std::vector<Edge> & block, size_t trimK, size_t minSize)
	{				
		size_t pos = 0;
		bool drop = false;
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
			size_t trimEnd = 0;
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
			if(trimEnd - trimStart + 1 >= minSize && trimEnd != 0 && trimStart != oo)
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
			else
			{
				drop = true;
			}
		}

		block.swap(ret);
		return drop;
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
				if(std::find_if(nowBlock.begin(), nowBlock.end(), boost::bind(&Edge::Overlap, boost::cref(edge[nowEdge]), _1)) == nowBlock.end())
				{
					nowBlock.push_back(edge[nowEdge]);
				}
			}

			while(TrimBlocks(nowBlock, trimK, minSize));
			for(size_t nowEdge = 0; nowEdge < nowBlock.size(); nowEdge++)
			{
				occur[nowBlock[nowEdge].GetChr()]++;				
				hit = hit || (visit[nowBlock[nowEdge].GetChr()].count(nowBlock[nowEdge].GetOriginalPosition()) != 0);
			}

			if(!hit && nowBlock.size() > 1 && (!sharedOnly || std::count(occur.begin(), occur.end(), 1) == rawSeq_.size()))
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
		/*
		block.clear();
		int b1[] = {+13, -12, +5, +7, +10, -16, +8, +19, +11, +6, -17, +21, -9, +3, -4, +18, +1, +14, +20, +16, +15, -23, +2};
		int b2[] = {+13, -14, +5, +7, +4, +11, +6, +9, +21, +17, +3, -19, -8, +22, -10, +1, +12, +20, +16, +15, +23, -18, +22, +2};
		size_t s1 = sizeof(b1) / sizeof(b1[0]);
		size_t s2 = sizeof(b2) / sizeof(b2[0]);
		for(size_t i = 0; i < s1; i++)
		{
			block.push_back(BlockInstance(b1[i], &(*originalChrList_)[0], 0, 0));
		}

		for(size_t i = 0; i < s2; i++)
		{
			block.push_back(BlockInstance(b2[i], &(*originalChrList_)[1], 0, 0));
		}*/

		GlueStripes(block);
		std::sort(block.begin(), block.end(), CompareBlocksNaturally);
	}
}
