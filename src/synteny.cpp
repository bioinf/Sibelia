//****************************************************************************
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
	
	bool BlockFinder::TrimBlocks(std::vector<Edge> & block, size_t trimK, size_t minSize)
	{	
		size_t pos = 0;
		bool drop = false;		
		std::vector<std::string> blockSeq(block.size());		
		for(size_t i = 0; i < block.size(); i++)
		{
			std::string::const_iterator begin = (*originalChrList_)[block[i].GetChr()].GetSequence().begin();
			blockSeq[i].assign(begin + block[i].GetOriginalPosition(), begin + block[i].GetOriginalPosition() + block[i].GetOriginalLength());
		}
		
		const size_t oo = UINT_MAX;
		IteratorProxyVector startKMer;
		IndexedSequence iseq(blockSeq, trimK, "");
		DNASequence & sequence = iseq.Sequence();
		BifurcationStorage & bifStorage = iseq.BifStorage();		
		iseq.ConstructChrIndex();
		std::vector<Edge> ret;
		for(size_t chr = 0; chr < block.size(); chr++)
		{
			StrandIterator begin = sequence.Begin(block[chr].GetDirection(), chr);
			StrandIterator end = sequence.End(block[chr].GetDirection(), chr);
			StrandIterator trimStart = end;
			StrandIterator trimEnd = end;
			StrandIterator otherTrimStart;
			StrandIterator otherTrimEnd;
			size_t minBifStart = oo;
			size_t minBifEnd = oo;
			size_t minStartSum = oo;
			size_t minEndSum = oo;
			for(StrandIterator it = begin; it != end; ++it)
			{												
				size_t bifId = bifStorage.GetBifurcation(it);
				if(bifId != BifurcationStorage::NO_BIFURCATION)
				{
					startKMer.clear();
					bifStorage.ListPositions(bifId, std::back_inserter(startKMer));
					for(size_t pos = 0; pos < startKMer.size(); pos++)
					{                                                  
						StrandIterator kmer = *startKMer[pos];
						size_t kmerChr = iseq.GetChr(kmer);
						if(chr != kmerChr)
						{
							StrandIterator kmerChrStart = sequence.Begin(block[kmerChr].GetDirection(), kmerChr);
							StrandIterator kmerChrEnd = sequence.End(block[kmerChr].GetDirection(), kmerChr);
							size_t kmerStartDist = IndexedSequence::StrandIteratorDistance(kmer, kmerChrStart);
							size_t kmerEndDist = IndexedSequence::StrandIteratorDistance(kmer, AdvanceBackward(kmerChrEnd, 1));
							size_t itStartDist = IndexedSequence::StrandIteratorDistance(it, begin);
							size_t itEndDist = IndexedSequence::StrandIteratorDistance(it, AdvanceBackward(end, 1));
							size_t nowStartSum = kmerStartDist + itStartDist;
							size_t nowEndSum = kmerEndDist + itEndDist;
							if(nowStartSum < minStartSum || (nowStartSum == minStartSum && bifId < minBifStart))
							{
								minBifStart = bifId;
								minStartSum = nowStartSum;
								trimStart = it;
								otherTrimStart = kmer;
							}              
							
							if(nowEndSum < minEndSum || (nowEndSum == minEndSum && bifId < minBifEnd))
							{
								minBifEnd = bifId;
								minEndSum = nowEndSum;
								trimEnd = it;
								otherTrimEnd = kmer;
							}
						}
					}
				}
			}
			
			if(minStartSum < oo && minEndSum < oo)
			{				
				size_t size = IndexedSequence::StrandIteratorDistance(trimStart, trimEnd) + trimK;
				if(size >= minSize)
				{	
					std::advance(trimEnd, trimK - 1);
					size_t start = block[chr].GetOriginalPosition() + std::min(trimStart.GetOriginalPosition(), trimEnd.GetOriginalPosition());
					size_t end = block[chr].GetOriginalPosition() + std::max(trimStart.GetOriginalPosition(), trimEnd.GetOriginalPosition()) + 1;				
					ret.push_back(Edge(block[chr].GetChr(), block[chr].GetDirection(), block[chr].GetStartVertex(), block[chr].GetEndVertex(),
						block[chr].GetActualPosition(), block[chr].GetActualLength(), start, end - start, block[chr].GetFirstChar()));					
				}
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
			IndexedSequence iseq(rawSeq_, originalPos_, k, tempDir_);
			ListEdges(iseq.Sequence(), iseq.BifStorage(), k, edge);
			visit.resize(iseq.Sequence().ChrNumber());
		}
		
		block.clear();
		int blockCount = 1;
		edge.erase(std::remove_if(edge.begin(), edge.end(), boost::bind(EdgeEmpty, _1, minSize)), edge.end());
		std::vector<std::pair<size_t, size_t> > group;
		GroupBy(edge, CompareEdgesNaturally, std::back_inserter(group));
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

		GlueStripes(block);
		std::sort(block.begin(), block.end(), CompareBlocksNaturally);
	}
}
