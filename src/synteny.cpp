//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wc++11-extensions"


namespace SyntenyFinder
{
	/*
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

	const char BlockFinder::POS_FREE = 0;
	const char BlockFinder::POS_OCCUPIED = 1;

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

	void BlockFinder::ResolveOverlap(EdgeIterator start, EdgeIterator end, size_t minSize, std::vector<Indicator> & overlap, std::vector<Edge> & nowBlock) const
	{
		nowBlock.clear();		
		std::set<ChrPos> localOverlap;
		for(; start != end; ++start)
		{
			size_t segEnd = 0;
			size_t bestStart = 0;
			size_t bestEnd = 0;			
			size_t chrNumber = start->GetChr();
			size_t end = start->GetOriginalPosition() + start->GetOriginalLength();
			for(size_t segStart = start->GetOriginalPosition(); segStart < end; segStart = segEnd)
			{
				for(segEnd = segStart; segEnd < end && overlap[chrNumber][segEnd] == POS_FREE; segEnd++)
				{
					ChrPos now(chrNumber, segEnd);
					if(localOverlap.count(now) > 0)
					{
						break;
					}
				}

				if(segEnd - segStart > bestEnd - bestStart)
				{
					bestStart = segStart;
					bestEnd = segEnd;
				}

				segEnd += segEnd == segStart ? 1 : 0;
			}

			if(bestEnd - bestStart >= minSize)
			{
				nowBlock.push_back(Edge(start->GetChr(), start->GetDirection(), start->GetStartVertex(), start->GetEndVertex(),
					start->GetActualPosition(), start->GetActualLength(), bestStart, bestEnd - bestStart, start->GetFirstChar()));
				for(size_t pos = bestStart; pos < bestEnd; pos++)
				{
					ChrPos now(chrNumber, pos);
					localOverlap.insert(now);
				}
			}
		}
	}

	namespace
	{
		const char EMPTY = ' ';	
		template<class First, class Second>
			bool LessFirst(const std::pair<First, Second> & a, const std::pair<First, Second> & b)
			{
				return a.first < b.first;
			}

		typedef std::pair<size_t, std::vector<StrandIterator> > IteratorEdgeGroup;
		bool CompareIteratorEdgeGroup(const IteratorEdgeGroup & a, const IteratorEdgeGroup & b)
		{
			return std::make_pair(a.second.size(), a.first) < std::make_pair(b.second.size(), b.first);
		}		
	}

	void BlockFinder::Extend(std::vector<StrandIterator> current, std::vector<size_t> & start, std::vector<size_t> & end, std::vector<Indicator> & overlap, std::set<ChrPos> & localOverlap, IndexedSequence & iseq, bool forward)
	{
		bool extend = true;
		for(size_t step = 0; extend; ++step)
		{
			char headingChar = EMPTY;			
			for(size_t i = 0; i < current.size() && extend; i++)
			{					
				extend = current[i].AtValidPosition();
				if(extend)
				{
					size_t chrId = iseq.GetChr(current[i]);
					if(headingChar == EMPTY)
					{
						headingChar = *current[i];
					}

					size_t nowPos = current[i].GetOriginalPosition();
					ChrPos nowChrPos(chrId, nowPos);
					extend = (headingChar == *current[i]) && (overlap[chrId][nowPos] == POS_FREE) && localOverlap.count(nowChrPos) == 0;
				}					
			}

			if(extend)
			{
				for(size_t i = 0; i < start.size(); i++)
				{
					size_t chrId = iseq.GetChr(current[i]);
					size_t nowPos = current[i].GetOriginalPosition();
					localOverlap.insert(ChrPos(chrId, nowPos));
					if(forward)
					{
						end[i] = nowPos;						
						++current[i];
					}
					else
					{
						start[i] = nowPos;
						--current[i];
					}
				}
			}
		}
	}
	
	void BlockFinder::GenerateSyntenyBlocks(size_t k, size_t trimK, size_t minSize, std::vector<BlockInstance> & block, bool sharedOnly, ProgressCallBack enumeration)
	{
		std::vector<Edge> edge;
		std::vector<Indicator> overlap(rawSeq_.size());
		for(size_t i = 0; i < rawSeq_.size(); i++)
		{
			overlap[i].resize(originalSize_[i], POS_FREE);
		}

		{
			IndexedSequence iseq(rawSeq_, originalPos_, k, tempDir_);
			ListEdges(iseq.Sequence(), iseq.BifStorage(), k, edge);			
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
			std::vector<Edge>::iterator firstEdge = edge.begin() + group[g].first;
			std::vector<Edge>::iterator lastEdge = edge.begin() + group[g].second;			 
			std::sort(firstEdge, lastEdge, CompareEdgesByDirection);
			if(lastEdge - firstEdge < 2 || std::find_if(firstEdge, lastEdge, Edge::PositiveEdge) == lastEdge)
			{
				continue;
			}
			
			std::vector<Edge> nowBlock;
			std::vector<size_t> occur(rawSeq_.size(), 0);
			ResolveOverlap(firstEdge, lastEdge, minSize, overlap, nowBlock);			
			while(TrimBlocks(nowBlock, trimK, minSize));
			for(size_t nowEdge = 0; nowEdge < nowBlock.size(); nowEdge++)
			{
				occur[nowBlock[nowEdge].GetChr()]++;				
			}
			
			if(nowBlock.size() > 1 && (!sharedOnly || std::count(occur.begin(), occur.end(), 1) == rawSeq_.size()))
			{
				for(size_t i = 0; i < nowBlock.size(); i++)
				{
					int strand = nowBlock[i].GetDirection() == FastaRecord::positive ? +1 : -1;
					size_t start = nowBlock[i].GetOriginalPosition();
					size_t end = start + nowBlock[i].GetOriginalLength();
					Indicator::iterator it = overlap[nowBlock[i].GetChr()].begin();
					std::fill(it + start, it + end, POS_OCCUPIED);
					block.push_back(BlockInstance(blockCount * strand, &(*originalChrList_)[nowBlock[i].GetChr()], start, end));
				}

				blockCount++;
			}
		}

		std::sort(block.begin(), block.end(), CompareBlocksNaturally);
	}*/
}
