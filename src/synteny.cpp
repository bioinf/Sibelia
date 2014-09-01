//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wc++11-extensions"

#include "blockbuilder.h"

namespace SyntenyFinder
{
	void BlockBuilder::ResolveOverlap(PotentialInstanceVector & v, size_t minSize, std::vector<std::vector<bool> > & overlap) const
	{
		typedef std::pair<size_t, size_t> ChrPos;
		PotentialInstanceVector ret;
		std::set<ChrPos> localOverlap;
		for(size_t start = 0; start < v.size(); ++start)
		{
			size_t segEnd = 0;
			size_t bestStart = 0;
			size_t bestEnd = 0;			
			size_t chrNumber = v[start].chr;
			size_t fullEnd = v[start].end;
			for(size_t segStart = v[start].end; segStart < end; segStart = segEnd)
			{
				for(segEnd = segStart; segEnd < end && overlap[chrNumber][segEnd] == false; segEnd++)
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
				ret.push_back(PotentialBlockInstance(v[start].chr, bestStart, bestEnd, v[start].dir));
				for(size_t pos = bestStart; pos < bestEnd; pos++)
				{
					ChrPos now(chrNumber, pos);
					localOverlap.insert(now);
				}
			}
		}

		v = ret;
	}
	

	bool BlockBuilder::CmpPotVector(const PotentialInstanceVector & a, const PotentialInstanceVector & b)
	{
		if(a.size() != b.size())
		{
			return a.size() > b.size();
		}

		size_t ac = 0;
		size_t bc = 0;
		for(size_t i = 0; i < a.size(); i++)
		{
			ac += 1 ? a[i].dir == FastaRecord::positive : 0;
			bc += 1 ? b[i].dir == FastaRecord::positive : 0;
		}

		return ac > bc;
	}

	bool BlockBuilder::CmpBifIt(const DeBruijnIndex::BifurcationIterator & a, const DeBruijnIndex::BifurcationIterator & b)
	{
		return a.GetOutMark() < b.GetOutMark();
	}
	
	void BlockBuilder::GenerateSyntenyBlocks(size_t k, size_t trimK, size_t minSize, std::vector<BlockInstance> & block) const
	{		
		std::vector<PotentialInstanceVector> potBlock;
		for(size_t bifId = 0; bifId < index_->GetBifurcationsNumber(); bifId++)
		{
			std::vector<std::pair<size_t, size_t> > group;
			std::vector<DeBruijnIndex::BifurcationIterator> it;
			index_->GetBifurcationInstances(bifId, it);
			GroupBy(it, CmpBifIt, std::back_inserter(group));
			for(size_t g = 0; g < group.size(); g++)
			{
				size_t positive = 0;
				PotentialInstanceVector v;
				for(size_t i = group[g].first; i < group[g].second; i++)
				{
					DeBruijnIndex::BifurcationIterator jt = it[i];
					if(FastaRecord::IsDefiniteBase(jt.GetOutMark()) && (jt + 1).AtEnd())
					{
						v.push_back(PotentialBlockInstance(jt.GetChromosomeId(), jt.GetProjection(), (jt + 1).GetEndingProjection(), jt.GetStrand()));
						if(v.back().dir != FastaRecord::positive)
						{
							std::swap(v.back().start, v.back().end);
						}
						else
						{
							++positive;
						}

						++v.back().end;
					}
				}

				if(v.size() > 1 && positive > 0)
				{
					potBlock.push_back(v);
				}
			}
		}
	
		std::sort(potBlock.begin(), potBlock.end(), CmpPotVector);
		std::vector<std::vector<bool> > overlap(originalChr_->size());
		
		for(size_t i = 0; i < overlap.size(); i++)
		{
			overlap[i].resize((*originalChr_)[i].GetSequence().size(), false);
		}
		
		block.clear();
		int blockCount = 1;
		for(size_t g = 0; g < potBlock.size(); g++)
		{
			ResolveOverlap(potBlock[g], minSize, overlap);			
			if(potBlock[g].size() > 1)
			{
				for(size_t i = 0; i < potBlock[g].size(); i++)
				{
					int strand = potBlock[g][i].dir == FastaRecord::positive ? +1 : -1;
					size_t start = potBlock[g][i].start;
					size_t end = potBlock[g][i].end;
					std::vector<bool>::iterator it = overlap[potBlock[g][i].chr].begin();
					std::fill(it + start, it + end, true);
					block.push_back(BlockInstance(blockCount * strand, &(*originalChr_)[potBlock[g][i].chr], start, end));
				}

				blockCount++;
			}
		}

		std::sort(block.begin(), block.end(), CompareBlocksNaturally);
	}
}
