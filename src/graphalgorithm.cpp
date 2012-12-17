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
		const size_t PROGRESS_STRIDE = 50;
	}
	
#ifdef _DEBUG	
	size_t BlockFinder::GetMustBeBifurcation(StrandIterator jt, size_t k)
	{		
		std::string buf(std::string(jt, AdvanceForward(jt, k)));
		KMerBifMap::iterator kt = idMap.find(buf);
		return kt == idMap.end() ? BifurcationStorage::NO_BIFURCATION : kt->second;
	}

	void BlockFinder::Test(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k)
	{	
		for(size_t strand = 0; strand < 2; strand++)
		{
			for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
			{
				StrandIterator begin = sequence.Begin((DNASequence::Direction)strand, chr);
				StrandIterator end = sequence.End((DNASequence::Direction)strand, chr);
				SlidingWindow<StrandIterator> window(begin, end, k);
				for(; window.Valid(); window.Move())
				{
					StrandIterator jt = window.GetBegin();
					size_t pos = sequence.GlobalIndex(jt);
					std::string buf(std::string(jt, AdvanceForward(jt, k)));
					size_t actualBifurcation = bifStorage.GetBifurcation(jt);
					size_t mustBeBifurcation = GetMustBeBifurcation(jt, k);
					assert(actualBifurcation == mustBeBifurcation);
				}	
			}
		}
	}
#endif
	
	size_t BlockFinder::SimplifyGraph(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, size_t minBranchSize, size_t maxIterations, ProgressCallBack callBack)
	{
		size_t count = 0;
		size_t totalBulges = 0;
		size_t iterations = 0;
		size_t totalProgress = 0;
		bool anyChanges = true;
		if(!callBack.empty())
		{
			callBack(totalProgress, start);
		}

		size_t threshold = (bifStorage.GetMaxId() * maxIterations) / PROGRESS_STRIDE;
		do
		{
			iterations++;
			for(size_t id = 0; id < bifStorage.GetMaxId(); id++)
			{			
				totalBulges += RemoveBulges(sequence, bifStorage, k, minBranchSize, id);
				if(++count >= threshold && !callBack.empty())
				{
					count = 0;
					totalProgress = std::min(totalProgress + 1, PROGRESS_STRIDE);
					callBack(totalProgress, run);
				}
			}
		}
		while((totalBulges > 0) && iterations < maxIterations);

		if(!callBack.empty())
		{
			callBack(PROGRESS_STRIDE, end);
		}
		
		return totalBulges;
	}
}
