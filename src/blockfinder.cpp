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
		const size_t PROGRESS_STRIDE = 50;
	}
	/*
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
			for(size_t id = 0; id <= bifStorage.GetMaxId(); id++)
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
	}*/
	

	/*
	size_t BlockFinder::PerformGraphSimplifications(size_t k, size_t minBranchSize, size_t maxIterations, ProgressCallBack f)
	{
		IndexedSequence iseq(rawSeq_, originalPos_, k, tempDir_, true);
		iseq_ = &iseq;
		DNASequence & sequence = iseq.Sequence();
		BifurcationStorage & bifStorage = iseq.BifStorage();		
		size_t ret = SimplifyGraph(sequence, bifStorage, k, minBranchSize, maxIterations, f);
		for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
		{
			originalPos_[chr].clear();
			rawSeq_[chr].clear();
			StrandIterator end = sequence.PositiveEnd(chr);
			for(StrandIterator it = sequence.PositiveBegin(chr); it != end; ++it)
			{
				rawSeq_[chr].push_back(*it);
				originalPos_[chr].push_back(static_cast<Pos>(it.GetOriginalPosition()));
			}
		}

		return ret;
	}*/	
}