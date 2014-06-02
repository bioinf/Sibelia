//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "blockfinder.h"

namespace SyntenyFinder
{
	const size_t BlockFinder::PROGRESS_STRIDE = 50;

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

		std::vector<size_t> segment;
		for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
		{
			size_t nowSegment = 0;
			std::set<size_t> passedBifId;			
			StrandIterator end = sequence.End(DNASequence::positive, chr);
			for(StrandIterator it = sequence.Begin(DNASequence::positive, chr); it != end; ++it)
			{
				size_t bifId = bifStorage.GetBifurcation(it);
				if(bifId != BifurcationStorage::NO_BIFURCATION)
				{
					if(passedBifId.count(bifId) > 0)
					{
						segment.push_back(nowSegment);
						nowSegment = 0;
					}

					passedBifId.insert(bifId);
				}

				nowSegment++;
			}

			segment.push_back(nowSegment);
		}

		std::cout << "Segments:" << std::endl;
		std::copy(segment.begin(), segment.end(), std::ostream_iterator<size_t>(std::cout, " "));

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
	}
	
	BlockFinder::BlockFinder(const std::vector<FASTARecord> & chrList):
		originalChrList_(&chrList)
	{
		Init(chrList);
	}

	BlockFinder::BlockFinder(const std::vector<FASTARecord> & chrList, const std::string & tempDir):
		originalChrList_(&chrList), tempDir_(tempDir)
	{
		Init(chrList);
	}

	void BlockFinder::Init(const std::vector<FASTARecord> & chrList)
	{
		rawSeq_.resize(chrList.size());
		originalPos_.resize(chrList.size());
		for(size_t i = 0; i < originalPos_.size(); i++)
		{
			rawSeq_[i] = chrList[i].GetSequence();
			originalSize_.push_back(rawSeq_[i].size());
			originalPos_[i].resize(chrList[i].GetSequence().size());
			std::generate(originalPos_[i].begin(), originalPos_[i].end(), Counter<Pos>());
		}
	}

	size_t BlockFinder::PerformGraphSimplifications(size_t k, size_t minBranchSize, size_t minPathLength, size_t maxIterations, ProgressCallBack f, size_t model, bool easy)
	{
		IndexedSequence iseq(rawSeq_, originalPos_, k, tempDir_, true, model);
		iseq_ = &iseq;
		DNASequence & sequence = iseq.Sequence();
		BifurcationStorage & bifStorage = iseq.BifStorage();
		size_t bulges = 0;
		if(easy)
		{
			bulges = SimplifyGraphEasily(sequence, bifStorage, k, minBranchSize, minPathLength, maxIterations, f);
		}
		else
		{
			bulges = SimplifyGraph(sequence, bifStorage, k, minBranchSize, maxIterations, f);
		}

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

		return bulges;
	}	
}