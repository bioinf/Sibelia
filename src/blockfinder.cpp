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

	void BlockFinder::PerformGraphSimplifications(size_t k, size_t minBranchSize, size_t maxIterations, ProgressCallBack f, size_t model)
	{
		IndexedSequence iseq(rawSeq_, originalPos_, k, tempDir_, true, model);
		if(model != IndexedSequence::NO_MODEL)
		{
			std::ofstream out("out/kgraph.dot");
			out << "digraph G" << std::endl << "{" << std::endl;
			out << "rankdir=LR" << std::endl;
			std::vector<Edge> edge;
			ListEdges(iseq.Sequence(), iseq.BifStorage(), k, edge);
			for(size_t i = 0; i < edge.size(); i++)
			{
				char buf[1 << 8];
				std::string color = edge[i].GetDirection() == DNASequence::positive ? "blue" : "red";
				int uchr = static_cast<int>(edge[i].GetChr());
				int uorpos = static_cast<int>(edge[i].GetOriginalPosition());
				int uorlength = static_cast<int>(edge[i].GetOriginalLength());
				int upos = static_cast<int>(edge[i].GetActualPosition());
				int ulength = static_cast<int>(edge[i].GetActualLength());
				out << edge[i].GetStartVertex() << " -> " << edge[i].GetEndVertex();
				sprintf(&buf[0], "[color=\"%s\", label=\"chr=%i pos=%i len=%i orpos=%i orlen=%i  ch='%c'\"];", color.c_str(), uchr, upos, ulength, uorpos, uorlength, edge[i].GetFirstChar());
				out << " " << buf << std::endl;
			}

			out << "}" << std::endl;
		}

		iseq_ = &iseq;
		DNASequence & sequence = iseq.Sequence();
		BifurcationStorage & bifStorage = iseq.BifStorage();		
		SimplifyGraph(sequence, bifStorage, k, minBranchSize, maxIterations, f);
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

	}	
}