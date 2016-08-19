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

	struct Range
	{
		size_t chr;
		size_t first;
		size_t second;
		Range() {}
		Range(size_t chr, size_t first, size_t second): chr(chr), first(first), second(second){}
	};

	static bool Intersect(Range r1, Range r2)
 	{
		if(r1.chr == r2.chr)
		{
 			size_t a = r1.first;
 			size_t b = r1.second;
 			size_t c = r2.first;
 			size_t d = r2.second;			
 			if ((a >= c && b < d) || (c >= a && d < b) || (a >= c && a < d) || (b >= c && b < d))
 			{
				return true;
			}
		}

 		return false;
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

		std::vector<size_t> interest;
		std::string buf;
  		std::vector<Range> interestRange;
  		std::ifstream interestFile("investigate/investigate.maf");
		if(!interestFile)
		{
			throw std::exception("Can't open the interest file");
		}

 		while (std::getline(interestFile, buf))
 		{
 			std::stringstream ss;
 			ss << buf;
 			std::vector<std::string> part;
 			while (ss >> buf)
 			{
 				part.push_back(buf);
 			}

			Range rng;
			rng.chr = atoi(&(*(part[1].begin() + 2))) - 1;
			rng.first = atoi(part[2].c_str());
			rng.second = rng.first + atoi(part[3].c_str());
			if(part[4] == "-")
			{
				size_t size = atoi(part[5].c_str());
				rng = Range(rng.chr, size - rng.second, size - rng.first);
			}

			interestRange.push_back(rng);
 		}

		std::vector<Edge> edge;
		ListEdges(sequence, bifStorage, k, edge);
		for(size_t i = 0; i < edge.size(); i++)
		{
			int uchr = static_cast<int>(edge[i].GetChr());
			int uorpos = static_cast<int>(edge[i].GetOriginalPosition());
			int uorlength = static_cast<int>(edge[i].GetOriginalLength());
			int upos = static_cast<int>(edge[i].GetActualPosition());
			int ulength = static_cast<int>(edge[i].GetActualLength());
			Range edgeRange(uchr, uorpos, uorpos + uorlength);
			for(size_t j = 0; j < interestRange.size(); j++)
			{
				if(Intersect(edgeRange, interestRange[j]))
				{
					interest.push_back(edge[i].GetStartVertex());
					interest.push_back(edge[i].GetEndVertex());
				}
			}
		}

		std::sort(interest.begin(), interest.end());
		interest.erase(std::unique(interest.begin(), interest.end()), interest.end());
		size_t threshold = (bifStorage.GetMaxId() * maxIterations) / PROGRESS_STRIDE;
		bool interestAffected = false;
		do
		{
			iterations++;			
			for(size_t id = 0; id <= bifStorage.GetMaxId(); id++)
			{
				totalBulges += RemoveBulges(sequence, bifStorage, k, minBranchSize, id, interest, interestAffected);
				if(++count >= threshold && !callBack.empty())
				{
					count = 0;
					totalProgress = std::min(totalProgress + 1, PROGRESS_STRIDE);
					callBack(totalProgress, run);
				}
			}

			if (interestAffected)
			{
				std::stringstream ss;
				ss << "investigate/dbg" << iterations << ".dot";
				std::ofstream out(ss.str().c_str());
				this->SerializeCondensedGraph(k, out, interest, sequence, bifStorage);
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
	}	
}