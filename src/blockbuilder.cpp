//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "blockbuilder.h"

namespace SyntenyFinder
{
	namespace
	{		
		template<class T>
			FancyIterator<T, char(*)(char), char> CompIterator(T it)
			{
				return CFancyIterator(it, FastaRecord::Translate, ' ');
			}
	}

	const size_t BlockBuilder::PROGRESS_STRIDE = 50;

	BlockBuilder::BlockBuilder(const std::vector<FastaRecord> * originalChr, const std::string & tempDir):
		originalChr_(originalChr), tempDir_(tempDir)
	{
	}

	void BlockBuilder::ConstructIndex(size_t k)
	{
		if(index_.get() == 0)
		{
			lastK_ = k;
			size_t bifCount;
			std::vector<size_t> originalChrSize;
			std::vector<DeBruijnIndex::ChrBifVector> bifurcation;
			std::vector<std::string> record(originalChr_->size());
			for(size_t i = 0; i < record.size(); i++)
			{
				originalChrSize.push_back(record[i].size());
				record[i].resize((*originalChr_)[i].GetSequence().size());
				virtualChrSize_.push_back(record[i].size());
				for(size_t j = 0; j < record[i].size(); j++)
				{
					char ch = (*originalChr_)[i].GetSequence()[j];
					bool change = !FastaRecord::IsDefiniteBase(ch) || FastaRecord::IsMaskedBase(ch);
					record[i][j] = change ? FastaRecord::DEFINITE_BASE[rand() % FastaRecord::DEFINITE_BASE.size()] : ch;
				}
			}
			
			bifCount = EnumerateBifurcationsSArray(record, k, tempDir_, bifurcation);
			index_.reset(new DeBruijnIndex(bifurcation, record, lastK_, bifCount, originalChrSize));
		}
		else
		{
			//Reconstruct for new K
		}
	}

	void BlockBuilder::WriteIndexToDot(std::ostream & out) const
	{
		out << "digraph G\n{\nrankdir=LR" << std::endl;
		for(size_t chr = 0; chr < originalChr_->size(); chr++)
		{
			for(size_t strand = 0; strand < 2; strand++)
			{
				FastaRecord::Direction dir = static_cast<FastaRecord::Direction>(strand);
				DeBruijnIndex::BifurcationIterator jt = index_->Begin(chr, dir);
				DeBruijnIndex::BifurcationIterator end = index_->End(chr, dir);
				if(jt != end)
				{
					DeBruijnIndex::BifurcationIterator it = jt++;
					for(; jt != end; ++it, ++jt)
					{
						out << it.GetBifurcationId() << " -> " << jt.GetBifurcationId() << "[";
						out << "color=" << (strand == 0 ? "blue" : "red") << " label=\"chr=" << chr;
						out << " pos=" << it.GetPosition() << ";" << jt.GetPosition();
						out << " proj=" << it.GetProjection() << ";" << jt.GetProjection();
						out << " mark=" << it.GetOutMark() << ";" << jt.GetInMark() << "\"]" << std::endl;
					}
				}
			}
		}

		out << "}" << std::endl;
	}

	size_t BlockBuilder::Simplify(size_t maxBranchSize, size_t maxIterations, ProgressCallBack callBack)
	{				
		size_t totalProgress = 0;
		bool anyChanges = true;
		if(!callBack.empty())
		{
			callBack(totalProgress, start);
		}

		size_t count = 0;
		size_t totalBulges;
		size_t iterations = 0;
		size_t threshold = (index_->GetBifurcationsNumber() * maxIterations) / PROGRESS_STRIDE;
		do
		{
			iterations++;
			totalBulges = 0;
			for(size_t bifId = 0; bifId < index_->GetBifurcationsNumber(); bifId++)
			{			
				totalBulges += RemoveBulges(maxBranchSize, bifId);
				if(++count >= threshold && !callBack.empty())
				{
					count = 0;
					totalProgress = std::min(totalProgress + 1, PROGRESS_STRIDE);
					callBack(totalProgress, run);
				}
			}

			index_->ApplyChanges();
		}
		while((totalBulges > 0) && iterations < maxIterations);

		if(!callBack.empty())
		{
			callBack(PROGRESS_STRIDE, end);
		}
		
		return totalBulges;
	}

	namespace
	{
		const size_t INVALID_BRANCH = -1;

		void PrintBranch(DeBruijnIndex::BifurcationIterator it, size_t distance, std::ostream & out)
		{
			out << (it.GetStrand() == FastaRecord::positive ? "+" : "-") << it.GetChromosomeId() << " ";
			for(size_t i = 0; i < distance + 1; i++, ++it)
			{
				out << "{Id=" << it.GetBifurcationId() << ", Pos=" << it.GetPosition() << "} ";
			}

			out << std::endl;
		}
	
		size_t MaxBifurcationMultiplicity(const DeBruijnIndex & index, DeBruijnIndex::BifurcationIterator it, size_t distance, std::vector<size_t> & interim)
		{
			size_t ret = 0;
			for(size_t step = 0; step < distance; step++, ++it)
			{
				if(step != 0 && step != distance - 1)
				{
					size_t pos = it.GetPosition() + step;
					interim.push_back(it.GetBifurcationId());
					ret = std::max(ret, index.CountInstances(it.GetBifurcationId()));
				}

				if(!it.IsValid())
				{
					return INVALID_BRANCH;
				}
			}

			std::sort(interim.begin(), interim.end());
			return ret;
		}
	}

	bool BlockBuilder::AnyBulges(std::vector<DeBruijnIndex::BifurcationIterator> bif, BulgedBranches & bulges, size_t maxBranchSize) const
	{
		bulges.clear();
		std::vector<char> startOutMark(bif.size());
		for(size_t i = 0; i < bif.size(); i++)
		{
			startOutMark[i] = bif[i].GetOutMark();
		}

		boost::unordered_map<size_t, BranchData> visit;
		for(size_t i = 0; i < bif.size(); i++)
		{	
			size_t maxBifMlp = 0;
			size_t startId = bif[i].GetBifurcationId();
			size_t startPos = (bif[i]++).GetPosition();
			for(size_t step = 1; !bif[i].AtEnd(); step++, ++bif[i])
			{					
				size_t pos = bif[i].GetPosition();
				size_t bifId = bif[i].GetBifurcationId();
				if(!bif[i].IsValid() || pos - startPos > maxBranchSize || bifId == startId)
				{
					break;
				}
					
				boost::unordered_map<size_t, BranchData>::iterator kt = visit.find(bifId);
				if(kt == visit.end())
				{
					BranchData bData(startOutMark[i], maxBifMlp);
					bData.branch.push_back(VisitData(i, step));
					visit[bifId] = bData;
				}
				else if(kt->second.endChar != startOutMark[i])
				{						
					kt->second.branch.push_back(VisitData(i, step));
					break;
				}
				else if(maxBifMlp > kt->second.maxBifMlp)
				{
					kt->second.maxBifMlp = maxBifMlp;
					kt->second.branch[0] = VisitData(i, step);
				}

				maxBifMlp = std::max(maxBifMlp, index_->CountInstances(bifId));
			}
		}

		for (boost::unordered_map<size_t, BranchData>::iterator kt = visit.begin(); kt !=visit.end(); ++kt)
		{
			if (kt->second.branch.size() > 1)
			{
				bulges.push_back(kt->second.branch);
			}
		}
			
		return !bulges.empty();
	}		

	bool BlockBuilder::Overlap(const std::vector<DeBruijnIndex::BifurcationIterator> & bif, VisitData sourceData, VisitData targetData) const
	{
		DeBruijnIndex::BifurcationIterator srcIt = bif[sourceData.kmerId];
		DeBruijnIndex::BifurcationIterator trgIt = bif[targetData.kmerId];
		if(srcIt.GetChromosomeId() != trgIt.GetChromosomeId())
		{
			return false;
		}

		std::pair<size_t, size_t> srcRange(srcIt.GetPositivePosition(), (srcIt + sourceData.distance).GetPositiveEndingPosition());
		std::pair<size_t, size_t> trgRange(trgIt.GetPositivePosition(), (trgIt + sourceData.distance).GetPositiveEndingPosition());	
		return RangesOverlap(srcRange, trgRange);
	}

	void BlockBuilder::CollapseBulge(const std::vector<DeBruijnIndex::BifurcationIterator> & bif, VisitData sourceData, VisitData targetData)
	{
		DeBruijnIndex::BifurcationIterator srcIt = bif[sourceData.kmerId];
		DeBruijnIndex::BifurcationIterator trgIt = bif[targetData.kmerId];
	#ifdef _DEBUG
		static size_t bulge = 0;
		std::cerr << "Bulge #" << bulge++ << std::endl;
		std::cerr << "Before: " << std::endl;
		std::cerr << "Source branch: " << std::endl;
		PrintBranch(srcIt, sourceData.distance, std::cerr);
		std::cerr << "Target branch: " << std::endl;
		PrintBranch(trgIt, targetData.distance, std::cerr);
	#endif
		index_->Replace(srcIt, srcIt + sourceData.distance, trgIt, trgIt + targetData.distance);
	}

	size_t BlockBuilder::RemoveBulges(size_t maxBranchSize, size_t bifId)
	{
		size_t ret = 0;
		std::vector<DeBruijnIndex::BifurcationIterator> bif;
		if(index_->GetBifurcationInstances(bifId, bif) < 2)
		{
			return ret;
		}
		
		BulgedBranches bulges;
		if(!AnyBulges(bif, bulges, maxBranchSize))
		{
			return ret;
		}

		for (size_t numBulge = 0; numBulge < bulges.size(); ++numBulge)
		{
			for (size_t idI = 0; idI < bulges[numBulge].size(); ++idI)
			{
				VisitData idata = bulges[numBulge][idI];
				for(size_t idJ = idI + 1; idJ < bulges[numBulge].size() && bif[idata.kmerId].IsValid(); ++idJ)
				{
					VisitData jdata = bulges[numBulge][idJ];
					if(Overlap(bif, idata, jdata) || !bif[jdata.kmerId].IsValid())
					{
						continue;
					}

					std::vector<size_t> interim[2];
					size_t imlp = MaxBifurcationMultiplicity(*index_, bif[idata.kmerId], idata.distance, interim[0]);
					size_t jmlp = MaxBifurcationMultiplicity(*index_, bif[jdata.kmerId], jdata.distance, interim[0]);					
					if(imlp != INVALID_BRANCH && jmlp != INVALID_BRANCH)
					{
						bool foundInterim = false;
						for(size_t idx = 0; idx < interim[0].size() && !foundInterim; idx++)
						{
							foundInterim = std::binary_search(interim[1].begin(), interim[1].end(), interim[0][idx]);
						}

						if(!foundInterim)
						{
							++ret;
							bool iless = imlp > jmlp || (imlp == jmlp && idata.kmerId < jdata.kmerId);
							if(iless)
							{
								CollapseBulge(bif, idata, jdata);
							}
							else
							{
								CollapseBulge(bif, jdata, idata);										
							}

							break;
						}
					}
				}
			}
		}

		return ret;
	}

	void BlockBuilder::GenerateBlocks(std::vector<BlockInstance> & ret, size_t minBlockSize) const
	{
	}
}