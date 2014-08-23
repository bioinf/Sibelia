//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "blockbuilder.h"
#include "indexedsequence.h"

namespace SyntenyFinder
{
	//namespace
	//{		
	//	template<class T>
	//		FancyIterator<T, char(*)(char), char> CompIterator(T it)
	//		{
	//			return CFancyIterator(it, FastaRecord::Translate, ' ');
	//		}
	//}

	//const size_t BlockBuilder::PROGRESS_STRIDE = 50;

	//BlockBuilder::BlockBuilder(const std::vector<FastaRecord> * originalChr, const std::string & tempDir):
	//	originalChr_(originalChr), tempDir_(tempDir)
	//{
	//}

	//void BlockBuilder::ConstructIndex(size_t k)
	//{
	//	if(index_.get() == 0)
	//	{
	//		lastK_ = k;
	//		size_t bifCount;
	//		std::vector<IndexedSequence::BifurcationInstance> bifurcation[2];
	//		{
	//			std::vector<std::string> record(originalChr_->size());
	//			for(size_t i = 0; i < record.size(); i++)
	//			{
	//				record[i].resize((*originalChr_)[i].GetSequence().size());
	//				virtualChrSize_.push_back(record[i].size());
	//				for(size_t j = 0; j < record[i].size(); j++)
	//				{
	//					char ch = (*originalChr_)[i].GetSequence()[j];
	//					bool change = !FastaRecord::IsDefiniteBase(ch) || FastaRecord::IsMaskedBase(ch);
	//					record[i][j] = change ? FastaRecord::DEFINITE_BASE[rand() % FastaRecord::DEFINITE_BASE.size()] : ch;
	//				}
	//			}
	//		
	//			bifCount = IndexedSequence::EnumerateBifurcationsSArray(record, k, tempDir_, bifurcation[0], bifurcation[1]);
	//		}

	//		index_.reset(new DeBruijnIndex(originalChr_->size(), bifCount));
	//		for(size_t strand = 0; strand < 2; strand++)
	//		{
	//			size_t nowBif = 0;
	//			FastaRecord::Direction dir = static_cast<FastaRecord::Direction>(strand);
	//			for(size_t chr = 0; chr < originalChr_->size(); chr++)
	//			{
	//				const std::string & seq = (*originalChr_)[chr].GetSequence();
	//				for(size_t pos = 0; pos < seq.size(); ++pos)
	//				{
	//					if(nowBif < bifurcation[strand].size() && chr == bifurcation[strand][nowBif].chr && pos == bifurcation[strand][nowBif].pos)
	//					{
	//						char mark = pos + k < seq.size() ? seq[pos + k] : IndexedSequence::SEPARATION_CHAR;
	//						index_->AddEdge(chr, pos, dir, bifurcation[strand][nowBif++].bifId, mark, pos);
	//					}
	//				}
	//			}
	//		}

	//	#ifdef _DEBUG
	//		for(size_t strand = 0; strand < 2; strand++)
	//		{
	//			FastaRecord::Direction dir = static_cast<FastaRecord::Direction>(strand);
	//			for(size_t i = 0; i < bifurcation[strand].size(); i++)
	//			{
	//				std::string kmer;
	//				IndexedSequence::BifurcationInstance bif = bifurcation[strand][i];
	//				const std::string & seq = (*originalChr_)[bif.chr].GetSequence();					
	//				if(dir == FastaRecord::positive)
	//				{
	//					kmer.assign(seq.begin() + bif.pos, seq.begin() + bif.pos + lastK_);
	//				}
	//				else
	//				{
	//					kmer.assign(CompIterator(seq.rbegin()), CompIterator(seq.rbegin()));
	//				}
	//				
	//				debugIndex_[kmer] = bif.bifId;
	//			}
	//		}

	//		virtualChr_.resize((*originalChr_).size());
	//		for(size_t chr = 0; chr < virtualChr_.size(); chr++)
	//		{
	//			virtualChr_[chr] = (*originalChr_)[chr].GetSequence();
	//		}
	//	#endif

	//	}
	//	else
	//	{
	//		//Reconstruct for new K
	//	}
	//}

	//void BlockBuilder::WriteIndexToDot(std::ostream & out) const
	//{
	//	out << "digraph G\n{\nrankdir=LR" << std::endl;
	//	for(size_t chr = 0; chr < originalChr_->size(); chr++)
	//	{
	//		for(size_t strand = 0; strand < 2; strand++)
	//		{
	//			FastaRecord::Direction dir = static_cast<FastaRecord::Direction>(strand);
	//			DeBruijnIndex::Edge prevEdge = index_->GetEdgeAtPosition(chr, 0, dir);
	//			for(size_t pos = 1; pos < virtualChrSize_[chr]; ++pos)
	//			{
	//				DeBruijnIndex::Edge nowEdge = index_->GetEdgeAtPosition(chr, pos, dir);
	//				if(nowEdge.Valid())
	//				{
	//					out << prevEdge.GetBifurcationId() << " -> " << nowEdge.GetBifurcationId() << "[";
	//					out << "color=" << (strand == 0 ? "blue" : "red") << " label=\"chr=" << chr;
	//					out << " pos=" << prevEdge.GetVirtualPosition() << ";" << pos;
	//					out << " proj=" << prevEdge.GetProjection() << ";" << nowEdge.GetProjection();
	//					out << " mark=" << prevEdge.GetMark() << "\"]" << std::endl;
	//					prevEdge = nowEdge;
	//				}
	//			}
	//		}
	//	}

	//	out << "}" << std::endl;
	//}

	//size_t BlockBuilder::Simplify(size_t minBranchSize, size_t maxIterations, ProgressCallBack callBack)
	//{				
	//	size_t totalProgress = 0;
	//	bool anyChanges = true;
	//	if(!callBack.empty())
	//	{
	//		callBack(totalProgress, start);
	//	}

	//	size_t count = 0;
	//	size_t totalBulges;
	//	size_t iterations = 0;
	//	size_t threshold = (index_->GetBifurcationsNumber() * maxIterations) / PROGRESS_STRIDE;
	//	do
	//	{
	//		iterations++;
	//		totalBulges = 0;
	//		for(size_t bifId = 0; bifId < index_->GetBifurcationsNumber(); bifId++)
	//		{			
	//			totalBulges += RemoveBulges(minBranchSize, bifId);
	//			if(++count >= threshold && !callBack.empty())
	//			{
	//				count = 0;
	//				totalProgress = std::min(totalProgress + 1, PROGRESS_STRIDE);
	//				callBack(totalProgress, run);
	//			}
	//		}
	//	}
	//	while((totalBulges > 0) && iterations < maxIterations);

	//	if(!callBack.empty())
	//	{
	//		callBack(PROGRESS_STRIDE, end);
	//	}
	//	
	//	return totalBulges;
	//}

	//namespace
	//{
	//	const char EMPTY = ' ';

	//	struct BranchData
	//	{
	//		BranchData() {}
	//		BranchData(char ch): endChar(ch) {}
	//		char endChar;
	//		std::vector<size_t> branchIds;
	//	};

	//	struct BifurcationMark
	//	{
	//		size_t bifId;
	//		size_t distance;
	//		BifurcationMark() {}
	//		BifurcationMark(size_t bifId, size_t distance): bifId(bifId),
	//			distance(distance) {}

	//		bool operator < (const BifurcationMark & compare) const
	//		{
	//			if(bifId != compare.bifId)
	//			{
	//				return bifId < compare.bifId;
	//			}

	//			return distance < compare.distance;
	//		}
	//	};

	//	void FillVisit(const DeBruijnIndex & index, DeBruijnIndex::Edge e, size_t minBranchSize, std::vector<BifurcationMark> & visit)
	//	{
	//		visit.clear();
	//		size_t start = e.GetBifurcationId();			
	//		for(size_t step = 1; step < minBranchSize; step++)
	//		{
	//			DeBruijnIndex::Edge nowEdge = index.GetEdgeAtPosition(e.GetChromosomeId(), e.GetPosition() + step, e.GetDirection());
	//			if(nowEdge.Valid())
	//			{
	//				size_t bifId = nowEdge.GetBifurcationId();
	//				if(bifId == start)
	//				{
	//					break;
	//				}

	//				visit.push_back(BifurcationMark(bifId, step));
	//			}
	//		}

	//		std::sort(visit.begin(), visit.end());
	//	}
	//	
	//	size_t MaxBifurcationMultiplicity(const DeBruijnIndex & index, DeBruijnIndex::Edge e, size_t distance)
	//	{
	//		size_t ret = 0;
	//		for(size_t step = 1; step < distance - 1; step++)
	//		{
	//			size_t pos = e.GetPosition() + step;
	//			DeBruijnIndex::Edge nowEdge = index.GetEdgeAtPosition(e.GetChromosomeId(), pos, e.GetDirection());
	//			if(nowEdge.Valid())
	//			{
	//				ret = std::max(ret, index.CountEdges(nowEdge.GetBifurcationId()));
	//			}
	//		}

	//		return ret;
	//	}

	//	typedef std::vector< std::vector<size_t> > BulgedBranches;
	//	
	//	bool AnyBulges(const DeBruijnIndex & index,
	//		std::vector<DeBruijnIndex::Edge> & edge,			
	//		BulgedBranches & bulges,
	//		size_t minBranchSize)
	//	{
	//		
	//		bulges.clear();
	//		boost::unordered_map<size_t, BranchData> visit;
	//		for(size_t i = 0; i < edge.size(); i++)
	//		{
	//			size_t start = edge[i].GetBifurcationId();
	//			for(size_t step = 1; step < minBranchSize; step++)
	//			{
	//				size_t pos = edge[i].GetPosition() + step;
	//				DeBruijnIndex::Edge nowEdge = index.GetEdgeAtPosition(edge[i].GetChromosomeId(), pos, edge[i].GetDirection()); 
	//				if(!nowEdge.Valid())
	//				{
	//					continue;
	//				}

	//				size_t bifId = nowEdge.GetBifurcationId();
	//				if(bifId == start)
	//				{
	//					break;
	//				}

	//				boost::unordered_map<size_t, BranchData>::iterator kt = visit.find(bifId);
	//				if(kt == visit.end())
	//				{
	//					BranchData bData(edge[i].GetMark());
	//					bData.branchIds.push_back(i);
	//					visit[bifId] = bData;
	//				}
	//				else if(kt->second.endChar != edge[i].GetMark())
	//				{
	//					kt->second.branchIds.push_back(i);
	//					break;
	//				}
	//			}
	//		}

	//		for (boost::unordered_map<size_t, BranchData>::iterator kt = visit.begin(); kt !=visit.end(); ++kt)
	//		{
	//			if (kt->second.branchIds.size() > 1)
	//			{
	//				bulges.push_back(kt->second.branchIds);		
	//			}
	//		}
	//		
	//		return !bulges.empty();
	//	}			
	//}

	////OPTIMIZE THIS!!!1111
	//bool BlockBuilder::Overlap(const std::vector<DeBruijnIndex::Edge> & edge, VisitData sourceData, VisitData targetData) const
	//{
	//	DeBruijnIndex::Edge sourceEdge = edge[sourceData.kmerId];
	//	DeBruijnIndex::Edge targetEdge = edge[targetData.kmerId];
	//	if(sourceEdge.GetChromosomeId() != targetEdge.GetChromosomeId())
	//	{
	//		return false;
	//	}

	//	std::vector<size_t> occur;			
	//	for(size_t i = 0; i < sourceData.distance + lastK_; i++)
	//	{
	//		size_t pos = sourceEdge.GetPosition() + i;
	//		occur.push_back(sourceEdge.GetPosition());
	//	}
	//		
	//	std::sort(occur.begin(), occur.end());
	//	for(size_t i = 0; i < targetData.distance + lastK_; i++)
	//	{
	//		size_t pos = targetEdge.GetPosition() + i;
	//		if(std::binary_search(occur.begin(), occur.end(), pos))
	//		{
	//			return true;
	//		}
	//	}
	//		
	//	return false;
	//}

	//#ifdef _DEBUG
	//	void BlockBuilder::PrintRaw(size_t chr0, size_t chr1, std::ostream & out) const
	//	{
	//		size_t chrId[] = {chr0, chr1};
	//		for(size_t i = 0; i < 2; i++)
	//		{
	//			size_t chr = chrId[i];
	//			out << "Sequence #" << chr << std::endl;
	//			for(size_t i = 0; i < virtualChr_[chr].size(); i++)
	//			{
	//				out << i % 10;
	//			}

	//			out << std::endl;
	//			std::copy(virtualChr_[chr].begin(), virtualChr_[chr].end(), std::ostream_iterator<char>(out));
	//			out << std::endl;				
	//			std::copy(CompIterator(virtualChr_[chr].begin()), CompIterator(virtualChr_[chr].end()), std::ostream_iterator<char>(out));
	//			out << std::endl;
	//		}
	//	}
	//	
	//	void BlockBuilder::PrintPath(DeBruijnIndex::Edge e, size_t distance, std::ostream & out) const
	//	{
	//		out << (e.GetDirection() == FastaRecord::positive ? "+" : "-") << e.GetPosition() << ' ';
	//		if(e.GetDirection() == FastaRecord::positive)
	//		{
	//			CopyN(virtualChr_[e.GetChromosomeId()].begin() + e.GetPosition(), distance + lastK_, std::ostream_iterator<char>(out));
	//		}
	//		else
	//		{				
	//			CopyN(CompIterator(virtualChr_[e.GetChromosomeId()].rbegin() + e.GetPosition()), distance + lastK_, std::ostream_iterator<char>(out));
	//		}

	//		std::cerr << std::endl;
	//	}
	//#endif
	//	
	//	void BlockBuilder::CollapseBulgeGreedily(std::vector<DeBruijnIndex::Edge> & edge, VisitData sourceData, VisitData targetData)
	//	{
	//		DeBruijnIndex::Edge sourceEdge = edge[sourceData.kmerId];
	//		DeBruijnIndex::Edge targetEdge = edge[targetData.kmerId];
	//	#ifdef _DEBUG

	//		static size_t bulge = 0;
	//		std::cerr << "Bulge #" << bulge++ << std::endl;
	//		std::cerr << "Before: " << std::endl;
	//		PrintRaw(sourceEdge.GetChromosomeId(), targetEdge.GetChromosomeId(), std::cerr);
	//		std::cerr << "Source branch: " << std::endl;
	//		PrintPath(sourceEdge, sourceData.distance, std::cerr);
	//		std::cerr << "Target branch: " << std::endl;
	//		PrintPath(targetEdge, targetData.distance, std::cerr);
	//		//bifStorage.Dump(sequence, k, std::cerr);
	//		//iseq_->Test();
	//	#endif
	//	/*	std::vector<std::pair<size_t, size_t> > lookForward;
	//		std::vector<std::pair<size_t, size_t> > lookBack;
	//		EraseBifurcations(sequence, bifStorage, k, startKMer, targetData, lookForward, lookBack);
	//		StrandIterator sourceIt = *startKMer[sourceData.kmerId];
	//		StrandIterator targetIt = *startKMer[targetData.kmerId];
	//		sequence.Replace(AdvanceForward(sourceIt, k),
	//			sourceData.distance,
	//			AdvanceForward(targetIt, k),
	//			targetData.distance,
	//			boost::bind(&BifurcationStorage::NotifyBefore, boost::ref(bifStorage), _1, _2),
	//			boost::bind(&BifurcationStorage::NotifyAfter, boost::ref(bifStorage), _1, _2));
	//		UpdateBifurcations(sequence, bifStorage, k, startKMer, sourceData, targetData, lookForward, lookBack);
	//		*/
	//	#ifdef _DEBUG/*
	//		std::cerr << "After: " << std::endl;
	//		BlockFinder::PrintRaw(sequence, std::cerr);
	//		std::cerr << "Source branch: " << std::endl;s
	//		BlockFinder::PrintPath(sequence, *startKMer[sourceData.kmerId], k, sourceData.distance, std::cerr);
	//		std::cerr << "Target branch: " << std::endl;
	//		BlockFinder::PrintPath(sequence, *startKMer[targetData.kmerId], k, sourceData.distance, std::cerr);
	//		bifStorage.Dump(sequence, k, std::cerr);
	//		iseq_->Test();
	//		std::cerr << DELIMITER << std::endl;*/
	//	#endif
	//	}

	//size_t BlockBuilder::RemoveBulges(size_t minBranchSize, size_t bifId)
	//{
	//	size_t ret = 0;
	//	std::vector<DeBruijnIndex::Edge> edge;
	//	if(index_->GetEdgesOfVertex(bifId, edge) < 2)
	//	{
	//		return ret;
	//	}
	//	
	//	BulgedBranches bulges;
	//	if(!AnyBulges(*index_, edge, bulges, minBranchSize))
	//	{
	//		return ret;
	//	}

	//	std::vector<BifurcationMark> visit;
	//	for (size_t numBulge = 0; numBulge < bulges.size(); ++numBulge)
	//	{
	//		for (size_t  idI = 0; idI < bulges[numBulge].size(); ++idI)
	//		{
	//			size_t kmerI = bulges[numBulge][idI];
	//			FillVisit(*index_, edge[kmerI], minBranchSize, visit);
	//			for(size_t  idJ = idI + 1; idJ < bulges[numBulge].size(); ++idJ)
	//			{
	//				size_t kmerJ = bulges[numBulge][idJ];
	//				if(edge[kmerI].GetMark() == edge[kmerJ].GetMark())
	//				{
	//					continue;
	//				}
	//				
	//				for(size_t step = 1; step < minBranchSize; step++)
	//				{
	//					size_t pos = edge[kmerJ].GetPosition() + step;
	//					DeBruijnIndex::Edge nowEdge = index_->GetEdgeAtPosition(edge[kmerJ].GetChromosomeId(), pos, edge[kmerJ].GetDirection());
	//					if(nowEdge.Valid())
	//					{
	//						size_t nowBif = nowEdge.GetBifurcationId();
	//						if(nowBif == bifId)
	//						{
	//							break;
	//						}

	//						std::vector<BifurcationMark>::iterator vt = std::lower_bound(visit.begin(), visit.end(), BifurcationMark(nowBif, 0));
	//						if(vt != visit.end() && vt->bifId == nowBif)
	//						{
	//							VisitData jdata(kmerJ, step);
	//							VisitData idata(kmerI, vt->distance);
	//							if(Overlap(edge, idata, jdata) || nowBif == bifId)
	//							{
	//								break;
	//							}

	//							++ret;
	//							size_t imlp = MaxBifurcationMultiplicity(*index_, edge[kmerI], idata.distance);
	//							size_t jmlp = MaxBifurcationMultiplicity(*index_, edge[kmerI], idata.distance);
	//							bool iless = imlp > jmlp || (imlp == jmlp && idata.kmerId < jdata.kmerId);
	//							if(iless)
	//							{
	//								DeBruijnIndex::Edge & replace = edge[jdata.kmerId];
	//								CollapseBulgeGreedily(edge, idata, jdata);
	//								replace = index_->GetEdgeAtPosition(replace.GetChromosomeId(), replace.GetPosition(), replace.GetDirection());
	//							}
	//							else
	//							{
	//								DeBruijnIndex::Edge & replace = edge[idata.kmerId];
	//								CollapseBulgeGreedily(edge, jdata, idata);										
	//								replace = index_->GetEdgeAtPosition(replace.GetChromosomeId(), replace.GetPosition(), replace.GetDirection());
	//								FillVisit(*index_, edge[kmerI], minBranchSize, visit);
	//							}

	//							break;
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}

	//	return ret;
	//}

	//void BlockBuilder::GenerateBlocks(std::vector<BlockInstance> & ret, size_t minBlockSize) const
	//{
	//}
}