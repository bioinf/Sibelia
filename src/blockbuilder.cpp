//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "blockbuilder.h"
#include "indexedsequence.h"

namespace SyntenyFinder
{
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
			std::vector<std::string> record(originalChr_->size());
			for(size_t i = 0; i < record.size(); i++)
			{
				record[i].resize((*originalChr_)[i].GetSequence().size());
				virtualChrSize_.push_back(record[i].size());
				for(size_t j = 0; j < record[i].size(); j++)
				{
					char ch = (*originalChr_)[i].GetSequence()[j];
					bool change = !FastaRecord::IsDefiniteBase(ch) || FastaRecord::IsMaskedBase(ch);
					record[i][j] = change ? FastaRecord::DEFINITE_BASE[rand() % FastaRecord::DEFINITE_BASE.size()] : ch;
				}
			}

			std::vector<IndexedSequence::BifurcationInstance> bifurcation[2];
			size_t bifCount = IndexedSequence::EnumerateBifurcationsSArray(record, k, tempDir_, bifurcation[0], bifurcation[1]);
			index_.reset(new DeBruijnIndex(originalChr_->size(), bifCount));
			for(size_t strand = 0; strand < 2; strand++)
			{
				size_t nowBif = 0;
				FastaRecord::Direction dir = static_cast<FastaRecord::Direction>(strand);
				for(size_t chr = 0; chr < originalChr_->size(); chr++)
				{
					for(size_t pos = 0; pos < (*originalChr_)[chr].GetSequence().size(); ++pos)
					{
						if(nowBif < bifurcation[strand].size() && chr == bifurcation[strand][nowBif].chr && pos == bifurcation[strand][nowBif].pos)
						{
							char mark = pos + k < record[chr].size() ? record[chr][pos + k] : IndexedSequence::SEPARATION_CHAR;
							index_->AddEdge(chr, pos, dir, bifurcation[strand][nowBif++].bifId, mark, pos);
						}
					}
				}
			}

		#ifdef _DEBUG
			for(size_t strand = 0; strand < 2; strand++)
			{
				size_t nowBif = 0;
				FastaRecord::Direction dir = static_cast<FastaRecord::Direction>(strand);
				std::vector<std::map<size_t, size_t> > posBif(originalChr_->size());				
				for(size_t chr = 0; chr < virtualChrSize_.size(); chr++)
				{
					for(size_t i = 0; i < bifurcation[strand].size(); i++)
					{
						posBif[chr][bifurcation[strand][i].pos] = bifurcation[strand][i].bifId;
					}
				}

				for(size_t chr = 0; chr < virtualChrSize_.size(); chr++)
				{
					for(size_t pos = 0; pos < virtualChrSize_[chr]; pos++)
					{
						DeBruijnIndex::Edge e = index_->GetEdgeAtPosition(chr, pos, dir);
						std::map<size_t, size_t>::const_iterator it = posBif[chr].find(pos);
						assert((e.Valid() && e.GetBifurcationId() == posBif[chr][pos]) || (!e.Valid() && it == posBif[chr].end()));
					}
				}
			}
		#endif

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
				DeBruijnIndex::Edge prevEdge = index_->GetEdgeAtPosition(chr, 0, dir);
				for(size_t pos = 1; pos < virtualChrSize_[chr]; ++pos)
				{
					DeBruijnIndex::Edge nowEdge = index_->GetEdgeAtPosition(chr, pos, dir);
					if(nowEdge.Valid())
					{
						out << prevEdge.GetBifurcationId() << " -> " << nowEdge.GetBifurcationId() << "[";
						out << "color=" << (strand == 0 ? "blue" : "red") << " label=\"chr=" << chr;
						out << " pos=" << prevEdge.GetVirtualPosition() << ";" << pos;
						out << " proj=" << prevEdge.GetProjection() << ";" << nowEdge.GetProjection();
						out << " mark=" << prevEdge.GetMark() << "\"]" << std::endl;
						prevEdge = nowEdge;
					}
				}
			}
		}

		out << "}" << std::endl;
	}

	size_t BlockBuilder::Simplify(size_t minBranchSize, size_t maxIterations, ProgressCallBack callBack)
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
				totalBulges += RemoveBulges(minBranchSize, bifId);
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

	size_t BlockBuilder::RemoveBulges(size_t minBranchSize, size_t bifId)
	{
		size_t ret = 0;
	/*	IteratorProxyVector startKMer;
		if(bifStorage.ListPositions(bifId, std::back_inserter(startKMer)) < 2)
		{
			return ret;
		}

		std::vector<char> endChar(startKMer.size(), EMPTY);
		for(size_t i = 0; i < startKMer.size(); i++)
		{
			if(ProperKMer(*startKMer[i], k + 1))
			{
				endChar[i] = *AdvanceForward(*startKMer[i], k);
			}
		}

		//std::vector<bool> isBulge(startKMer.size(), false);
		BulgedBranches bulges;
		if(!AnyBulges(sequence, bifStorage, k, startKMer, endChar, bulges, minBranchSize))
		{
			return ret;
		}

		std::vector<BifurcationMark> visit;
		for (size_t numBulge = 0; numBulge < bulges.size(); ++numBulge)
		{

			//for(size_t kmerI = 0; kmerI < startKMer.size(); kmerI++)
			for (size_t  idI = 0; idI < bulges[numBulge].size(); ++idI)
			{
				size_t kmerI = bulges[numBulge][idI];
				//if(!startKMer[kmerI].Valid() || !isBulge[kmerI])
				if(!startKMer[kmerI].Valid())
				{
					//if (!isBulge[kmerI]) std::cerr << "optimised!\n";
					continue;
				}

				FillVisit(sequence, bifStorage, *startKMer[kmerI], minBranchSize, visit);
				//for(size_t kmerJ = kmerI + 1; kmerJ < startKMer.size(); kmerJ++)
				for(size_t  idJ = idI + 1; idJ < bulges[numBulge].size(); ++idJ)
				{
					size_t kmerJ = bulges[numBulge][idJ];
					//if(!startKMer[kmerJ].Valid() || endChar[kmerI] == endChar[kmerJ] || !isBulge[kmerJ])
					if(!startKMer[kmerJ].Valid() || endChar[kmerI] == endChar[kmerJ])
					{
						//if (!isBulge[kmerJ]) std::cerr << "optimised!\n";
						continue;
					}

					StrandIterator kmer = ++StrandIterator(*startKMer[kmerJ]);
					for(size_t step = 1; kmer.AtValidPosition() && step < minBranchSize; ++kmer, step++)
					{
						size_t nowBif = bifStorage.GetBifurcation(kmer);
						if(nowBif != BifurcationStorage::NO_BIFURCATION)
						{
							if(nowBif == bifId)
							{
								break;
							}

							std::vector<BifurcationMark>::iterator vt = std::lower_bound(visit.begin(), visit.end(), BifurcationMark(nowBif, 0));
							if(vt != visit.end() && vt->bifId == nowBif)
							{
								VisitData jdata(kmerJ, step);
								VisitData idata(kmerI, vt->distance);
								if(Overlap(k, startKMer, idata, jdata) || nowBif == bifId)
								{
									break;
								}

								++ret;
								size_t imlp = MaxBifurcationMultiplicity(bifStorage, *startKMer[kmerI], idata.distance);
								size_t jmlp = MaxBifurcationMultiplicity(bifStorage, *startKMer[kmerJ], jdata.distance);
								bool iless = imlp > jmlp || (imlp == jmlp && idata.kmerId < jdata.kmerId);
								if(iless)
								{
									endChar[jdata.kmerId] = endChar[idata.kmerId];
									CollapseBulgeGreedily(sequence, bifStorage, k, startKMer, idata, jdata);
								}
								else
								{
									endChar[idata.kmerId] = endChar[jdata.kmerId];
									CollapseBulgeGreedily(sequence, bifStorage, k, startKMer, jdata, idata);
									FillVisit(sequence, bifStorage, *startKMer[kmerI], minBranchSize, visit);
								}

								break;
							}
						}
					}
				}
			}
		}*/

		return ret;
	}

	void BlockBuilder::GenerateBlocks(std::vector<BlockInstance> & ret, size_t minBlockSize) const
	{
	}
}