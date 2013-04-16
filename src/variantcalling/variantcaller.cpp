//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "variantcaller.h"
#include <seqan/align.h>
#undef min
#undef max

namespace SyntenyFinder
{
	typedef std::vector<BlockInstance>::const_iterator BLCIterator;

	namespace 
	{
		const size_t oo = UINT_MAX;
		bool NextVertex(IndexedSequence & iseq, StrandIterator & reference, StrandIterator referenceEnd, StrandIterator & assembly, StrandIterator assemblyEnd, bool strict)
		{
			size_t minSum = oo;						
			IteratorProxyVector startKMer;
			StrandIterator nowReference = reference;
			StrandIterator nextReference = referenceEnd;
			StrandIterator nextAssembly = assemblyEnd;	
			size_t assemblyChr = iseq.GetChr(assembly);
			for(size_t step = 0; step < minSum && nowReference.AtValidPosition(); ++nowReference, ++step)
			{
				size_t bifId = iseq.BifStorage().GetBifurcation(nowReference);
				if(bifId != BifurcationStorage::NO_BIFURCATION)
				{
					startKMer.clear();
					iseq.BifStorage().ListPositions(bifId, std::back_inserter(startKMer));
					for(size_t i = 0; i < startKMer.size(); i++)
					{
						StrandIterator kmer = *startKMer[i];
						if(iseq.GetChr(kmer) == assemblyChr && IndexedSequence::StrandIteratorPosGEqual(kmer, assembly))
						{
							size_t sum = step + IndexedSequence::StrandIteratorDistance(assembly, kmer);
							if((sum > 0 || !strict) && sum < minSum)
							{
								nextReference = nowReference;
								nextAssembly = kmer;
								minSum = sum;
							}
						}
					}
				}				
			}

			reference = nextReference;
			assembly = nextAssembly;
			return minSum != oo;
		}

		size_t Traverse(StrandIterator reference, StrandIterator assembly, size_t bound = UINT_MAX)
		{
			size_t ret = 0;
			for(;ret < bound && reference.AtValidPosition() && assembly.AtValidPosition() && *reference == *assembly; ++ret)
			{
				++reference;
				++assembly;
			}

			return ret;
		}

		typedef std::vector<BlockInstance>::iterator BLIterator;
		bool SignedBlocksIdEquals(const BlockInstance & a, const BlockInstance & b)
		{
			return a.GetSignedBlockId() == b.GetSignedBlockId();
		}

		void Reverse(BLCIterator begin, BLCIterator end, BLIterator out)
		{
			std::copy(begin, end, out);
			std::for_each(out, out + (end - begin), boost::bind(&BlockInstance::Reverse, _1));
			std::reverse(out, out + (end - begin));
		}
		
		BLCIterator Apply(BLCIterator begin, BLCIterator end, BLCIterator patternBegin, BLCIterator patternEnd)
		{
			size_t size = patternEnd - patternBegin;
			for(BLCIterator it = begin; it < end - size + 1; ++it)
			{
				if(std::equal(patternBegin, patternEnd, it, SignedBlocksIdEquals))
				{
					return it;
				}
			}

			return end;
		}
	}

	bool VariantCaller::SearchInReference(const std::string & stringPattern) const
	{		
		const std::string & sequence = (*chr_)[refSeqId_].GetSequence();
		std::vector<sauchar_t> pattern(stringPattern.begin(), stringPattern.end());
		saidx_t count = sa_search(reinterpret_cast<const sauchar_t*>(sequence.c_str()), sequence.size(), &pattern[0], pattern.size(), &referenceSuffixArray_[0], referenceSuffixArray_.size(), &indexOut_[0]);
		return count == 0;
	}

	
	bool VariantCaller::SearchInAssembly(const std::string & stringPattern) const
	{
		std::vector<sauchar_t> pattern(stringPattern.begin(), stringPattern.end());
		saidx_t count = sa_search(reinterpret_cast<const sauchar_t*>(assemblySequence_.c_str()), assemblySequence_.size(), &pattern[0], pattern.size(), &assemblySuffixArray_[0], assemblySuffixArray_.size(), &indexOut_[0]);
		return count == 0;
	}

	bool VariantCaller::ConfirmVariant(StrandIterator referenceStart, StrandIterator referenceEnd, StrandIterator assemblyStart, StrandIterator assemblyEnd) const
	{
		size_t confirmDist = trimK_ * 2;
		for(size_t i = 0; i < confirmDist; i++)
		{
			StrandIterator nextReferenceStart = AdvanceBackward(referenceStart, 1);
			StrandIterator nextAssemblyStart = AdvanceBackward(assemblyStart, 1);
			if(nextReferenceStart.AtValidPosition() && assemblyStart.AtValidPosition() && *referenceStart == *assemblyStart)
			{
				referenceStart = nextReferenceStart;
				assemblyStart = nextAssemblyStart;
			}
			else
			{
				break;
			}
		}
		
		size_t forward = Traverse(referenceEnd, assemblyEnd, confirmDist);
		std::advance(referenceEnd, forward);
		std::advance(assemblyEnd, forward);		
		std::string pattern(assemblyStart, assemblyEnd);
		return SearchInReference(pattern);
	}

	VariantCaller::VariantCaller(const std::vector<FASTARecord> & chr, size_t refSeqId, const std::vector<std::vector<BlockInstance> > & history, size_t trimK, size_t minBlockSize):
		chr_(&chr), refSeqId_(refSeqId), history_(history), trimK_(trimK), minBlockSize_(minBlockSize)
	{		
		const std::string & reference = (*chr_)[refSeqId_].GetSequence();							
		std::stringstream ss;
		for(size_t i = 0; i < chr_->size(); i++)
		{
			ss << (*chr_)[i].GetSequence() << "$";
		}
		
		assemblySequence_ = ss.str();
		assemblySuffixArray_.resize(assemblySequence_.size());	
		referenceSuffixArray_.resize(reference.size());	
		divsufsort(reinterpret_cast<const sauchar_t*>(assemblySequence_.c_str()), &assemblySuffixArray_[0], static_cast<saidx_t>(assemblySequence_.size()));
		divsufsort(reinterpret_cast<const sauchar_t*>(reference.c_str()), &referenceSuffixArray_[0], static_cast<saidx_t>(reference.size()));
		indexOut_.resize(reference.size() + assemblySequence_.size());
	}

	void VariantCaller::AlignBulgeBranches(size_t blockId, StrandIterator referenceBegin, StrandIterator referenceEnd, StrandIterator assemblyBegin, StrandIterator assemblyEnd, const FASTARecord & assemblySequence, std::vector<Variant> & variantList) const
	{
		using namespace seqan;
		typedef String<char> TSequence;
		typedef Align<TSequence, ArrayGaps>	TAlign;
		typedef Row<TAlign>::Type TRow;
		typedef Iterator<TRow>::Type TIterator;
		typedef Position<TAlign>::Type TPosition;
		
		size_t contextEnd = Traverse(referenceEnd, assemblyEnd, trimK_);
		std::advance(referenceEnd, contextEnd);
		std::advance(assemblyEnd, contextEnd);
		std::string referenceContext(referenceBegin, referenceEnd);
		std::string assemblyContext(assemblyBegin, assemblyEnd);
		std::string referenceBuf(referenceBegin, referenceEnd);
		std::string assemblyBuf(assemblyBegin, assemblyEnd);
		TSequence referenceSeq(referenceBuf);
		TSequence assemblySeq(assemblyBuf);
		StrandIterator saveReferenceBegin = referenceBegin;		
		Variant wholeVariant(referenceBegin.GetOriginalPosition(), blockId, referenceBuf, assemblyBuf, assemblySequence, referenceContext, assemblyContext);
		if(referenceBuf == assemblyBuf || !ConfirmVariant(referenceBegin, referenceEnd, assemblyBegin, assemblyEnd))
		{
			return;
		}		

		if(referenceBuf.size() == 0 || assemblyBuf.size() == 0)
		{
			if(referenceBuf.size() + assemblyBuf.size() > 0)
			{
				variantList.push_back(wholeVariant);
			}

			return;
		}
		
		if(referenceBuf.size() * assemblyBuf.size() > (1 << 20))
		{
			variantList.push_back(wholeVariant);
			return;
		}		

		std::vector<Variant> tempList;
		TAlign align;
		resize(rows(align), 2); 
		assignSource(row(align,0), referenceSeq);
		assignSource(row(align,1), assemblySeq);
		int score = globalAlignment(align, Score<int>(1, -1, -1, -1), Hirschberg());		
		TPosition colBegin = beginPosition(cols(align));
		TPosition colEnd = endPosition(cols(align));
		TIterator it1 = iter(row(align, 0), colBegin);
		TIterator it1End = iter(row(align, 0), colEnd);
		TIterator it2 = iter(row(align, 1), colBegin);
		while(it1 != it1End)
		{
			if(!isGap(it1) && !isGap(it2))
			{
				if(*it1 != *it2)
				{					
					tempList.push_back(Variant(referenceBegin.GetOriginalPosition(), blockId, std::string(1, *it1), std::string(1, *it2), assemblySequence, referenceContext, assemblyContext));
				}

				++referenceBegin;
				++it1;
				++it2;
			}
			else
			{
				std::string variantAllele;
				std::string referenceAllele;
				size_t pos = referenceBegin.GetOriginalPosition();
				for(;it1 != it1End; ++it1, ++it2)
				{
					if(isGap(it1))
					{
						variantAllele += *it2;
					}
					else if(isGap(it2))
					{
						++referenceBegin;						
						referenceAllele += *it1;
					}
					else
					{
						break;
					}
				}

				tempList.push_back(Variant(pos, blockId, referenceAllele, variantAllele, assemblySequence, referenceContext, assemblyContext));
			}
		}
		
		if(tempList.size() > 1)
		{
			variantList.push_back(wholeVariant);
		}
		else
		{
			std::copy(tempList.begin(), tempList.end(), std::back_inserter(variantList));
		}
	}

	void VariantCaller::AlignSyntenyBlocks(const BlockInstance & reference, const BlockInstance & assembly, std::vector<Variant> & variantList) const
	{
		size_t pos = reference.GetStart();
		std::vector<std::string> blockSeq(2);
		std::vector<std::vector<Pos> > originalPos(2);
		BlockInstance block[] = {reference, assembly};
		for(size_t i = 0; i < 2; i++)
		{
			std::string::const_iterator begin = block[i].GetChrInstance().GetSequence().begin();
			blockSeq[i].assign(begin + block[i].GetStart(), begin + block[i].GetEnd());
			originalPos[i].resize(block[i].GetEnd() - block[i].GetStart());
			std::generate(originalPos[i].begin(), originalPos[i].end(), Counter<Pos>(static_cast<Pos>(block[i].GetStart())));
		}

		IndexedSequence iseq(blockSeq, originalPos, trimK_, "");
		DNASequence & sequence = iseq.Sequence();
		BifurcationStorage & bifStorage = iseq.BifStorage();
		iseq.ConstructChrIndex();
		StrandIterator referenceIt = sequence.Begin(reference.GetDirection(), 0);
		StrandIterator referenceEnd = sequence.End(reference.GetDirection(), 0);
		StrandIterator assemblyIt = sequence.Begin(assembly.GetDirection(), 1);
		StrandIterator assemblyEnd = sequence.End(assembly.GetDirection(), 1);
		if(NextVertex(iseq, referenceIt, referenceEnd, assemblyIt, assemblyEnd, false))
		{			
			AlignBulgeBranches(reference.GetBlockId(), sequence.Begin(reference.GetDirection(), 0), referenceIt, sequence.Begin(assembly.GetDirection(), 1), assemblyIt, assembly.GetChrInstance(), variantList);
			while(referenceIt != referenceEnd)
			{
				size_t dist = Traverse(referenceIt, assemblyIt) - trimK_;
				std::advance(referenceIt, dist);
				std::advance(assemblyIt, dist);
				StrandIterator nextReferenceIt = referenceIt;
				StrandIterator nextAssemblyIt = assemblyIt;
				NextVertex(iseq, nextReferenceIt, referenceEnd, nextAssemblyIt, assemblyEnd, true);
				AlignBulgeBranches(reference.GetBlockId(), referenceIt, nextReferenceIt, assemblyIt, nextAssemblyIt, assembly.GetChrInstance(), variantList);
				referenceIt = nextReferenceIt;
				assemblyIt = nextAssemblyIt;				
			}
		}
	}

	void VariantCaller::CallVariants(std::vector<Variant> & variantList) const
	{
		variantList.clear();
		const int NO_BLOCK = 0;
		const size_t NO_STAGE = -1;
		typedef std::pair<size_t, int> CoverUnit;
		std::map<CoverUnit, size_t> posInReference;
		std::vector<std::vector<CoverUnit> > cover(chr_->size());
		for(size_t chr = 0; chr < cover.size(); chr++)
		{
			cover[chr].assign((*chr_)[chr].GetSequence().size(), CoverUnit(NO_STAGE, NO_BLOCK));
		}
		
		for(size_t stage = 0; stage < history_.size(); stage++)
		{
			std::vector<IndexPair> group;
			std::vector<BlockInstance> & blockList = history_[stage];
			GroupBy(blockList, compareById, std::back_inserter(group));
			for(std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
			{
				size_t referenceStart = 0;
				size_t inReference = 0;
				size_t inAssembly = 0;
				for(size_t i = it->first; i < it->second; i++) 
				{
					inReference += blockList[i].GetChrId() == refSeqId_ ? 1 : 0;
					inAssembly += blockList[i].GetChrId() != refSeqId_ ? 1 : 0;
				}

				CoverUnit unit(stage, blockList[it->first].GetBlockId());
				if(inReference >= 1 && inAssembly >= 1)
				{
					size_t refPos;
					for(size_t i = it->first; i < it->second; i++)
					{
						size_t start = blockList[i].GetStart();
						size_t end = blockList[i].GetEnd();
						if(blockList[i].GetChrId() == refSeqId_)
						{
							refPos = blockList[i].GetEnd();
						}

						std::vector<CoverUnit>::iterator jt = cover[blockList[i].GetChrId()].begin();
						std::fill(jt + start, jt + end, unit);
						if(blockList[i].GetChrId() == refSeqId_)
						{
							refPos = end;
						}
					}

					if(inReference == 1)
					{
						posInReference[unit] = refPos;
					}
				}
			
				if(inReference == 1 && inAssembly == 1)
				{
					if(stage == history_.size() - 1)
					{
						if(blockList[it->first].GetChrId() != refSeqId_)
						{
							std::swap(blockList[it->first], blockList[it->first + 1]);
						}

						if(blockList[it->first].GetDirection() != DNASequence::positive)
						{					
							blockList[it->first].Reverse();
							blockList[it->first + 1].Reverse();
						}

						AlignSyntenyBlocks(blockList[it->first], blockList[it->first + 1], variantList);
					}
				}
			}
		}		

		for(size_t chr = 0; chr < cover.size(); chr++)
		{
			for(size_t pos = 0; pos < cover[chr].size(); )
			{
				if(cover[chr][pos].first == NO_STAGE)
				{
					size_t end = pos;
					for(; end < cover[chr].size() && cover[chr][end].first == NO_STAGE; end++);
					if(end - pos > minBlockSize_)
					{						
						const std::string & sequence = (*chr_)[chr].GetSequence();
						if(chr == refSeqId_)
						{					
							std::string referenceAllele(sequence.begin() + pos, sequence.begin() + end);
							if(!SearchInAssembly(referenceAllele))
							{
								variantList.push_back(Variant(pos, Variant::UNKNOWN_BLOCK, referenceAllele, "", (*chr_)[chr], "", ""));
							}
						}
						else
						{							
							std::string::const_iterator jt = (*chr_)[chr].GetSequence().begin();
							std::string assemblyAllele(jt + pos, jt + end);
							size_t refPos = (pos > 0 && posInReference.count(cover[chr][pos - 1]) > 0) ? posInReference[cover[chr][pos - 1]] : Variant::UNKNOWN_POS;
							if(!SearchInReference(assemblyAllele))
							{
								variantList.push_back(Variant(refPos, Variant::UNKNOWN_BLOCK, assemblyAllele, "", (*chr_)[chr], "", ""));
							}
						}
					}

					pos = end;
				}
				else
				{
					pos++;
				}
			}
		}

		std::sort(variantList.begin(), variantList.end());
	}

	void VariantCaller::GetHistory(std::vector<std::vector<BlockInstance> > & history) const
	{
		history = history_;
	}
}