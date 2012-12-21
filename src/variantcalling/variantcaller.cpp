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
	namespace 
	{
		bool StrandIteratorGEqual(StrandIterator a, StrandIterator b)
		{
			if(a.GetDirection() == DNASequence::positive && b.GetDirection() == DNASequence::positive)
			{
				return a.GetOriginalPosition() >= b.GetOriginalPosition();
			}

			if(a.GetDirection() == DNASequence::negative && b.GetDirection() == DNASequence::negative)
			{
				return a.GetOriginalPosition() <= b.GetOriginalPosition();
			}

			return false;
		}

		size_t StrandIteratorDistance(StrandIterator start, StrandIterator end)
		{
			size_t min = std::min(start.GetOriginalPosition(), end.GetOriginalPosition());
			size_t max = std::max(start.GetOriginalPosition(), end.GetOriginalPosition());
			return max - min;
		}
	}

	VariantCaller::VariantCaller(size_t refSeqId, const std::vector<BlockInstance> & blockList, size_t trimK):
		refSeqId_(refSeqId), blockList_(blockList), trimK_(trimK)
	{

	}

	void VariantCaller::AlignBulgeBranches(StrandIterator referenceBegin, StrandIterator referenceEnd, StrandIterator assemblyBegin, StrandIterator assemblyEnd, std::vector<Variant> & variantList) const
	{
		using namespace seqan;
		typedef String<char>				TSequence;
		typedef Align<TSequence, ArrayGaps>	TAlign;
		typedef Row<TAlign>::Type			TRow;
		typedef Iterator<TRow>::Type		TIterator;
		typedef Position<TAlign>::Type		TPosition;

		std::string buf1(referenceBegin, referenceEnd);
		std::string buf2(assemblyBegin, assemblyEnd);
		TSequence seq1(buf1);
		TSequence seq2(buf2);

		if(buf1.size() * buf2.size() > 1E5 || buf1.size() == 0 || buf2.size() == 0)
		{
			return;
		}		

		TAlign align;
		resize(rows(align), 2); 
		assignSource(row(align,0),seq1);
		assignSource(row(align,1),seq2);
		int score = globalAlignment(align,Score<int>(1,-1,-1,-1), Hirschberg());		
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
					variantList.push_back(Variant(referenceBegin.GetOriginalPosition(), std::string(1, *it1), std::string(1, *it2)));
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

				variantList.push_back(Variant(pos, variantAllele, referenceAllele));
			}
		}

	}

	//To be refactored, code duplication with BlockFinder::GenerateSyntenyBlocks
	void VariantCaller::AlignSyntenyBlocks(const BlockInstance & reference, const BlockInstance & assembly, std::vector<Variant> & variantList) const
	{
		size_t pos = reference.GetStart();
		std::vector<size_t> chrStart(2);
		std::vector<std::string> blockSeq(2);
		std::vector<std::vector<Pos> > chrPos(2);
		BlockInstance block[] = {reference, assembly};
		for(size_t i = 0; i < 2; i++)
		{
			std::string::const_iterator begin = block[i].GetChrInstance().GetSequence().begin();
			blockSeq[i].assign(begin + block[i].GetStart(), begin + block[i].GetEnd());
			for(size_t j = 0; j < blockSeq[i].size(); j++)
			{
				chrPos[i].push_back(static_cast<Pos>(pos++));
			}

			chrStart[i] = chrPos[i][0];
		}
		
		const size_t oo = UINT_MAX;
		std::auto_ptr<DNASequence> sequence;
		std::auto_ptr<BifurcationStorage> bifStorage;
		{
			std::vector<std::vector<BlockFinder::BifurcationInstance> > bifurcation(2);	
			size_t maxId = BlockFinder::EnumerateBifurcationsSArrayInRAM(blockSeq, trimK_, bifurcation[0], bifurcation[1]);
			bifStorage.reset(new BifurcationStorage(maxId));
			sequence.reset(new DNASequence(blockSeq, chrPos, true));
			BlockFinder::ConstructBifStorage(*sequence, bifurcation, *bifStorage);
		}

		size_t refChr = 0;
		IteratorProxyVector startKMer;
		StrandIterator referenceStart = sequence->Begin(reference.GetDirection(), 0);
		StrandIterator referenceEnd = sequence->End(reference.GetDirection(), 0);
		StrandIterator assemblyStart = sequence->Begin(reference.GetDirection(), 1);
		StrandIterator assemblyEnd = sequence->End(reference.GetDirection(), 1);
		while(referenceStart != referenceEnd)
		{
			size_t minSum = oo;
			StrandIterator nextAssemblyStart = assemblyStart;
			StrandIterator nextReferenceStart = referenceStart;
			StrandIterator referenceProbe = referenceStart;
			for(size_t referenceStep = 0; referenceProbe != referenceEnd && referenceStep < minSum; ++referenceProbe, ++referenceStep)
			{
				size_t bifId = bifStorage->GetBifurcation(referenceProbe);
				if(bifId != BifurcationStorage::NO_BIFURCATION)
				{
					startKMer.clear();
					bifStorage->ListPositions(bifId, std::back_inserter(startKMer));
					for(size_t pos = 0; pos < startKMer.size(); pos++)
					{
						StrandIterator kmer = *startKMer[pos];
						size_t bifChr = std::upper_bound(chrStart.begin(), chrStart.end(), kmer.GetOriginalPosition()) - chrStart.begin() - 1;
						if(bifChr != refChr)
						{							
							size_t assemblyStep = StrandIteratorDistance(kmer, assemblyStart);
							if(assemblyStep + referenceStep < minSum && StrandIteratorGEqual(kmer, assemblyStart))
							{								
								nextAssemblyStart = kmer;
								nextReferenceStart = referenceProbe;
								minSum = assemblyStep + referenceStep;
							}
						}
					}
				}
			}

			if(minSum == oo)
			{
				nextAssemblyStart = assemblyEnd;
				nextReferenceStart = referenceEnd;
			}

			AlignBulgeBranches(referenceStart, nextReferenceStart, assemblyStart, nextAssemblyStart, variantList);
			referenceStart = nextReferenceStart;
			assemblyStart = nextAssemblyStart;
			while(referenceStart != referenceEnd && assemblyStart != assemblyEnd && *referenceStart == *assemblyStart)
			{
				++referenceStart;
				++assemblyStart;
			}
		}
	}

	void VariantCaller::CallVariants(std::vector<Variant> & variantList) const
	{
		variantList.clear();
		std::vector<IndexPair> group;
		GroupBy(blockList_, compareById, std::back_inserter(group));
		for(std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
		{
			size_t inReference = 0;
			size_t inAssembly = 0;std::cout << it - group.begin() << std::endl;
			for(size_t i = it->first; i < it->second; i++)
			{
				inReference += blockList_[i].GetChrId() == refSeqId_ ? 1 : 0;
				inAssembly += blockList_[i].GetChrId() != refSeqId_ ? 1 : 0;
			}

			if(inReference == 1 && inAssembly == 1)
			{
				if(blockList_[it->first].GetChrId() != refSeqId_)
				{
					std::swap(blockList_[it->first], blockList_[it->first + 1]);
				}

				AlignSyntenyBlocks(blockList_[it->first], blockList_[it->first + 1], variantList);
			}
		}

		std::sort(variantList.begin(), variantList.end());
	}
}