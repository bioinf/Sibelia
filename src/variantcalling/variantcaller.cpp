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
		const size_t oo = UINT_MAX;
		bool NextVertex(IndexedSequence & iseq, StrandIterator & reference, StrandIterator referenceEnd, StrandIterator & assembly, StrandIterator assemblyEnd, bool strict)
		{
			size_t minSum = oo;						
			IteratorProxyVector startKMer;
			StrandIterator nowReference = reference;
			StrandIterator nextReference = referenceEnd;
			StrandIterator nextAssembly = assemblyEnd;	
			size_t assemblyChr = iseq.GetChr(assembly);
			for(size_t step = 0; step < minSum && nowReference.AtValidPosition(); ++nowReference)
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

		size_t Traverse(StrandIterator reference, StrandIterator assembly)
		{
			size_t ret = 0;
			for(;reference.AtValidPosition() && assembly.AtValidPosition() && *reference == *assembly; ++ret)
			{
				++reference;
				++assembly;
			}

			return ret;
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
		if(buf1.size() == 0 || buf2.size() == 0)
		{
			if(buf1.size() + buf2.size() > 0)
			{
				variantList.push_back(Variant(referenceBegin.GetOriginalPosition(), buf1, buf2));
			}

			return;
		}

		if(buf1.size() * buf2.size() > 1E6)
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

				variantList.push_back(Variant(pos, referenceAllele, variantAllele));
			}
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
			AlignBulgeBranches(sequence.Begin(reference.GetDirection(), 0), referenceIt,
				sequence.Begin(assembly.GetDirection(), 1), assemblyIt, variantList);
			while(referenceIt != referenceEnd)
			{
				size_t dist = Traverse(referenceIt, assemblyIt) - trimK_;
				std::advance(referenceIt, dist);
				std::advance(assemblyIt, dist);
				StrandIterator nextReferenceIt = referenceIt;
				StrandIterator nextAssemblyIt = assemblyIt;
				NextVertex(iseq, nextReferenceIt, referenceEnd, nextAssemblyIt, assemblyEnd, true);
				AlignBulgeBranches(referenceIt, nextReferenceIt, assemblyIt, nextAssemblyIt, variantList);
				referenceIt = nextReferenceIt;
				assemblyIt = nextAssemblyIt;
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
			size_t inAssembly = 0;
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