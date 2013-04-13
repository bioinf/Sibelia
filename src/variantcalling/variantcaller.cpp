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

	bool VariantCaller::ConfirmVariant(StrandIterator referenceStart, StrandIterator referenceEnd, StrandIterator assemblyStart, StrandIterator assemblyEnd) const
	{
		size_t confirmDist = trimK_ * 2;
		for(size_t i = 0; i < confirmDist; i++)
		{
			--referenceStart;
			--assemblyStart;
			if(!referenceStart.AtValidPosition() || !assemblyStart.AtValidPosition() || *referenceStart != *assemblyStart)
			{
				++referenceStart;
				++assemblyStart;
			}
		}
		
		size_t forward = Traverse(referenceEnd, assemblyEnd, confirmDist);
		std::advance(referenceEnd, forward);
		std::advance(assemblyEnd, forward);
		const std::string & sequence = refSeq_.GetSequence();		
		std::vector<sauchar_t> pattern(assemblyStart, assemblyEnd);
		saidx_t count = sa_search(reinterpret_cast<const sauchar_t*>(sequence.c_str()), sequence.size(), &pattern[0], pattern.size(), &suffixArray_[0], suffixArray_.size(), &indexOut_[0]);
		return count == 0;
	}

	VariantCaller::VariantCaller(const FASTARecord & refSeq, const std::vector<BlockInstance> & blockList, size_t trimK):
		refSeq_(refSeq), blockList_(blockList), trimK_(trimK)
	{
		const std::string & sequence = refSeq.GetSequence();
		suffixArray_.resize(sequence.size());
		indexOut_.resize(sequence.size());
		divsufsort(reinterpret_cast<const sauchar_t*>(sequence.c_str()), &suffixArray_[0], static_cast<saidx_t>(sequence.size()));
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
		bool collinear = referenceBegin.GetDirection() == assemblyBegin.GetDirection();
		Variant wholeVariant(referenceBegin.GetOriginalPosition(), blockId, collinear, referenceBuf, assemblyBuf, assemblySequence, referenceContext, assemblyContext);
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
		int score = globalAlignment(align, Score<int>(1,-1,-1,-1), Hirschberg());		
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
					tempList.push_back(Variant(referenceBegin.GetOriginalPosition(), blockId, collinear, std::string(1, *it1), std::string(1, *it2), assemblySequence, referenceContext, assemblyContext));
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

				tempList.push_back(Variant(pos, blockId, collinear, referenceAllele, variantAllele, assemblySequence, referenceContext, assemblyContext));
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
			AlignBulgeBranches(reference.GetBlockId(), sequence.Begin(reference.GetDirection(), 0), referenceIt,
				sequence.Begin(assembly.GetDirection(), 1), assemblyIt, assembly.GetChrInstance(), variantList);
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
		std::vector<IndexPair> group;
		GroupBy(blockList_, compareById, std::back_inserter(group));
		size_t refSeqId = refSeq_.GetId();
		for(std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
		{
			size_t inReference = 0;
			size_t inAssembly = 0;
			for(size_t i = it->first; i < it->second; i++)
			{
				inReference += blockList_[i].GetChrId() == refSeqId ? 1 : 0;
				inAssembly += blockList_[i].GetChrId() != refSeqId ? 1 : 0;
			}

			if(inReference == 1 && inAssembly == 1)
			{
				if(blockList_[it->first].GetChrId() != refSeqId)
				{
					std::swap(blockList_[it->first], blockList_[it->first + 1]);
				}

				if(blockList_[it->first].GetDirection() != DNASequence::positive)
				{					
					blockList_[it->first].Reverse();
					blockList_[it->first + 1].Reverse();
				}

				AlignSyntenyBlocks(blockList_[it->first], blockList_[it->first + 1], variantList);
			}
		}

		std::sort(variantList.begin(), variantList.end());
	}

	void VariantCaller::CallRearrangements(std::vector<Reversal> & reversal, std::vector<Translocation> & translocation) const
	{
		std::vector<IndexPair> group;
		std::vector<BlockInstance> reversed;
		std::vector<BlockInstance>::iterator listBegin = blockList_.begin();
		GroupBy(blockList_, compareByChrId, std::back_inserter(group));
		std::sort(listBegin + group[0].first, listBegin + group[0].second, CompareBlocksNaturally);
		const BLCIterator referenceBegin = listBegin + group[0].first;
		const BLCIterator referenceEnd = listBegin + group[0].second;
		for(size_t i = 1; i < group.size(); i++)
		{
			std::vector<BlockInstance>::iterator begin = listBegin + group[i].first;
			std::vector<BlockInstance>::iterator end = listBegin + group[i].second;
			if(begin == end)
			{
				continue;
			}

			std::sort(begin, end, CompareBlocksNaturally);
			size_t size = end - begin; 
			reversed.resize(size);			
			Reverse(begin, end, reversed.begin());
			if(Apply(referenceBegin, referenceEnd, begin, end) == referenceEnd && Apply(referenceBegin, referenceEnd, reversed.begin(), reversed.end()) == referenceEnd)
			{
				size_t fragmentSize = size;
				bool detected = false;
				do
				{
					std::vector<BlockInstance> reversedFragment(fragmentSize);
					for(size_t pos = 0; pos < size - fragmentSize + 1 && !detected; ++pos)
					{
						std::pair<BLCIterator, BLCIterator> range[] = 
						{
							std::make_pair(begin + pos, begin + pos + fragmentSize),
							std::make_pair(reversed.begin() + pos, reversed.begin() + pos + fragmentSize)
						};

						std::pair<BLCIterator, BLCIterator> bound[] = 
						{
							std::make_pair(begin, end),
							std::make_pair(reversed.begin(), reversed.end())
						};

						for(size_t p = 0; p < 2 && !detected; p++)
						{
							Reverse(range[p].first, range[p].second, reversedFragment.begin());
							BLCIterator place = Apply(referenceBegin, referenceEnd, reversedFragment.begin(), reversedFragment.end());
 							if(place != referenceEnd)
							{
								int prevReferenceId = place == referenceBegin ? 0 : (place - 1)->GetSignedBlockId();
								int prevAssemblyId = range[p].first == bound[p].first ? 0 : (range[p].first - 1)->GetSignedBlockId();
								int nextReferenceId = place + fragmentSize >= referenceEnd ? 0 : (place + fragmentSize)->GetSignedBlockId();
								int nextAssemblyId = range[p].first + fragmentSize >= bound[p].second ? 0 : (range[p].first + fragmentSize)->GetSignedBlockId(); 
								if((prevReferenceId == prevAssemblyId || prevAssemblyId == 0) && (nextReferenceId == nextAssemblyId || nextAssemblyId == 0))
								{									
									reversal.push_back(Reversal(place->GetStart(), (place + fragmentSize - 1)->GetEnd()));
									detected = true;
								}
							}
						}
					}
				}
				while(--fragmentSize > 0 && !detected);

				if(!detected)
				{
					size_t fragmentSize = size;
					do
					{
						std::vector<BlockInstance> reversedFragment(fragmentSize);
						for(size_t pos = 0; pos < size - fragmentSize + 1 && !detected; ++pos)
						{
							std::pair<BLCIterator, BLCIterator> range[] = 
							{
								std::make_pair(begin + pos, begin + pos + fragmentSize),
								std::make_pair(reversed.begin() + pos, reversed.begin() + pos + fragmentSize)
							};

							for(size_t p = 0; p < 2 && !detected; p++)
							{
								BLCIterator place = Apply(referenceBegin, referenceEnd, range[p].first, range[p].second);
								if(place != referenceEnd)
								{
									detected = true;
									translocation.push_back(Translocation(place->GetStart(), (place + fragmentSize - 1)->GetEnd(), 0));
								}
							}
						}
					}
					while(--fragmentSize > 0 && !detected);
				}
			}
		}
	}

	void VariantCaller::GetBlockList(std::vector<BlockInstance> & blockList) const
	{
		blockList = blockList_;
	}
}