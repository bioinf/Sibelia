//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wc++11-extensions"

#include "postprocessor.h"

namespace SyntenyFinder
{
	namespace
	{
		const size_t MAX_CORRECTION_RANGE = 1 << 10;

		struct Stripe
		{
			int firstBlock;
			int secondBlock;
			Stripe() {}
			Stripe(int firstBlock, int secondBlock): firstBlock(firstBlock), secondBlock(secondBlock) {}
			bool operator < (const Stripe & toCompare) const
			{
				return firstBlock < toCompare.firstBlock;
			}
		};

		std::string ReadInReverseDirection(const std::string & str)
		{
			std::string::const_reverse_iterator it1 = str.rbegin();
			std::string::const_reverse_iterator it2 = str.rend();
			return std::string(CFancyIterator(it1, DNASequence::Translate, ' '), CFancyIterator(it2, DNASequence::Translate, ' '));
		}
	}

	void Postprocessor::GlueStripes(std::vector<BlockInstance> & block)
	{
		std::vector<std::vector<BlockInstance> > perm(chr_->size());
		for(size_t i = 0; i < block.size(); i++)
		{
			perm[block[i].GetChrId()].push_back(block[i]);
		}

		for(size_t i = 0; i < perm.size(); i++)
		{
			std::sort(perm[i].begin(), perm[i].end(), compareByStart);
		}

		int sentinel = INT_MAX >> 1;
		bool glue = false;
		do
		{
			std::vector<Stripe> stripe;
			for(size_t chr = 0; chr < perm.size(); chr++)
			{
				for(size_t i = 0; i < perm[chr].size(); i++)
				{
					int bid = perm[chr][i].GetSignedBlockId();
					if(bid > 0)
					{
						int nextBid = i < perm[chr].size() - 1 ? perm[chr][i + 1].GetSignedBlockId() : sentinel;
						stripe.push_back(Stripe(bid, nextBid));
					}
					else
					{
						int prevBid = i > 0 ? perm[chr][i - 1].GetSignedBlockId() : -sentinel;
						stripe.push_back(Stripe(-bid, -prevBid));
					}
				}
			}

			size_t now = 0;
			size_t next = 0;
			std::sort(stripe.begin(), stripe.end());
			for(; now < stripe.size(); now = next)
			{
				glue = true;
				for(; next < stripe.size() && stripe[next].firstBlock == stripe[now].firstBlock; next++)
				{
					if(stripe[next].secondBlock != stripe[now].secondBlock || stripe[next].secondBlock == sentinel || Abs(stripe[next].secondBlock) == stripe[next].firstBlock)
					{
						glue = false;
					}
				}

				if(glue)
				{
					typedef std::vector<Stripe>::iterator It;
					std::pair<It, It> range = std::equal_range(stripe.begin(), stripe.end(), Stripe(Abs(stripe[now].secondBlock), 0));
					if(range.second - range.first != next - now)
					{
						glue = false;
					}
					else
					{
						break;
					}
				}
			}

			if(glue)
			{
				assert(next - now > 1);
				int glueBid = stripe[now].firstBlock;
				for(size_t chr = 0; chr < perm.size(); chr++)
				{
					for(size_t i = 0; i < perm[chr].size(); i++)
					{
						int bid = perm[chr][i].GetBlockId();
						if(bid == glueBid)
						{
							bid = perm[chr][i].GetSignedBlockId();
							if(bid > 0)
							{
								BlockInstance & a = perm[chr][i];
								BlockInstance & b = perm[chr][i + 1];
								a = BlockInstance(a.GetSignedBlockId(), &a.GetChrInstance(), a.GetStart(), b.GetEnd());
								perm[chr].erase(perm[chr].begin() + i + 1);
							}
							else
							{
								BlockInstance & a = perm[chr][--i];
								BlockInstance & b = perm[chr][i + 1];
								a = BlockInstance(b.GetSignedBlockId(), &a.GetChrInstance(), a.GetStart(), b.GetEnd());
								perm[chr].erase(perm[chr].begin() + i + 1);
							}
						}
					}
				}
			}
		}
		while(glue);

		block.clear();
		std::vector<int> oldId;
		for(size_t chr = 0; chr < perm.size(); chr++)
		{
			for(size_t i = 0; i < perm[chr].size(); i++)
			{
				block.push_back(perm[chr][i]);
				oldId.push_back(perm[chr][i].GetBlockId());
			}
		}

		std::sort(oldId.begin(), oldId.end());
		oldId.erase(std::unique(oldId.begin(), oldId.end()), oldId.end());
		for(std::vector<BlockInstance>::iterator it = block.begin(); it != block.end(); ++it)
		{
			int sign = it->GetSignedBlockId() > 0 ? +1 : -1;
			size_t newId = std::lower_bound(oldId.begin(), oldId.end(), it->GetBlockId()) - oldId.begin() + 1;
			*it = BlockInstance(static_cast<int>(newId) * sign, &it->GetChrInstance(), it->GetStart(), it->GetEnd());
		}
	}

	typedef std::vector<BlockInstance>::const_iterator BLCIterator;

	Postprocessor::Postprocessor(const std::vector<FASTARecord> & chr, size_t minBlockSize):
		chr_(&chr), minBlockSize_(minBlockSize), correctionRange_(std::min(minBlockSize, MAX_CORRECTION_RANGE))
	{
	}

	const BlockInstance* Postprocessor::PreviousBlock(const BlockInstance & block, const std::vector<BlockInstance> & blockList)
	{
		const BlockInstance* ret = 0;
		size_t start = block.GetStart();
		for(size_t i = 0; i < blockList.size(); i++)
		{
			if(blockList[i] != block && blockList[i].GetChrId() == block.GetChrId() && blockList[i].GetEnd() <= start)
			{
				if(ret == 0 || start - blockList[i].GetEnd() < start - ret->GetEnd())
				{
					ret = &blockList[i];
				}
			}
		}

		return ret;
	}

	const BlockInstance* Postprocessor::NextBlock(const BlockInstance & block, const std::vector<BlockInstance> & blockList)
	{
		const BlockInstance* ret = 0;
		size_t end = block.GetEnd();
		for(size_t i = 0; i < blockList.size(); i++)
		{
			if(blockList[i] != block && blockList[i].GetChrId() == block.GetChrId() && blockList[i].GetStart() >= end)
			{
				if(ret == 0 || blockList[i].GetStart() - end < ret->GetStart() - end)
				{
					ret = &blockList[i];
				}
			}
		}

		return ret;
	}
	
	std::pair<size_t, size_t> Postprocessor::DetermineLeftProbableBoundaries(std::vector<BlockInstance> & blockList, size_t blockid)
	{
		std::pair<size_t, size_t> ret;
		BlockInstance & block = blockList[blockid];
		size_t nowStart = block.GetStart();		
		ret.second = block.GetStart() + correctionRange_;
		const BlockInstance * previousBlock = PreviousBlock(block, blockList);
		if(previousBlock != 0)
		{	
			size_t previousEnd = previousBlock->GetEnd();
			ret.first = std::max(previousEnd, nowStart - correctionRange_) + 1;			
		}
		else
		{
			ret.first = nowStart >= correctionRange_ ? nowStart - correctionRange_ + 1: 0;
		}
				
		return ret;
	}

	std::pair<size_t, size_t> Postprocessor::DetermineRightProbableBoundaries(std::vector<BlockInstance> & blockList, size_t blockid)
	{
		std::pair<size_t, size_t> ret;
		BlockInstance & block = blockList[blockid];
		size_t nowEnd = block.GetEnd();
		ret.first = block.GetEnd() - correctionRange_ + 1;
		const BlockInstance * nextBlock = NextBlock(block, blockList);
		if(nextBlock != 0)
		{			
			size_t nextStart = nextBlock->GetStart();
			ret.second = std::min(nextStart, nowEnd + correctionRange_);
		}
		else
		{
			size_t chrSize = block.GetChrInstance().GetSequence().size();
			ret.second = nowEnd + correctionRange_ < chrSize ? nowEnd + correctionRange_ : chrSize;
		}
				
		return ret;
	}

	void Postprocessor::GetBoundariesSequence(const BlockInstance & block, std::pair<size_t, size_t> leftBoundaries, std::pair<size_t, size_t> rightBoundaries, std::string & start, std::string & end)
	{		
		const std::string::const_iterator & chr = block.GetChrInstance().GetSequence().begin();
		if(block.GetDirection() == DNASequence::positive)
		{			
			start.assign(chr + leftBoundaries.first, chr + leftBoundaries.second);
			end.assign(chr + rightBoundaries.first, chr + rightBoundaries.second);
		}
		else
		{
			start.assign(chr + rightBoundaries.first, chr + rightBoundaries.second);
			end.assign(chr + leftBoundaries.first, chr + leftBoundaries.second);
			start = ReadInReverseDirection(start);
			end = ReadInReverseDirection(end);
		}
	}

	void Postprocessor::LocalAlignment(const std::string & sequence1, const std::string & sequence2, std::pair<size_t, size_t> & coord1, std::pair<size_t, size_t> & coord2)
	{
		using namespace seqan;
		typedef String<char>				TSequence;	// sequence type
		typedef Align<TSequence, ArrayGaps>	TAlign;		// align type
		typedef Row<TAlign>::Type			TRow;
		typedef Iterator<TRow>::Type		TIterator;
		typedef Position<TAlign>::Type		TPosition;
		
		TSequence seq1 = sequence1;
		TSequence seq2 = sequence2;
		TAlign align;
		resize(rows(align), 2); 
		assignSource(row(align, 0), seq1);
		assignSource(row(align, 1), seq2);
		localAlignment(align, Score<int>(38, -25, -25));
		coord1.first = clippedBeginPosition(row(align, 0));
		coord1.second = clippedEndPosition(row(align, 0));
		coord2.first = clippedBeginPosition(row(align, 1));
		coord2.second = clippedEndPosition(row(align, 1));
	}

	void Postprocessor::UpdateBlockBoundaries(BlockInstance & block, std::pair<size_t, size_t> leftBoundaries, std::pair<size_t, size_t> rightBoundaries, std::pair<size_t, size_t> startAlignmentCoords, std::pair<size_t, size_t> endAlignmentCoords)
	{
		if(block.GetDirection() == DNASequence::positive)
		{
			size_t newStart = leftBoundaries.first + startAlignmentCoords.first;
			size_t newEnd = rightBoundaries.first + endAlignmentCoords.second;
			block = BlockInstance(block.GetSignedBlockId(), &block.GetChrInstance(), newStart, newEnd);
		}
		else
		{			
			size_t newStart = leftBoundaries.second - endAlignmentCoords.second;
			size_t newEnd = rightBoundaries.second - startAlignmentCoords.first;
			block = BlockInstance(block.GetSignedBlockId(), &block.GetChrInstance(), newStart, newEnd);
		}
	}

	void Postprocessor::CorrectBlocksBoundaries(std::vector<BlockInstance> & blockList, size_t referenceBlock, size_t assemblyBlock)
	{		
		std::pair<size_t, size_t> referenceStartCoord;
		std::pair<size_t, size_t> referenceEndCoord;
		std::pair<size_t, size_t> assemblyStartCoord;
		std::pair<size_t, size_t> assemblyEndCoord;
		std::string referenceStart;
		std::string referenceEnd;
		std::string assemblyStart;
		std::string assemblyEnd;
		std::pair<size_t, size_t> referenceLeftBoundaries = DetermineLeftProbableBoundaries(blockList, referenceBlock);
		std::pair<size_t, size_t> referenceRightBoundaries = DetermineRightProbableBoundaries(blockList, referenceBlock);
		std::pair<size_t, size_t> assemblyLeftBoundaries = DetermineLeftProbableBoundaries(blockList, assemblyBlock);
		std::pair<size_t, size_t> assemblyRightBoundaries = DetermineRightProbableBoundaries(blockList, assemblyBlock);
		GetBoundariesSequence(blockList[referenceBlock], referenceLeftBoundaries, referenceRightBoundaries, referenceStart, referenceEnd);
		GetBoundariesSequence(blockList[assemblyBlock], assemblyLeftBoundaries, assemblyRightBoundaries, assemblyStart, assemblyEnd);		
		LocalAlignment(referenceStart, assemblyStart, referenceStartCoord, assemblyStartCoord);
		LocalAlignment(referenceEnd, assemblyEnd, referenceEndCoord, assemblyEndCoord);
		UpdateBlockBoundaries(blockList[referenceBlock], referenceLeftBoundaries, referenceRightBoundaries, referenceStartCoord, referenceEndCoord);
		UpdateBlockBoundaries(blockList[assemblyBlock], assemblyLeftBoundaries, assemblyRightBoundaries, assemblyStartCoord, assemblyEndCoord);
	}

	bool Postprocessor::ImproveBlockBoundaries(std::vector<BlockInstance> & blockList)
	{
		bool ret = false;
		std::vector<IndexPair> group;
		std::vector<BlockInstance> newBlockList(blockList);
		GroupBy(blockList, compareById, std::back_inserter(group));		
		for(std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
		{		
			for(size_t i = it->first; i < it->second; i++)
			{
				std::vector<BlockInstance> tempBlock(2, blockList[i]);
				for(size_t j = it->first; j < it->second; j++)
				{
					if(i != j)
					{
						tempBlock[1] = blockList[j];
						CorrectBlocksBoundaries(tempBlock, 0, 1);
						if(tempBlock[0].GetLength() > newBlockList[i].GetLength())
						{
							ret = true;
							newBlockList[i] = tempBlock[0];
						}
					}
				}
			}
		}

		blockList = newBlockList;
		return true;
	}

	void Postprocessor::MatchRepeats(std::vector<BlockInstance> & blockList, const std::set<size_t> & referenceSequenceId)
	{
		bool glue = true;
		std::map<int, size_t> inReference;
		std::map<int, std::vector<size_t> > multiplicity;
		while(glue)
		{
			glue = false;
			inReference.clear();
			multiplicity.clear();
			std::sort(blockList.begin(), blockList.end(), CompareBlocksNaturally);
			for(std::vector<BlockInstance>::iterator it = blockList.begin(); it != blockList.end(); ++it)
			{
				if(multiplicity.find(it->GetBlockId()) == multiplicity.end())
				{
					inReference.insert(std::make_pair(it->GetBlockId(), 0));
					multiplicity.insert(std::make_pair(it->GetBlockId(), std::vector<size_t>()));					
				}

				multiplicity[it->GetBlockId()].push_back(it - blockList.begin());
				if(referenceSequenceId.count(it->GetChrId()) > 0)
				{
					inReference[it->GetBlockId()]++;
				}
			}

			for(std::map<int, std::vector<size_t> >::iterator it = multiplicity.begin(); it != multiplicity.end() && !glue; ++it)
			{
				if(it->second.size() == 2 && inReference[it->first] == 1)
				{
					std::vector<size_t> nextBlockPos;
					for(std::vector<size_t>::iterator posIt = it->second.begin(); posIt != it->second.end(); ++posIt)
					{
						int pos = static_cast<int>(*posIt);
						int nextPos = pos + blockList[pos].GetSign();
						if(nextPos >= 0 && nextPos < static_cast<int>(blockList.size()) && blockList[pos].GetChrId() == blockList[nextPos].GetChrId())
						{
							nextBlockPos.push_back(nextPos);
						}
					}

					if(nextBlockPos.size() == 2)
					{
						BlockInstance startBlock[2] = {blockList[it->second[0]], blockList[it->second[1]]};
						BlockInstance nextBlock[2] = {blockList[nextBlockPos[0]], blockList[nextBlockPos[1]]};
						int startMlp = startBlock[0].GetSign() * startBlock[1].GetSign();
						int nextMlp = nextBlock[0].GetSign() * nextBlock[1].GetSign();
						if(nextBlock[0].GetBlockId() == nextBlock[1].GetBlockId() && startMlp == nextMlp)
						{
							glue = true;
							for(size_t i = 0; i < 2; i++)
							{
								blockList.erase(std::remove_if(blockList.begin(), blockList.end(), boost::bind(&BlockInstance::operator==, startBlock[i], _1)), blockList.end());
								blockList.erase(std::remove_if(blockList.begin(), blockList.end(), boost::bind(&BlockInstance::operator==, nextBlock[i], _1)), blockList.end());
								size_t start = std::min(startBlock[i].GetStart(), nextBlock[i].GetStart());
								size_t end = std::max(startBlock[i].GetEnd(), nextBlock[i].GetEnd());
								blockList.push_back(BlockInstance(startBlock[i].GetSignedBlockId(), &startBlock[i].GetChrInstance(), start, end));
							}
						}
					}
				}
			}
		}

		std::vector<BlockInstance> toDelete;
		for(std::map<int, std::vector<size_t> >::iterator it = multiplicity.begin(); it != multiplicity.end(); ++it)
		{
			if(it->second.size() == 1)
			{
				toDelete.push_back(blockList[it->second[0]]);
			}
		}

		for(std::vector<BlockInstance>::iterator it = toDelete.begin(); it != toDelete.end(); ++it)
		{
			blockList.erase(std::remove_if(blockList.begin(), blockList.end(), boost::bind(&BlockInstance::operator==, *it, _1)), blockList.end());
		}
		
		std::vector<int> oldId;
		for(std::vector<BlockInstance>::iterator it = blockList.begin(); it != blockList.end(); ++it)
		{
			oldId.push_back(it->GetBlockId());
		}

		std::sort(oldId.begin(), oldId.end());
		oldId.erase(std::unique(oldId.begin(), oldId.end()), oldId.end());
		for(std::vector<BlockInstance>::iterator it = blockList.begin(); it != blockList.end(); ++it)
		{
			int sign = it->GetSign();
			size_t newId = std::lower_bound(oldId.begin(), oldId.end(), it->GetBlockId()) - oldId.begin() + 1;
			*it = BlockInstance(static_cast<int>(newId) * sign, &it->GetChrInstance(), it->GetStart(), it->GetEnd());
		}
	}
}