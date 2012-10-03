#include "blockfinder.h"

namespace SyntenyFinder
{
	BlockFinder::BlockFinder(const std::vector<FASTARecord> & chrList): chrList_(chrList)
	{
		originalPos_.resize(chrList.size());
		for(size_t i = 0; i < originalPos_.size(); i++)
		{
			originalPos_[i].resize(chrList[i].sequence.size());
			for(size_t j = 0; j < originalPos_[i].size(); j++)
			{
				originalPos_[i][j] = static_cast<Pos>(j);
			}
		}
	}

	void BlockFinder::PerformGraphSimplifications(size_t k, size_t minBranchSize, size_t maxIterations, ProgressCallBack f)
	{
		BifurcationStorage bifStorage;
		DNASequence sequence(chrList_, originalPos_);
		EnumerateBifurcations(sequence, bifStorage, k, f);
		SimplifyGraph(sequence, bifStorage, k, minBranchSize, maxIterations, f);
		for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
		{
			originalPos_[chr].clear();
			chrList_[chr].sequence.clear();
			StrandIterator end = sequence.PositiveEnd(chr);
			for(StrandIterator it = sequence.PositiveBegin(chr); it != end; ++it)
			{
				chrList_[chr].sequence.push_back(*it);
				originalPos_[chr].push_back(static_cast<Pos>(it.GetOriginalPosition()));
			}
		}
	}

	void BlockFinder::GenerateSyntenyBlocks(size_t k, std::vector<BlockInstance> & block, ProgressCallBack enumeration) const
	{
		block.clear();
		BifurcationStorage bifStorage;
		DNASequence sequence(chrList_, originalPos_);
		EnumerateBifurcations(sequence, bifStorage, k, enumeration);
		ConvertEdgesToBlocks(sequence, bifStorage, k, block);
	}
}