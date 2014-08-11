//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "blockbuilder.h"
#include "indexedsequence.h"

namespace SyntenyFinder
{
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
				for(size_t j = 0; j < record[i].size(); j++)
				{
					bool change = !FastaRecord::IsDefiniteBase(record[i][j]) || FastaRecord::IsMaskedBase(record[i][j]);
					record[i][j] = change ? FastaRecord::DEFINITE_BASE[rand() % FastaRecord::DEFINITE_BASE.size()] : record[i][j];
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
							index_->AddEdge(chr, pos, dir, nowBif, record[chr][pos + k], pos);
						}
					}
				}			
			}
		}
		else
		{
			//Reconstruct for new K
		}
	}		

	void BlockBuilder::Simplify(size_t minBranchSize, size_t maxIterations, ProgressCallBack f)
	{
	}

	void BlockBuilder::GenerateBlocks(std::vector<BlockInstance> & ret, size_t minBlockSize) const
	{
	}
}