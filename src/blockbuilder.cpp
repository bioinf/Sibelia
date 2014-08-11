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
		index_(0), originalChr_(originalChr), tempDir_(tempDir)
	{
	}

	void BlockBuilder::ConstructIndex(size_t k)
	{
		if(index_ == 0)
		{
			std::vector<std::string> record(originalChr_->size());
			for(size_t i = 0; i < record.size(); i++)
			{
				record[i].resize((*originalChr_)[i].GetSequence().size());
				for(size_t j = 0; j < record[i].size(); j++)
				{
					record[i][j] = FastaRecord::IsDefiniteBase(record[i][j]) ? record[i][j] : FastaRecord::DEFINITE_BASE[rand() % FastaRecord::DEFINITE_BASE.size()];
				}
			}

			std::vector<IndexedSequence::BifurcationInstance> bifurcation[2];
			IndexedSequence::EnumerateBifurcationsSArray(record, k, tempDir_, bifurcation[0], bifurcation[1]);
		}
	}

		/*
			for(size_t strand = 0; strand < 2; strand++)
			{
				size_t nowBif = 0;
				DNASequence::Direction dir = static_cast<DNASequence::Direction>(strand);
				for(size_t chr = 0; chr < sequence_->ChrNumber(); chr++)
				{
					size_t pos = 0;
					StrandIterator end = sequence_->End(dir, chr);
					for(DNASequence::StrandIterator it = sequence_->Begin(dir, chr); it != end; ++it, ++pos)
					{
						if(nowBif < bifurcation[strand].size() && chr == bifurcation[strand][nowBif].chr && pos == bifurcation[strand][nowBif].pos)
						{
							bifStorage_->AddPoint(it, bifurcation[strand][nowBif++].bifId);
						}
					}
				}			
			}
		}*/

	void BlockBuilder::Simplify(size_t minBranchSize, size_t maxIterations, ProgressCallBack f)
	{
	}

	void BlockBuilder::GenerateBlocks(std::vector<BlockInstance> & ret, size_t minBlockSize) const
	{
	}
}