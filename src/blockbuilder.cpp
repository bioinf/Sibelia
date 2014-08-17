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
				virtualChrSize_.push_back(record[i].size());
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
							char mark = pos + k < record[chr].size() ? record[chr][pos + k] : IndexedSequence::SEPARATION_CHAR;
							index_->AddEdge(chr, pos, dir, nowBif++, mark, pos);
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

	void BlockBuilder::WriteIndexToDot(std::ostream & out) const
	{
		out << "digraph G\n{\nrankdir=LR" << std::endl;
		for(size_t chr = 0; chr < originalChr_->size(); chr++)
		{
			for(size_t strand = 0; strand < 2; strand++)
			{
				FastaRecord::Direction dir = static_cast<FastaRecord::Direction>(strand);
				DeBruijnIndex::Edge prevEdge = index_->GetEdgeAtPosition(chr, 0, dir);
				for(size_t pos = 1; pos < virtualChrSize_.size(); ++pos)
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

	void BlockBuilder::Simplify(size_t minBranchSize, size_t maxIterations, ProgressCallBack f)
	{

	}

	void BlockBuilder::GenerateBlocks(std::vector<BlockInstance> & ret, size_t minBlockSize) const
	{
	}
}