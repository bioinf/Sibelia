//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "indexedsequence.h"

namespace SyntenyFinder
{
	IndexedSequence::IndexedSequence(const std::vector<std::string> & record, std::vector<std::vector<Pos> > & originalPos, size_t k, const std::string & tempDir, bool clear): k_(k)
	{
		size_t maxId;
		std::vector<std::vector<BifurcationInstance> > bifurcation(2);	
		if(tempDir.size() == 0)
		{
			maxId = EnumerateBifurcationsSArrayInRAM(record, bifurcation[0], bifurcation[1]);
		}
		else
		{
			maxId = EnumerateBifurcationsSArray(record, tempDir, bifurcation[0], bifurcation[1]);
		}

		bifStorage_.reset(new BifurcationStorage(maxId));
		sequence_.reset(new DNASequence(record, originalPos, clear));
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

	#ifdef _DEBUG
		bifStorage_->FormDictionary(idMap_, k_);
	#endif
	}

#ifdef _DEBUG	
	size_t IndexedSequence::GetMustBeBifurcation(StrandIterator jt)
	{		
		std::string buf(std::string(jt, AdvanceForward(jt, k_)));
		KMerBifMap::iterator kt = idMap_.find(buf);
		return kt == idMap_.end() ? BifurcationStorage::NO_BIFURCATION : kt->second;
	}

	void IndexedSequence::Test()
	{	
		for(size_t strand = 0; strand < 2; strand++)
		{
			for(size_t chr = 0; chr < sequence_->ChrNumber(); chr++)
			{
				StrandIterator begin = sequence_->Begin((DNASequence::Direction)strand, chr);
				StrandIterator end = sequence_->End((DNASequence::Direction)strand, chr);
				SlidingWindow<StrandIterator> window(begin, end, k_);
				for(; window.Valid(); window.Move())
				{
					StrandIterator jt = window.GetBegin();
					size_t pos = sequence_->GlobalIndex(jt);
					std::string buf(std::string(jt, AdvanceForward(jt, k_)));
					size_t actualBifurcation = bifStorage_->GetBifurcation(jt);
					size_t mustBeBifurcation = GetMustBeBifurcation(jt);
					assert(actualBifurcation == mustBeBifurcation);
				}	
			}
		}
	}
#endif

	DNASequence& IndexedSequence::Sequence()
	{
		return *sequence_;
	}

	const DNASequence& IndexedSequence::Sequence() const
	{
		return *sequence_;
	}

	BifurcationStorage& IndexedSequence::BifStorage()
	{
		return *bifStorage_;
	}

	const BifurcationStorage& IndexedSequence::BifStorage() const
	{
		return *bifStorage_;
	}
}