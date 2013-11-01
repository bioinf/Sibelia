//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "indexedsequence.h"

namespace SyntenyFinder
{
	IndexedSequence::IndexedSequence(const std::vector<std::string> & record, std::vector<std::vector<Pos> > & originalPos, size_t k, const std::string & tempDir, bool clear, size_t model): k_(k)
	{
		Init(record, originalPos, k, tempDir, clear, model);
	}

	IndexedSequence::IndexedSequence(const std::vector<std::string> & record, size_t k, const std::string & tempDir, size_t model): k_(k)
	{
		std::vector<std::vector<Pos> > originalPos(record.size());
		for(size_t i = 0; i < originalPos.size(); i++)
		{
			originalPos[i].resize(record[i].size());
			std::generate(originalPos[i].begin(), originalPos[i].end(), Counter<Pos>());
		}

		Init(record, originalPos, k, tempDir, false, model);
	}

	void IndexedSequence::Init(std::vector<std::string> record, std::vector<std::vector<Pos> > & originalPos, size_t k, const std::string & tempDir, bool clear, size_t model)
	{
		size_t maxId;
		for(size_t i = 0; i < record.size(); i++)
		{
			for(size_t j = 0; j < record[i].size(); j++)
			{
				record[i][j] = IsDefiniteBase(record[i][j]) ? record[i][j] : DEFINITE_BASE[rand() % DEFINITE_BASE.size()];
			}
		}

		std::vector<std::vector<BifurcationInstance> > bifurcation(2);
		if(model == IndexedSequence::NO_MODEL)
		{
			if(tempDir.size() == 0)
			{
				maxId = EnumerateBifurcationsSArrayInRAM(record, bifurcation[0], bifurcation[1]);
			}
			else
			{
				maxId = EnumerateBifurcationsSArray(record, tempDir, bifurcation[0], bifurcation[1]);
			}
		}
		else
		{
			//maxId = 
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

	void IndexedSequence::ConstructChrIndex()
	{
		chrIndex_.clear();
		for(size_t chr = 0; chr < sequence_->ChrNumber(); chr++)
		{
			StrandIterator end = sequence_->PositiveEnd(chr);
			for(StrandIterator begin = sequence_->PositiveBegin(chr); begin != end; ++begin)
			{
				chrIndex_.push_back(std::make_pair(begin, chr));
			}
		}

		std::sort(chrIndex_.begin(), chrIndex_.end());
	}
	
	size_t IndexedSequence::GetChr(StrandIterator it) const
	{
		IteratorChrPair look(it, 0);
		std::vector<IteratorChrPair>::const_iterator jt = std::lower_bound(chrIndex_.begin(), chrIndex_.end(), look);
		return jt->second;
	}

	bool IndexedSequence::StrandIteratorPosGEqual(StrandIterator a, StrandIterator b)
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

	size_t IndexedSequence::StrandIteratorDistance(StrandIterator start, StrandIterator end)
	{
		size_t min = std::min(start.GetOriginalPosition(), end.GetOriginalPosition());
		size_t max = std::max(start.GetOriginalPosition(), end.GetOriginalPosition());
		return max - min;
	}
}