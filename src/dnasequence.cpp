//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "dnasequence.h"

namespace SyntenyFinder
{
	inline std::string ConstructComplementarityTable()
	{
		std::string ret(1 << (sizeof(char) * 8), ' ');
		for(size_t i = 0; i < ret.size(); i++)
		{
			ret[i] = static_cast<char>(i);
		}

		ret['a'] = 't';
		ret['t'] = 'a';
		ret['g'] = 'c';
		ret['c'] = 'g';
		ret['A'] = 'T';
		ret['T'] = 'A';
		ret['G'] = 'C';
		ret['C'] = 'G';
		return ret;
	}
	
	const char DNASequence::UNKNOWN_BASE = 'n';
	const char DNASequence::SEPARATION_CHAR = '$';
	const std::string DNASequence::alphabet("agctn$");
	const std::string DNASequence::complementary_(ConstructComplementarityTable());

	char DNASequence::Translate(char ch)
	{
		return complementary_[ch];
	}

	DNASequence::StrandIterator::StrandIterator()
	{
		
	}

	DNASequence::StrandIterator::StrandIterator(GenericIterator * it): it_(it)
	{
		
	}

	DNASequence::StrandIterator::StrandIterator(SequencePosIterator it, Direction direction)		
	{
		GenericIterator * genericIt = 0;
		if(direction == positive)
		{
			genericIt = new ForwardIterator(it);
		}
		else
		{
			genericIt = new BackwardIterator(SequenceNegIterator(it));
		}

		it_ = std::auto_ptr<GenericIterator>(genericIt);
	}

	DNASequence::StrandIterator::StrandIterator(const StrandIterator & toCopy): it_(toCopy.it_->Clone())
	{
	}

	void DNASequence::StrandIterator::Swap(StrandIterator & toSwap)
	{
		GenericIterator * me = it_.release();
		GenericIterator * he = toSwap.it_.release();
		it_ = std::auto_ptr<GenericIterator>(he);
		toSwap.it_ = std::auto_ptr<GenericIterator>(me);
	}

	DNASequence::StrandIterator& DNASequence::StrandIterator::operator = (const StrandIterator & toCopy)
	{
		StrandIterator temp(toCopy);
		Swap(temp);
		return *this;
	}

	size_t DNASequence::StrandIterator::GetElementId() const
	{
		return reinterpret_cast<size_t>(it_->GetNaked());
	}

	bool DNASequence::StrandIterator::AtValidPosition() const
	{
		return **this != DNASequence::SEPARATION_CHAR;
	}

	bool DNASequence::StrandIterator::operator < (const StrandIterator & toCompare) const
	{
		if(GetDirection() == toCompare.GetDirection())
		{
			return it_->GetNaked() < toCompare.it_->GetNaked();
		}

		return GetDirection() == positive;
	}

	bool DNASequence::StrandIterator::operator == (const StrandIterator & toCompare) const
	{
		return it_->Equal(*toCompare.it_);
	}

	bool DNASequence::StrandIterator::operator != (const StrandIterator & toCompare) const
	{
		return !(*this == toCompare);
	}

	DNASequence::StrandIterator& DNASequence::StrandIterator::operator++()
	{
		assert(AtValidPosition());
		it_->MoveForward();
		return *this;
	}

	DNASequence::StrandIterator DNASequence::StrandIterator::operator++(int)
	{
		assert(AtValidPosition());
		StrandIterator ret(*this);
		it_->MoveForward();
		return ret;
	}

	size_t DNASequence::StrandIterator::GetOriginalPosition() const
	{
		return it_->GetNaked()->pos;
	}

	DNASequence::StrandIterator& DNASequence::StrandIterator::operator--()
	{
		it_->MoveBackward();
		return *this;
	}

	DNASequence::StrandIterator DNASequence::StrandIterator::operator--(int)
	{
		StrandIterator ret(*this);
		it_->MoveBackward();
		return ret;
	}

	DNASequence::Direction DNASequence::StrandIterator::GetDirection() const
	{
		return it_->GetDirection();
	}

	DNASequence::StrandIterator DNASequence::StrandIterator::Invert() const
	{
		return StrandIterator(it_->Invert());
	}

	char DNASequence::StrandIterator::operator * () const
	{
		return it_->Spell();
	}

	const DNASequence::DNACharacter* DNASequence::StrandIterator::GetNaked() const
	{
		return it_->GetNaked();
	}

	DNASequence::StrandIterator DNASequence::PositiveBegin(size_t chr) const
	{
		Sequence & ref = const_cast<Sequence&>(sequence_);
		return StrandIterator(new ForwardIterator(posBegin_[chr]));
	}

	DNASequence::StrandIterator DNASequence::PositiveEnd(size_t chr) const
	{
		return StrandIterator(new ForwardIterator(posEnd_[chr]));
	}

	DNASequence::StrandIterator DNASequence::NegativeBegin(size_t chr) const
	{
		return StrandIterator(new BackwardIterator(SequenceNegIterator(posEnd_[chr])));
	}

	DNASequence::StrandIterator DNASequence::NegativeEnd(size_t chr) const
	{
		return StrandIterator(new BackwardIterator(SequenceNegIterator(posBegin_[chr])));
	}

	DNASequence::StrandIterator DNASequence::Begin(Direction direction, size_t chr) const
	{
		return direction == positive ? PositiveBegin(chr) : NegativeBegin(chr);
	}
	DNASequence::StrandIterator DNASequence::End(Direction direction, size_t chr) const
	{
		return direction == positive ? PositiveEnd(chr) : NegativeEnd(chr);
	}

	DNASequence::SequencePosIterator DNASequence::StrandIterator::Base() const
	{
		return it_->Base();
	}

	DNASequence::DNASequence(const std::vector<FASTARecord> & record) 
	{
		sequence_.push_back(DNACharacter(SEPARATION_CHAR, -1));
		for(size_t chr = 0; chr < record.size(); chr++)
		{
			SequencePosIterator chrPosBegin = --sequence_.end();
			for(size_t pos = 0; pos < record[chr].sequence.size(); pos++)
			{
				sequence_.push_back(DNACharacter(record[chr].sequence[pos], Pos(pos)));
			}

			sequence_.push_back(DNACharacter(SEPARATION_CHAR, Pos(record[chr].sequence.size())));
			posBegin_.push_back(++chrPosBegin);
			posEnd_.push_back(--sequence_.end());
		}
	}

	DNASequence::DNASequence(const std::vector<FASTARecord> & record, const std::vector<std::vector<Pos> > & original) 
	{
		sequence_.push_back(DNACharacter(SEPARATION_CHAR, -1));
		for(size_t chr = 0; chr < record.size(); chr++)
		{
			SequencePosIterator chrPosBegin = --sequence_.end();
			for(size_t pos = 0; pos < record[chr].sequence.size(); pos++)
			{
				sequence_.push_back(DNACharacter(record[chr].sequence[pos], original[chr][pos]));
			}

			sequence_.push_back(DNACharacter(SEPARATION_CHAR, Pos(record[chr].sequence.size())));
			posBegin_.push_back(++chrPosBegin);
			posEnd_.push_back(--sequence_.end());
		}
	}

	size_t DNASequence::TotalSize() const
	{
		return sequence_.size();
	}

	size_t DNASequence::ChrNumber() const
	{
		return posBegin_.size();
	}

	char DNASequence::StrandIterator::TranslateChar(char ch) const
	{
		return it_->TranslateChar(ch);
	}

	void DNASequence::EraseN(StrandIterator now, size_t count)
	{
		SequencePosIterator it = now.Base();
		if(now.GetDirection() == negative)
		{
			it = AdvanceForward(now, count).Base();
		}

		for(size_t i = 0; i < count; i++)
		{
			it = sequence_.erase(it);
		}
	}

	void DNASequence::CopyN(StrandIterator source, size_t count, StrandIterator target)
	{
		for(size_t i = 0; i < count; i++)
		{
			DNACharacter * ptr = const_cast<DNACharacter*>(target.GetNaked());
			ptr->actual = target.TranslateChar(*source);
			++target;
			++source;
		}
	}

	void DNASequence::Replace(StrandIterator source, size_t sourceDistance, 
			StrandIterator target, size_t targetDistance,
			const boost::function<void (const StrandIterator&)> & alarmBefore,
			const boost::function<void (const StrandIterator&)> & alarmAfter)
	{
		size_t pos = target.GetNaked()->pos;
		StrandIterator saveTarget = target;
		SequencePosIterator it = target.Base();
		if(target.GetDirection() == negative)
		{
			it = AdvanceForward(target, targetDistance).Base();
			source = AdvanceForward(source, sourceDistance - 1);
		}

		for(size_t i = 0; i < sourceDistance; i++)
		{
			if(i == 0 && it != sequence_.begin())
			{
				alarmBefore(StrandIterator(new BackwardIterator(SequenceNegIterator(it))));
			}

			sequence_.insert(it, DNACharacter(target.TranslateChar(*source), Pos(pos)));
			if(i == 0 && it != sequence_.begin())
			{
				--it;
				alarmAfter(StrandIterator(new BackwardIterator(SequenceNegIterator(it++))));
			}

			if(target.GetDirection() == positive)
			{
				++source;
			}
			else
			{
				--source;
			}			
		}

		for(size_t i = 0; i < targetDistance; i++)
		{						
			it = sequence_.erase(it);			
		}
	}

	std::pair<size_t, size_t> DNASequence::SpellOriginal(StrandIterator it1, StrandIterator it2) const
	{
		--it2;
		size_t start = std::min(it1.GetOriginalPosition(), it2.GetOriginalPosition());
		size_t end = std::max(it1.GetOriginalPosition(), it2.GetOriginalPosition());
		return std::make_pair(start, end + 1);
	}

	void DNASequence::Clear()
	{
		sequence_.clear();
		posBegin_.clear();
		posEnd_.clear();
	}
}