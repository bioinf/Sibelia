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
	const char DNASequence::DELETED_CHAR = -1;
	const Pos DNASequence::DELETED_POS = -1;
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
	
	DNASequence::DNASequence(const std::vector<FASTARecord> & record): sequence_(DNACharacter(DELETED_CHAR, DELETED_POS))
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

	DNASequence::DNASequence(const std::vector<FASTARecord> & record, const std::vector<std::vector<Pos> > & original):
		sequence_(DNACharacter(DELETED_CHAR, DELETED_POS))
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

	void DNASequence::Replace(StrandIterator source,
			size_t sourceDistance, 
			StrandIterator target,
			size_t targetDistance,
			Sequence::notify_func before,
			Sequence::notify_func after)
	{	
		size_t pos = 0;
		std::vector<size_t> oldPos;
		for(StrandIterator jt = target; pos < targetDistance; pos++, ++jt)
		{
			oldPos.push_back(jt.GetOriginalPosition());
		}

		SequencePosIterator begin = target.Base();
		SequencePosIterator end = AdvanceForward(target, targetDistance).Base();
		if(target.GetDirection() == positive)
		{
			begin = sequence_.erase(begin, end);
			begin = sequence_.insert(begin, source, AdvanceForward(source, sourceDistance), before, after);
			target = StrandIterator(begin, positive);
		}
		else
		{	
			std::swap(begin, end);
			source = AdvanceForward(source, sourceDistance).Invert();
			begin = sequence_.erase(begin, end);
			begin = sequence_.insert(begin, source, AdvanceForward(source, sourceDistance), before, after);
			target = StrandIterator(AdvanceForward(begin, sourceDistance), negative);
		}

		pos = 0;
		size_t record = 0;
		for(StrandIterator jt = target; pos < sourceDistance; pos++, ++jt)
		{
			size_t nowPos = record < oldPos.size() ? oldPos[record++] : oldPos.back();
			jt.it_->GetNaked()->pos = static_cast<Pos>(nowPos);
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
		posBegin_.clear();
		posEnd_.clear();		
		sequence_.erase(sequence_.begin(), sequence_.end());
	}
}