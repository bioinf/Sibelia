//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "dnasequence.h"

namespace SyntenyFinder
{
	
	DNASequence::StrandIterator::StrandIterator(SequencePosIterator it, Direction direction):
		it_(it), direction_(direction)
	{
	}

	DNASequence::SequencePosIterator DNASequence::StrandIterator::Base() const
	{
		return it_;
	}	
	
	char DNASequence::StrandIterator::TranslateChar(char ch) const
	{
		if(direction_ == positive)
		{
			return ch;
		}
		
		return complementary_[ch];
	}

	DNASequence::StrandIterator::StrandIterator(const StrandIterator & toCopy):
		it_(toCopy.it_), direction_(toCopy.direction_)
	{
	}

	void DNASequence::StrandIterator::Swap(StrandIterator & toSwap)
	{
		std::swap(it_, toSwap.it_);
		std::swap(direction_, toSwap.direction_);
	}

	DNASequence::StrandIterator& DNASequence::StrandIterator::operator = (const StrandIterator & toCopy)
	{
		StrandIterator temp(toCopy);
		Swap(temp);
		return *this;
	}

	size_t DNASequence::StrandIterator::GetElementId() const
	{
		return reinterpret_cast<size_t>(&*it_);
	}

	bool DNASequence::StrandIterator::AtValidPosition() const
	{
		return **this != DNASequence::SEPARATION_CHAR;
	}

	bool DNASequence::StrandIterator::operator < (const StrandIterator & toCompare) const
	{
		if(GetDirection() == toCompare.GetDirection())
		{
			return GetElementId() < toCompare.GetElementId();
		}

		return GetDirection() == positive;
	}

	bool DNASequence::StrandIterator::operator == (const StrandIterator & toCompare) const
	{
		return (it_ == toCompare.it_) && (direction_ == toCompare.direction_);
	}

	bool DNASequence::StrandIterator::operator != (const StrandIterator & toCompare) const
	{
		return !(*this == toCompare);
	}

	DNASequence::StrandIterator& DNASequence::StrandIterator::operator++()
	{
		assert(AtValidPosition());
		if(direction_ == positive)
		{
			++it_;
		}
		else
		{
			--it_;
		}

		return *this;
	}

	DNASequence::StrandIterator DNASequence::StrandIterator::operator++(int)
	{
		assert(AtValidPosition());
		StrandIterator ret(*this);
		if(direction_ == positive)
		{
			++it_;
		}
		else
		{
			--it_;
		}

		return ret;
	}

	size_t DNASequence::StrandIterator::GetOriginalPosition() const
	{
		return it_->pos;
	}

	DNASequence::StrandIterator& DNASequence::StrandIterator::operator--()
	{
		if(direction_ == negative)
		{
			++it_;
		}
		else
		{
			--it_;
		}

		return *this;
	}

	DNASequence::StrandIterator DNASequence::StrandIterator::operator--(int)
	{
		StrandIterator ret(*this);
		if(direction_ == negative)
		{
			++it_;
		}
		else
		{
			--it_;
		}

		return ret;
	}

	DNASequence::Direction DNASequence::StrandIterator::GetDirection() const
	{
		return direction_;
	}

	DNASequence::StrandIterator DNASequence::StrandIterator::Invert() const
	{
		if(direction_ == positive)
		{
			return StrandIterator(AdvanceBackward(it_, 1), negative);
		}

		return StrandIterator(AdvanceForward(it_, 1), positive);
	}

	char DNASequence::StrandIterator::operator * () const
	{
		if(direction_ == positive)
		{
			return it_->actual;
		}

		return complementary_[it_->actual];
	}	
}