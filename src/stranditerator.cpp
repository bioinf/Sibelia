//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "dnasequence.h"

namespace SyntenyFinder
{/*
	namespace
	{
		typedef boost::uint64_t uint64;
		const uint64 ONE = 1;
	}

	const Size DNASequence::StrandIterator::INFO_BITS = 2;
	
	Size DNASequence::StrandIterator::PositionMask()
	{
		uint64 numBits = 8 * sizeof(Size) - 1;
		uint64 allBitsMask = (ONE << numBits) - ONE;
		uint64 infoBitsMask = (ONE << uint64(INFO_BITS)) - ONE;
		infoBitsMask <<= numBits - INFO_BITS;
		uint64 ret = allBitsMask ^ infoBitsMask;
 		return static_cast<Size>(ret);
	}

	bool DNASequence::StrandIterator::GetInfoBit(size_t bit) const
	{
		uint64 numBits = 8 * sizeof(Size) - 1;
		uint64 bitMask = ONE << (numBits - uint64(bit));
		return (it_.meta() & bitMask) != 0;
	}

	void DNASequence::StrandIterator::SetInfoBit(size_t bit, bool value) const
	{
		uint64 numBits = 8 * sizeof(Size) - 1;
		uint64 bitMask = ~(ONE << (numBits - uint64(bit)));
		it_.meta() = it_.meta() & bitMask;

		if(value)
		{
			bitMask = ~bitMask;
			Size res = static_cast<Size>(uint64(it_.meta()) | bitMask);
			it_.meta() = res;
		}
	}

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

	DNASequence::StrandIterator::StrandIterator()
	{
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
		return GetElementId() < toCompare.GetElementId();
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
		return it_.meta() & PositionMask();
	}

	void DNASequence::StrandIterator::SetOriginalPosition(size_t position) const
	{
		it_.meta() = static_cast<Size>(position) & PositionMask();
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

	FastaRecord::Direction DNASequence::StrandIterator::GetDirection() const
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
	}	*/
}