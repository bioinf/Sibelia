#include "dnasequence.h"

namespace SyntenyBuilder
{
	inline std::string ConstructComplementarityTable()
	{
		std::string ret(1 << 8, ' ');
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
	
	const char DNASequence::EMPTY_CHARACTER = -1;
	const std::string DNASequence::alphabet("agct");	
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

	DNASequence::StrandIterator::StrandIterator(const StrandIterator & it): it_(it.it_->Clone())
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
		it_->MoveForward();
		return *this;
	}

	DNASequence::StrandIterator DNASequence::StrandIterator::operator++(int)
	{
		StrandIterator ret(*this);
		it_->MoveForward();
		return ret;
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

	char DNASequence::StrandIterator::operator * () const
	{
		return it_->Spell();
	}

	DNASequence::StrandIterator DNASequence::PositiveBegin() const
	{
		return StrandIterator(new ForwardIterator(sequence_.begin()));
	}

	DNASequence::StrandIterator DNASequence::PositiveEnd() const
	{
		return StrandIterator(new ForwardIterator(sequence_.end()));
	}

	DNASequence::StrandIterator DNASequence::NegativeBegin() const
	{
		return StrandIterator(new BackwardIterator(sequence_.rbegin()));
	}

	DNASequence::StrandIterator DNASequence::NegativeEnd() const
	{
		return StrandIterator(new BackwardIterator(sequence_.rend()));
	}

	DNASequence::DNASequence(const std::string & sequence)
	{
		for(size_t i = 0; i < sequence.size(); i++)
		{
			sequence_.push_back(DNACharacter(sequence[i], sequence[i]));
		}
	}

	size_t DNASequence::Size() const
	{
		return sequence_.size();
	}
}