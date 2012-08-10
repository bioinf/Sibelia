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
	
	const size_t DNASequence::NO_POS = -1;
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

	DNASequence::SequencePosIterator DNASequence::StrandIterator::Base() const
	{
		return it_->Base();
	}

	DNASequence::DNASequence(const std::string & sequence): original_(sequence)
	{
		for(size_t i = 0; i < sequence.size(); i++)
		{
			sequence_.push_back(DNACharacter(sequence[i], i));
		}
	}

	size_t DNASequence::Size() const
	{
		return sequence_.size();
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
			StrandIterator target, size_t targetDistance)
	{
		SequencePosIterator it = target.Base();
		if(target.GetDirection() == negative)
		{
			it = AdvanceForward(target, targetDistance).Base();
		}

		for(size_t i = 0; i < sourceDistance; i++)
		{
			sequence_.insert(it, DNACharacter(target.TranslateChar(*source), source.GetNaked()->pos));
			++source;
			if(i != targetDistance - 1)
			{
				++target;
			}
		}

		for(size_t i = 0; i < targetDistance; i++)
		{
			it = sequence_.erase(it);
		}
	}
}