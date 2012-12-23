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

	DNASequence::StrandIterator DNASequence::PositiveBegin(size_t chr) const
	{
		return StrandIterator(posBegin_[chr], positive);
	}

	DNASequence::StrandIterator DNASequence::PositiveEnd(size_t chr) const
	{
		return StrandIterator(posEnd_[chr], positive);
	}

	DNASequence::StrandIterator DNASequence::NegativeBegin(size_t chr) const
	{
		return PositiveEnd(chr).Invert();
	}

	DNASequence::StrandIterator DNASequence::NegativeEnd(size_t chr) const
	{
		return PositiveBegin(chr).Invert();
	}

	DNASequence::StrandIterator DNASequence::Begin(Direction direction, size_t chr) const
	{
		return direction == positive ? PositiveBegin(chr) : NegativeBegin(chr);
	}
	DNASequence::StrandIterator DNASequence::End(Direction direction, size_t chr) const
	{
		return direction == positive ? PositiveEnd(chr) : NegativeEnd(chr);
	}

	DNASequence::DNASequence(const std::vector<std::string> & record, std::vector<std::vector<Pos> > & original, bool clear):
		sequence_(DNACharacter(DELETED_CHAR))
	{
		sequence_.push_back(DNACharacter(SEPARATION_CHAR));
		for(size_t chr = 0; chr < record.size(); chr++)
		{
			SequencePosIterator chrPosBegin = --sequence_.end();
			for(size_t pos = 0; pos < record[chr].size(); pos++)
			{
				sequence_.push_back(DNACharacter(record[chr][pos]));
				StrandIterator(--sequence_.end(), positive).SetOriginalPosition(original[chr][pos]);
			}

			if(clear)
			{
				original[chr].clear();
				std::vector<Pos> temp;
				original[chr].swap(temp);
			}

			sequence_.push_back(DNACharacter(SEPARATION_CHAR));
			StrandIterator(--sequence_.end(), positive).SetOriginalPosition(record[chr].size());			
			posBegin_.push_back(++chrPosBegin);
			posEnd_.push_back(--sequence_.end());
		}

		std::for_each(posBegin_.begin(), posBegin_.end(), boost::bind(&DNASequence::SubscribeIterator, boost::ref(*this), _1));
		std::for_each(posEnd_.begin(), posEnd_.end(), boost::bind(&DNASequence::SubscribeIterator, boost::ref(*this), _1));
	}
	
	size_t DNASequence::TotalSize() const
	{
		return sequence_.size();
	}

	size_t DNASequence::ChrNumber() const
	{
		return posBegin_.size();
	}
	
	void DNASequence::SubscribeIterator(SequencePosIterator & it)
	{
		iteratorStore_.insert(&it);
	}

	void DNASequence::UnsubscribeIterator(SequencePosIterator & it)
	{		
		for(IteratorRange range = iteratorStore_.equal_range(&it); range.first != range.second; ++range.first)
		{
			if(*range.first == &it)
			{
				range.first = iteratorStore_.erase(range.first);
			}
		}
	}

	void DNASequence::NotifyBefore(SequencePosIterator begin, SequencePosIterator end, NotifyFunction before)
	{	
		if(before)
		{
			SequencePosIterator nowBegin = begin;
			while(nowBegin != end)
			{
				SequencePosIterator nowEnd = nowBegin;
				for(; nowEnd != end && nowEnd->actual != SEPARATION_CHAR; ++nowEnd);
				StrandIterator pbegin(nowBegin, positive);
				StrandIterator pend(nowEnd, positive);
				before(pbegin, pend);
				before(pend.Invert(), pbegin.Invert());
				nowBegin = nowEnd != end ? ++nowEnd : nowEnd;
			}			
		}

		for(; begin != end; ++begin)
		{			
			toReplace_.push_back(iteratorStore_.find(&begin));
		}
	}

	void DNASequence::NotifyAfter(SequencePosIterator begin, SequencePosIterator end, NotifyFunction after)
	{	
		if(after)
		{
			SequencePosIterator nowBegin = begin;
			while(nowBegin != end)
			{
				SequencePosIterator nowEnd = nowBegin;
				for(; nowEnd != end && nowEnd->actual != SEPARATION_CHAR; ++nowEnd);
				StrandIterator pbegin(nowBegin, positive);
				StrandIterator pend(nowEnd, positive);
				after(pbegin, pend);
				after(pend.Invert(), pbegin.Invert());
				nowBegin = nowEnd != end ? ++nowEnd : nowEnd;
			}			
		}

		size_t pos = 0;
		std::vector<boost::reference_wrapper<SequencePosIterator> > resubscribe;
		for(; begin != end; ++begin, ++pos)
		{
			if(toReplace_[pos] != iteratorStore_.end())
			{
				IteratorPlace it = toReplace_[pos];
				**it = begin;
				resubscribe.push_back(boost::ref(**it));
				it = iteratorStore_.erase(it);
			}
		}

		toReplace_.clear();
		assert(posEnd_.back() == --sequence_.end());
		std::for_each(resubscribe.begin(), resubscribe.end(), boost::bind(&DNASequence::SubscribeIterator, boost::ref(*this), _1));
	}

	DNASequence::SequencePosIterator DNASequence::ReplaceDirect(StrandIterator source,
			size_t sourceDistance, 
			SequencePosIterator target,
			size_t targetDistance,
			NotifyFunction before,
			NotifyFunction after)
	{	
		StrandIterator save = StrandIterator(target, positive);
		size_t firstPos = save.GetOriginalPosition();
		size_t lastPos = AdvanceForward(save, targetDistance).GetOriginalPosition();
		for(size_t i = 0; i < std::min(sourceDistance, targetDistance); i++)
		{
			*target = *source;
			++target;
			++source;
		}

		if(sourceDistance < targetDistance)
		{
			SequencePosIterator targetEnd = AdvanceForward(target, targetDistance - sourceDistance);
			target = sequence_.erase(target, targetEnd);			
		}
		else if(sourceDistance != targetDistance)
		{			
			Sequence::notify_func seqBefore = boost::bind(&DNASequence::NotifyBefore, boost::ref(*this), _1, _2, before);
			Sequence::notify_func seqAfter = boost::bind(&DNASequence::NotifyAfter, boost::ref(*this), _1, _2, after);
			StrandIterator sourceEnd = AdvanceForward(source, sourceDistance - targetDistance);
			std::string buf(source, sourceEnd);			
			target = sequence_.insert(target, buf.begin(), buf.end(), seqBefore, seqAfter);
			target = AdvanceForward(target, sourceDistance - targetDistance);
		}

		double acc = static_cast<double>(firstPos);
		double ssize = double(targetDistance) / sourceDistance;
		for(size_t step = 0; step < sourceDistance; step++, ++save, acc += ssize)
		{
			size_t pos = std::min(lastPos, size_t(acc));
			save.SetOriginalPosition(pos);
		}

		return target;
	}

	void DNASequence::Replace(StrandIterator source,
			size_t sourceDistance, 
			StrandIterator target,
			size_t targetDistance,
			NotifyFunction before,
			NotifyFunction after)
	{	
		if(target.GetDirection() == positive)
		{
			SequencePosIterator begin = target.Base();
			begin = ReplaceDirect(source, sourceDistance, begin, targetDistance, before, after);
			target = AdvanceBackward(StrandIterator(begin, positive), sourceDistance);
		}
		else
		{
			source = AdvanceForward(source, sourceDistance).Invert();
			SequencePosIterator begin = AdvanceForward(target, targetDistance).Invert().Base();
			begin = ReplaceDirect(source, sourceDistance, begin, targetDistance, before, after);
			target = StrandIterator(--begin, negative);
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

	size_t DNASequence::GlobalIndex(StrandIterator it) const
	{
		return std::distance(++(const_cast<DNASequence*>(this)->sequence_.begin()), it.Base());
	}
}