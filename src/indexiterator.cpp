#include "indexiterator.h"

namespace SyntenyBuilder
{
	size_t IteratorUtility::AdvanceForward(const std::string * sequence_, size_t & pos_, char deleted_, size_t value)
	{
		if(pos_ < sequence_->size())
		{
			for(pos_ += value; pos_ < sequence_->size() && (*sequence_)[pos_] == deleted_; pos_++);
		}
		else if(pos_ == IndexIterator::NPOS && value == 1)
		{
			pos_ = 0;
		}

		return pos_;
	}

	size_t IteratorUtility::AdvanceBackward(const std::string * sequence_, size_t & pos_, char deleted_, size_t value)
	{
		if(pos_ > 0 && pos_ != IndexIterator::NPOS)
		{
			for(pos_ -= value; pos_ >= 0 && (*sequence_)[pos_] == deleted_; pos_--);
		}
		else if(value != 0)
		{
			pos_ = IndexIterator::NPOS;
		}

		return pos_;
	}
}