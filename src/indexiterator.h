#ifndef _INDEX_ITERATOR_H_
#define _INDEX_ITERATOR_H_

#include "common.h"

namespace SyntenyBuilder
{
	enum Direction
	{
		forward,
		reverse
	};

	struct IteratorUtility
	{
		static int AdvanceForward(const std::string * sequence_, int & pos_, char deleted_, int value = 1)
		{
			if(pos_ < static_cast<int>(sequence_->size()))
			{
				for(pos_ += value; pos_ < static_cast<int>(sequence_->size()) && (*sequence_)[pos_] == deleted_; pos_++);
			}

			return pos_;
		}

		static int AdvanceBackward(const std::string * sequence_, int & pos_, char deleted_, int value = 1)
		{
			if(pos_ >= 0)
			{
				for(pos_ -= value; pos_ >= 0 && (*sequence_)[pos_] == deleted_; pos_--);
			}

			return pos_;
		}
	};

	class IndexIterator;
	class IndexConstIterator;

	class IndexIterator: public std::iterator<std::bidirectional_iterator_tag, char, void>
	{
	public:
		IndexIterator(): sequence_(0), pos_(std::string::npos) {}
		IndexIterator(std::string & sequence, int pos, char deleted, Direction dir = forward):
			sequence_(&sequence), pos_(pos), deleted_(deleted)
		{
			if(pos_ != static_cast<int>(sequence.size()) && pos_ != std::string::npos)
			{
				if(dir == forward)
				{
					IteratorUtility::AdvanceForward(sequence_, pos_, deleted_, 0);
				}
				else
				{
					IteratorUtility::AdvanceBackward(sequence_, pos_, deleted_, 0);
				}
			}
		}

		char& operator * () const
		{
			return (*sequence_)[pos_];
		}

		IndexIterator& operator ++ ()
		{
			IteratorUtility::AdvanceForward(sequence_, pos_, deleted_);
			return *this;
		}

		IndexIterator operator ++ (int)
		{
			IndexIterator ret(*sequence_, pos_, deleted_);
			IteratorUtility::AdvanceForward(sequence_, pos_, deleted_);
			return ret;
		}

		IndexIterator& operator -- ()
		{
			IteratorUtility::AdvanceBackward(sequence_, pos_, deleted_);
			return *this;
		}

		IndexIterator operator -- (int)
		{
			IndexIterator ret(*sequence_, pos_, deleted_);
			IteratorUtility::AdvanceBackward(sequence_, pos_, deleted_);
			return ret;
		}

		int GetPosition() const
		{
			return pos_;
		}

		bool operator == (IndexIterator it) const
		{
			return it.sequence_ == sequence_ && pos_ == it.pos_ && deleted_ == it.deleted_;
		}

		bool operator != (IndexIterator it) const
		{
			return !(*this == it);
		}

		bool Valid() const
		{
			return sequence_ != 0 && pos_ >= 0 && pos_ < static_cast<int>(sequence_->size());
		}

	private:
		std::string * sequence_;
		int pos_;
		char deleted_;
		friend class IndexConstIterator;
	};

	class IndexConstIterator: public std::iterator<std::bidirectional_iterator_tag, char, void>
	{
	public:
		IndexConstIterator(): sequence_(0), pos_(std::string::npos) {}
		IndexConstIterator(const std::string & sequence, int pos, char deleted, Direction dir = forward):
			sequence_(&sequence), pos_(pos), deleted_(deleted)
		{
			if(pos_ != static_cast<int>(sequence.size()) && pos_ != std::string::npos)
			{
				if(dir == forward)
				{
					IteratorUtility::AdvanceForward(sequence_, pos_, deleted_, 0);
				}
				else
				{
					IteratorUtility::AdvanceBackward(sequence_, pos_, deleted_, 0);
				}
			}
		}

		IndexConstIterator(IndexIterator other): sequence_(other.sequence_),
			pos_(other.pos_), deleted_(other.deleted_)
		{
		}

		char operator * () const
		{
			return (*sequence_)[pos_];
		}

		IndexConstIterator& operator ++ ()
		{
			IteratorUtility::AdvanceForward(sequence_, pos_, deleted_);
			return *this;
		}

		IndexConstIterator operator ++ (int)
		{
			IndexConstIterator ret(*sequence_, pos_, deleted_);
			IteratorUtility::AdvanceForward(sequence_, pos_, deleted_);
			return ret;
		}

		IndexConstIterator& operator -- ()
		{
			IteratorUtility::AdvanceBackward(sequence_, pos_, deleted_);
			return *this;
		}

		IndexConstIterator operator -- (int)
		{
			IndexConstIterator ret(*sequence_, pos_, deleted_);
			IteratorUtility::AdvanceBackward(sequence_, pos_, deleted_);
			return ret;
		}

		int GetPosition() const
		{
			return pos_;
		}

		bool operator == (IndexConstIterator it) const
		{
			return it.sequence_ == sequence_ && pos_ == it.pos_ && deleted_ == it.deleted_;
		}

		bool operator != (IndexConstIterator it) const
		{
			return !(*this == it);
		}

		bool Valid() const
		{
			return sequence_ != 0 && pos_ >= 0 && pos_ < static_cast<int>(sequence_->size());
		}

	private:
		const std::string * sequence_;
		int pos_;
		char deleted_;
	};
	

	inline IndexIterator MakeRightEnd(std::string & str, char deleted)
	{
		return IndexIterator(str, str.size(), deleted);
	}

	inline IndexIterator MakeLeftEnd(std::string & str, char deleted)
	{
		return IndexIterator(str, std::string::npos, deleted);
	}

	inline IndexConstIterator MakeLeftEnd(const std::string & str, char deleted)
	{
		return IndexConstIterator(str, std::string::npos, deleted);
	}

	inline IndexConstIterator MakeRightEnd(const std::string & str, char deleted)
	{
		return IndexConstIterator(str, str.size(), deleted);
	}

	template<class Iterator>
		Iterator Advance(Iterator it2, size_t step)
		{
			std::advance(it2, step);
			return it2;
		}

	template<class Iterator>
		Iterator AdvanceBackward(Iterator it2, size_t step)
		{
			for(size_t i = 0; i < step; i++)
			{
				--it2;
			}

			return it2;
		}	
}

#endif
