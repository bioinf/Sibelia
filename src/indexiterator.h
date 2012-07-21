#ifndef _INDEX_ITERATOR_H_
#define _INDEX_ITERATOR_H_

#include "common.h"

namespace SyntenyBuilder
{
	class IndexIterator;

	enum Direction
	{
		forward,
		reverse
	};

	struct IteratorUtility
	{
		static size_t AdvanceForward(const std::string * sequence_, size_t & pos_, char deleted_, size_t value = 1);
		static size_t AdvanceBackward(const std::string * sequence_, size_t & pos_, char deleted_, size_t value = 1);
	};

	class IndexIterator: public std::iterator<std::bidirectional_iterator_tag, char, size_t>
	{
	public:
		IndexIterator(): sequence_(0), pos_(NPOS) {}
		IndexIterator(std::string & sequence, size_t pos, char deleted, Direction dir = forward):
			sequence_(&sequence), pos_(pos), deleted_(deleted)
		{
			if(pos_ != sequence.size() && pos_ != NPOS)
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

		size_t GetPosition() const
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
			return sequence_ != 0 && pos_ != NPOS && pos_ < sequence_->size() && (*sequence_)[pos_] != deleted_;
		}

		static const size_t NPOS = -1;
	private:
		std::string * sequence_;
		size_t pos_;
		char deleted_;
		friend class IndexConstIterator;
	};

	class IndexConstIterator: public std::iterator<std::bidirectional_iterator_tag, char, size_t>
	{
	public:
		IndexConstIterator(): sequence_(0), pos_(NPOS) {}
		IndexConstIterator(const std::string & sequence, size_t pos, char deleted, Direction dir = forward):
			sequence_(&sequence), pos_(pos), deleted_(deleted)
		{
			if(pos_ != static_cast<size_t>(sequence.size()) && pos_ != NPOS)
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

		size_t GetPosition() const
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
			return sequence_ != 0 && pos_ != NPOS && pos_ < sequence_->size() && (*sequence_)[pos_] != deleted_;
		}

		static const size_t NPOS = -1;
	private:
		const std::string * sequence_;
		size_t pos_;
		char deleted_;
	};

	inline IndexIterator MakeRightEnd(std::string & str, char deleted)
	{
		return IndexIterator(str, static_cast<size_t>(str.size()), deleted);
	}

	inline IndexIterator MakeLeftEnd(std::string & str, char deleted)
	{
		return IndexIterator(str, IndexIterator::NPOS, deleted);
	}

	inline IndexConstIterator MakeLeftEnd(const std::string & str, char deleted)
	{
		return IndexConstIterator(str, IndexIterator::NPOS, deleted);
	}

	inline IndexConstIterator MakeRightEnd(const std::string & str, char deleted)
	{
		return IndexConstIterator(str, static_cast<size_t>(str.size()), deleted);
	}

	template<class Iterator1, class Iterator2>
		void CopyN(Iterator1 it, size_t count, Iterator2 out)
		{
			for(size_t i = 0; i < count; i++)
			{
				*out++ = *it++;
			}
		}

	template<class Iterator>
		Iterator AdvanceForward(Iterator it2, size_t step)
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
