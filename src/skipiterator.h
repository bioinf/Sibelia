#include "common.h"

#ifndef _SKIP_ITERATOR_H
#define _SKIP_ITERATOR_H

namespace SyntenyBuilder
{
	template<class Iterator>
		class SkipIterator: public std::iterator<std::forward_iterator_tag, typename std::iterator_traits<Iterator>::value_type, void> 
		{
		public:
			typedef typename std::iterator_traits<Iterator>::value_type ValueType;
			SkipIterator(const Iterator & it, const Iterator & end, const ValueType & toSkip): it_(it), end_(end), toSkip_(toSkip)
			{
				if(it_ != end && *it_ == toSkip_)
				{
					++(*this);
				}
			}

			bool Valid() const
			{
				return it_ != end_;
			}

			Iterator CurrentPosition() const
			{
				return it_;
			}

			ValueType& operator *()
			{
				return *it_;
			}

			ValueType* operator ->()
			{
				return &*it_;
			}

			SkipIterator& operator ++ ()
			{
				while(++it_ != end_ && *it_ == toSkip_)
				{

				}

				return *this;
			}

			SkipIterator operator ++ (int)
			{
				SkipIterator ret(*this);
				while(++it_ != end_ && *it_ == toSkip_)
				{

				}

				return ret;
			}

			bool operator == (const SkipIterator & comp)
			{
				return it_ == comp.it_;
			}

			bool operator != (const SkipIterator & comp)
			{
				return it_ != comp.it_;
			}

		private:
			Iterator it_;
			Iterator end_;		
			ValueType toSkip_;
		};

	template<class Iterator>
		class SkipConstIterator: public std::iterator<std::forward_iterator_tag, typename std::iterator_traits<Iterator>::value_type, void> 
		{
		public:
			typedef typename std::iterator_traits<Iterator>::value_type ValueType;
			SkipConstIterator(const Iterator & it, const Iterator & end, const ValueType & toSkip): it_(it), end_(end), toSkip_(toSkip)
			{
				if(it_ != end && *it_ == toSkip_)
				{
					++(*this);
				}
			}

			bool Valid() const
			{
				return it_ != end_;
			}

			Iterator CurrentPosition() const
			{
				return it_;
			}

			const ValueType& operator *() const
			{
				return *it_;
			}

			const ValueType* operator ->() const
			{
				return &*it_;
			}

			SkipConstIterator& operator ++ ()
			{
				while(++it_ != end_ && *it_ == toSkip_)
				{

				}

				return *this;
			}

			SkipConstIterator operator ++ (int)
			{
				SkipConstIterator ret(*this);
				while(++it_ != end_ && *it_ == toSkip_)
				{

				}

				return ret;
			}

			bool operator == (const SkipConstIterator & comp)
			{
				return it_ == comp.it_;
			}

			bool operator != (const SkipConstIterator & comp)
			{
				return it_ != comp.it_;
			}

		private:
			Iterator it_;
			Iterator end_;		
			ValueType toSkip_;
		};
}
#endif