//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _COMMON_H_
#define _COMMON_H_

#include <ctime>
#include <cstdio>
#include <string>
#include <vector>
#include <cctype>
#include <cstdlib>
#include <cerrno>
#include <set>
#include <map>
#include <list>
#include <memory>
#include <cstring>
#include <fstream>
#include <sstream>
#include <numeric>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <functional>
#include <boost/bind.hpp>
#include <boost/cstdint.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>

#include <signal.h>
#include <divsufsort.h>

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
	TypeName(const TypeName&); \
	void operator=(const TypeName&)

extern const std::string VERSION;
extern const std::string DELIMITER;

namespace SyntenyFinder
{
	typedef uint64_t Pos;
	typedef uint64_t Size;
	const size_t MAX_INPUT_SIZE = size_t(1) << size_t(63);

	template<class Iterator1, class Iterator2>
		void CopyN(Iterator1 it, size_t count, Iterator2 out)
		{
			for(size_t i = 0; i < count; i++)
			{
				*out++ = *it++;
			}
		}

	template<class Iterator>
		Iterator AdvanceForward(Iterator it, size_t step)
		{
			std::advance(it, step);
			return it;
		}

	template<class Iterator>
		Iterator AdvanceBackward(Iterator it, size_t step)
		{
			for(size_t i = 0; i < step; i++)
			{
				--it;
			}

			return it;
		}

	template<class Iterator>
		Iterator AdvanceBackward(Iterator it, Iterator lowerBound, size_t step)
		{
			for(size_t i = 0; i < step && it != lowerBound; i++)
			{
				--it;
			}

			return it;
		}

	template<class Iterator>
		Iterator AdvanceForward(Iterator it, Iterator upperBound, size_t step)
		{
			for(size_t i = 0; i < step && it != upperBound; i++)
			{
				++it;
			}

			return it;
		}

	template<class Iterator, class F, class ReturnType>
		struct FancyIterator: public std::iterator<std::forward_iterator_tag, ReturnType>
		{
		public:
			FancyIterator& operator++()
			{
				++it;
				return *this;
			}

			FancyIterator operator++(int)
			{
				FancyIterator ret(*this);
				++(*this);
				return ret;
			}

			bool operator == (FancyIterator toCompare) const
			{
				return it == toCompare.it;
			}

			bool operator != (FancyIterator toCompare) const
			{
				return !(*this == toCompare);
			}

			ReturnType operator * () 
			{
				return f(*it);
			}

			FancyIterator() {}
			FancyIterator(Iterator it, F f): it(it), f(f) {}

		private:
			F f;
			Iterator it;
		};

	template<class Iterator, class F, class ReturnType>
		FancyIterator<Iterator, F, ReturnType> CFancyIterator(Iterator it, F f, ReturnType)
		{
			return FancyIterator<Iterator, F, ReturnType>(it, f);
		}

	typedef std::pair<size_t, size_t> IndexPair;
	template<class T, class F, class It>
		void GroupBy(std::vector<T> & store, F pred, It out)
		{
			sort(store.begin(), store.end(), pred);
			for(size_t now = 0; now < store.size(); )
			{
				size_t prev = now;
				for(; now < store.size() && !pred(store[prev], store[now]); now++);
				*out++ = std::make_pair(prev, now);
			}
		}

	template<class T>
		struct Counter
		{
		public:
			Counter(T state = 0): state_(state) {}
			T operator ()()
			{
				return state_++;
			}

		private:
			T state_;
		};

	inline bool RangesOverlap(std::pair<size_t, size_t> a, std::pair<size_t, size_t> b)
	{
		throw 1;
		return true;
	}
}

#endif
