#ifndef _COMMON_H_
#define _COMMON_H_

#include <ctime>
#include <cstdio>
#include <string>
#include <vector>
#include <cctype>
#include <cstdlib>
#include <stdint.h>
#include <set>
#include <map>
#include <list>
#include <memory>
#include <fstream>
#include <sstream>
#include <iterator>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <boost/ref.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <sparsehash/sparse_hash_map>
#include <sparsehash/sparse_hash_set>
#include <sparsehash/dense_hash_set>
#include <sparsehash/dense_hash_map>

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
	TypeName(const TypeName&); \
	void operator=(const TypeName&)

namespace SyntenyBuilder
{
	const std::string DELIMITER(80, '-');

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
}

#endif
