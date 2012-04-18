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
#include <google/sparse_hash_map>
#include <google/sparse_hash_set>
#include <google/dense_hash_set>
#include <google/dense_hash_map>

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
	TypeName(const TypeName&); \
	void operator=(const TypeName&)

#define _DEBUG_

namespace SyntenyBuilder
{
}

#endif
