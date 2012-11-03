#include "unrolledlist.h"

namespace SyntenyFinder
{
	struct T
	{
		int x;
		T() {}
		T(int x): x(x) {}
	};

	void Dumb()
	{
		typedef unrolled_list<T, 1000> UList;
		UList ulist(T(0));
		UList::reverse_iterator it;
		*it = T(3);
	}
}