#ifndef _OUTPUT_GENERATOR_H_
#define _OUTPUT_GENERATOR_H_

#include "blockfinder.h"

namespace SyntenyFinder
{	
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

	class OutputGenerator
	{
	public:
		typedef std::vector<BlockInstance> BlockList;
		OutputGenerator(const ChrList & chrList, const BlockList & blockList): chrList_(chrList), blockList_(blockList) {}
		void GenerateReport(std::ostream & out) const;
		void ListBlocksIndices(std::ostream & out) const;
		void ListBlocksSequences(std::ostream & out) const;
		void ListChromosomesAsPermutations(std::ostream & out) const;
	private:
		ChrList chrList_;
		mutable BlockList blockList_;
		void ListChrs(std::ostream & out) const;
	};
}

#endif