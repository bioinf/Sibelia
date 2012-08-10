#ifndef _BIFURCATION_STORAGE_H_
#define _BIFURCATION_STORAGE_H_

#include "dnasequence.h"

namespace SyntenyBuilder
{
	typedef DNASequence::StrandIterator StrandIterator;

	class BifurcationStorage
	{
	public:
		static const size_t NO_BIFURCATION;

		void Clear();
		size_t GetMaxId() const;
		BifurcationStorage(): maxId_(0) {}
		void Dump(std::ostream & out) const;
		void ErasePoint(DNASequence::StrandIterator it);
		void AddPoint(DNASequence::StrandIterator it, size_t bifId);
		size_t CountBifurcations(size_t bifId) const;
		size_t GetBifurcation(DNASequence::StrandIterator it) const;		

		template<class Iterator>
			size_t ListPositions(size_t bifId, Iterator out) const
			{
				size_t ret = 0;
				for(size_t strand = 0; strand < 2; strand++)
				{
					std::pair<CBifMapIterator, CBifMapIterator> range = 
						bifurcationPos_[strand].equal_range(bifId);
					for(;range.first != range.second; ++range.first)
					{
						ret++;
						*out++ = range.first->second;
					}
				}

				return ret;
			}

	private:
		typedef boost::unordered_multimap<size_t, StrandIterator> BifurcationPos;
		typedef BifurcationPos::iterator BifMapIterator;
		typedef BifurcationPos::const_iterator CBifMapIterator;
		struct IteratorLess
		{
		public:
			bool operator () (BifMapIterator it1, BifMapIterator it2) const
			{
				return it1->second < it2->second;
			}
		};

		typedef std::set<BifurcationPos::iterator, IteratorLess> PosBifurcation;

		size_t maxId_;
		mutable BifurcationPos temp_;
		BifurcationPos bifurcationPos_[2];
		PosBifurcation posBifurcation_[2];
	};
}

#endif
