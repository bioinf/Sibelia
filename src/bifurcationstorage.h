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
		void Dump(std::ostream & out) const;
		void ErasePoint(DNASequence::StrandIterator it);
		void AddPoint(DNASequence::StrandIterator it, size_t bifId);
		size_t GetBifurcation(DNASequence::StrandIterator it) const;

		template<class Iterator>
			size_t ListPositions(size_t bifId, Iterator out, const DNASequence & seq) const
			{
				size_t ret = 0;
				typedef boost::function<StrandIterator (size_t)> Transformer;
				Transformer trans[2] = 
				{
					boost::bind(&DNASequence::PositiveByIndex, boost::cref(seq), _1),
					boost::bind(&DNASequence::NegativeByIndex, boost::cref(seq), _1)
				};

				for(size_t strand = 0; strand < 2; strand++)
				{
					std::pair<CBifMapIterator, CBifMapIterator> range = 
						bifurcationPos[strand].equal_range(bifId);
					for(CBifMapIterator it = range.first; it != range.second; ++it, ++ret)
					{
						assert(GetBifurcation(trans[strand](it->second)) == bifId);
						*out++ = trans[strand](it->second);
					}
				}

				return ret;
			}

	private:
		typedef boost::unordered_multimap<size_t, size_t> BifurcationPos;
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

		mutable BifurcationPos temp;
		BifurcationPos bifurcationPos[2];
		PosBifurcation posBifurcation[2];
	};
	
}

#endif
