//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file COPYING for details.
//****************************************************************************


#ifndef _BIFURCATION_STORAGE_H_
#define _BIFURCATION_STORAGE_H_

#include "dnasequence.h"

namespace SyntenyFinder
{
	typedef DNASequence::StrandIterator StrandIterator;

	class BifurcationStorage
	{
	public:
		typedef Size BifurcationId;
		static const BifurcationId NO_BIFURCATION;

		void Clear();
		size_t GetMaxId() const;
		size_t TotalElements() const;
		BifurcationStorage();
		void Dump(const DNASequence & sequence, size_t k, std::ostream & out) const;
		void ErasePoint(DNASequence::StrandIterator it);
		void AddPoint(DNASequence::StrandIterator it, size_t bifId);
		size_t CountBifurcations(size_t bifId) const;
		size_t GetBifurcation(DNASequence::StrandIterator it) const;		

		template<class Iterator>
			size_t ListPositions(size_t inBifId, Iterator out) const
			{
				size_t ret = 0;
				BifurcationId bifId = static_cast<BifurcationId>(inBifId);
				for(size_t strand = 0; strand < 2; strand++)
				{
					std::pair<CBifMapIterator, CBifMapIterator> range = bifurcationPos_[strand].equal_range(bifId);
					for(;range.first != range.second; ++range.first)
					{
						ret++;
						*out++ = StrandIterator(range.first->second, static_cast<DNASequence::Direction>(strand));
					}
				}

				return ret;
			}

	private:
		typedef DNASequence::SequencePosIterator BaseIterator;
		typedef boost::unordered_multimap<BifurcationId, BaseIterator> BifurcationPos;
		typedef BifurcationPos::iterator BifMapIterator;
		typedef BifurcationPos::const_iterator CBifMapIterator;

		struct IteratorEqual
		{
		public:
			bool operator () (BifMapIterator it1, BifMapIterator it2) const
			{
				return it1->second == it2->second;
			}
		};

		struct IteratorHash
		{
		public:
			size_t operator () (BifMapIterator it) const
			{
				return reinterpret_cast<size_t>(&*it->second);
			}
		};

		typedef boost::unordered_set<BifurcationPos::iterator, IteratorHash, IteratorEqual> PosBifurcation;

		BifurcationId maxId_;
		mutable BifurcationPos temp_;
		std::vector<BifurcationPos> bifurcationPos_;
		std::vector<PosBifurcation> posBifurcation_;
	};
}

#endif
