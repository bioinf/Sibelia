//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************


#ifndef _BIFURCATION_STORAGE_H_
#define _BIFURCATION_STORAGE_H_

#include "dnasequence.h"

namespace SyntenyFinder
{/*
	typedef DNASequence::StrandIterator StrandIterator;

	class BifurcationStorage
	{
	public:
		typedef Size BifurcationId;
		typedef DNASequence::SequencePosIterator PositiveIterator;
		typedef DNASequence::SequencePosIterator BaseIterator;
		typedef boost::container::slist<BaseIterator> IteratorList;
		typedef IteratorList::iterator IteratorPtr;
		typedef std::vector<IteratorList> ListVector;		
		typedef std::vector<IteratorPtr> IteratorVector;
		static const BifurcationId NO_BIFURCATION;

		void Clear();
		size_t GetMaxId() const;
		size_t TotalElements() const;
		size_t GetEmpty() const;
		BifurcationStorage(size_t maxId);
		void Cleanup();
		void Dump(const DNASequence & sequence, size_t k, std::ostream & out) const;
		void ErasePoint(DNASequence::StrandIterator it);
		void AddPoint(DNASequence::StrandIterator it, size_t bifId);
		size_t CountBifurcations(size_t bifId) const;
		size_t GetBifurcation(DNASequence::StrandIterator it) const;
		void NotifyBefore(StrandIterator begin, StrandIterator end);
		void NotifyAfter(StrandIterator begin, StrandIterator end);
		void FormDictionary(boost::unordered_map<std::string, size_t> & dict, size_t k) const;

		class IteratorProxy
		{
		public:
			IteratorProxy() {}
			IteratorProxy(IteratorPtr ptr, DNASequence::Direction direction): ptr_(ptr), direction_(direction)
			{
			}

			bool Valid() const;
			StrandIterator operator * () const;
		private:
			DNASequence::Direction direction_;
			IteratorPtr ptr_;			
		};

		template<class Iterator>
			size_t ListPositions(size_t inBifId, Iterator out)
			{
				size_t ret = 0;				
				for(size_t strand = 0; strand < 2; strand++)
				{
					for(IteratorList::iterator it = bifurcationPos_[strand][inBifId].begin(); it != bifurcationPos_[strand][inBifId].end(); ++it, ++ret)
					{
						*out++ = IteratorProxy(it, static_cast<DNASequence::Direction>(strand));
					}
				}

				return ret;
			}

	private:					
		static const size_t POSITIVE_BIT;
		static const size_t NEGATIVE_BIT;

		struct IteratorPtrHash
		{
		public:
			size_t operator () (IteratorPtr it) const
			{
				return reinterpret_cast<size_t>(&(**it));
			}
		};

		template<class T1, class T2>
			struct IteratorPtrEqual
			{
			public:
				bool operator () (const T1 & it1, const T2 & it2) const
				{
					return *it1 == *it2;
				}
			};

		struct BifurcationRecord
		{
			size_t pos;
			IteratorPtr ptrIt;
			BifurcationId bifId;
			BifurcationRecord() {}
			BifurcationRecord(size_t pos, IteratorPtr ptrIt, BifurcationId bifId): pos(pos), ptrIt(ptrIt), bifId(bifId) {}
		};

		BifurcationId ErasePointInternal(DNASequence::StrandIterator it, IteratorPtr & ret);
		
		typedef boost::unordered_set<IteratorPtr, IteratorPtrHash, IteratorPtrEqual<IteratorPtr, IteratorPtr> > IteratorMap;

		BifurcationId maxId_;
		ListVector bifurcationPos_[2];
		IteratorMap posBifurcation_[2];

		size_t nowInvalid_;		
		std::vector<std::vector<BifurcationRecord> > invalid_;
		std::vector<std::pair<IteratorList*, IteratorPtr> > toClear_;

		IteratorMap::iterator LookUp(StrandIterator it);
		IteratorMap::const_iterator LookUp(StrandIterator it) const;
	};*/
}

#endif
