//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
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
		typedef DNASequence::SequencePosIterator PositiveIterator;
		typedef DNASequence::SequenceNegIterator NegativeIterator;
		typedef DNASequence::SequencePosIterator BaseIterator;
		typedef boost::shared_ptr<BaseIterator> IteratorPtr;
		static const BifurcationId NO_BIFURCATION;

		void Clear();
		size_t GetMaxId() const;
		size_t TotalElements() const;
		BifurcationStorage(size_t maxId);
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
			IteratorProxy(IteratorPtr ptr, DNASequence::Direction direction): ptr_(ptr), direction_(direction) {}
			bool Valid() const;
			StrandIterator GetIterator() const;
		private:
			IteratorPtr ptr_;
			DNASequence::Direction direction_;
		};

		template<class Iterator>
			size_t ListPositions(size_t inBifId, Iterator out) const
			{
				size_t ret = 0;
				BifurcationId bifId = static_cast<BifurcationId>(inBifId);
				for(size_t strand = 0; strand < 2; strand++)
				{
					std::transform(bifurcationPos_[strand][bifId].begin(), bifurcationPos_[strand][bifId].end(), out,
						boost::bind(boost::value_factory<IteratorProxy>(), _1, static_cast<DNASequence::Direction>(pos)));
				}

				return ret;
			}

	private:
		typedef std::pair<size_t, BifurcationId> BifurcationRecord;		
		typedef std::vector<IteratorPtr> IteratorVector;
		typedef std::vector<IteratorVector> BifurcationStore;

		struct IteratorPtrHash
		{
		public:
			size_t operator () (IteratorPtr it) const
			{
				return reinterpret_cast<size_t>(&(**it));
			}
		};

		typedef boost::unordered_map<IteratorPtr, BifurcationId, IteratorPtrHash> IteratorMap;

		BifurcationId maxId_;
		BifurcationStore bifurcationPos_[2];
		IteratorMap posBifurcation_[2];

		size_t nowInvalid_;
		std::vector<std::vector<BifurcationRecord> > invalid_;
	};
}

#endif
