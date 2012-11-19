//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _DNA_SEQUENCE_H_
#define _DNA_SEQUENCE_H_

#include "fasta.h"
#include "common.h"
#include "unrolledlist.h"

#pragma warning(disable:4355)

namespace SyntenyFinder
{
	class DNASequence
	{
	public:
		enum Direction
		{
			positive,
			negative
		};

		struct DNACharacter
		{
			char actual;
			DNACharacter() {}
			DNACharacter(char actual): actual(actual) {}
			bool operator == (const DNACharacter & a) const
			{
				return actual == a.actual;
			}

			bool operator != (const DNACharacter & a) const
			{
				return !(*this == a);
			}
		};

		typedef unrolled_list<DNACharacter, Size, 25> Sequence;
		typedef Sequence::iterator SequencePosIterator;
		typedef Sequence::reverse_iterator SequenceNegIterator;
		typedef Sequence::chunk_size PaddingInt;

		class StrandIterator: public std::iterator<std::bidirectional_iterator_tag, char, size_t>
		{
		public:
			void MakeInverted();
			void Swap(StrandIterator & toSwap);
			char operator * () const;
			StrandIterator Invert() const;
			Direction GetDirection() const;
			StrandIterator& operator ++ ();
			StrandIterator operator ++ (int);
			StrandIterator& operator -- ();
			StrandIterator operator -- (int);
			size_t GetElementId() const;
			size_t GetOriginalPosition() const;
			PaddingInt& GetPadding();
			const PaddingInt& GetPadding() const;
			SequencePosIterator Base() const;
			char TranslateChar(char ch) const;
			bool AtValidPosition() const;
			bool operator < (const StrandIterator & comp) const;
			bool operator == (const StrandIterator & comp) const;
			bool operator != (const StrandIterator & comp) const;
			StrandIterator();
			StrandIterator(const StrandIterator & toCopy);
			StrandIterator(SequencePosIterator base, Direction direction);
			StrandIterator& operator = (const StrandIterator & toCopy);
		private:
			friend class DNASequence;
			SequencePosIterator it_;
			Direction direction_;
		};

		typedef boost::function<void(StrandIterator, StrandIterator)> NotifyFunction;
		
		void Clear();
		size_t TotalSize() const;
		size_t ChrNumber() const;
		static char Translate(char ch);
		StrandIterator PositiveBegin(size_t chr) const;
		StrandIterator PositiveEnd(size_t chr) const;
		StrandIterator NegativeBegin(size_t chr) const;
		StrandIterator NegativeEnd(size_t chr) const;
		StrandIterator Begin(Direction, size_t chr) const;
		StrandIterator End(Direction, size_t chr) const;
		void Replace(StrandIterator source,
			size_t sourceDistance, 
			StrandIterator target,
			size_t targetDistance,
			NotifyFunction before = 0,
			NotifyFunction after = 0);
		explicit DNASequence(const std::vector<FASTARecord> & record);
		DNASequence(const std::vector<FASTARecord> & record, std::vector<std::vector<Pos> > & original, bool clear = false);
		std::pair<size_t, size_t> SpellOriginal(StrandIterator it1, StrandIterator it2) const;
		size_t GlobalIndex(StrandIterator it) const;		
		static const char UNKNOWN_BASE;
		static const std::string alphabet;			
	private:
		DISALLOW_COPY_AND_ASSIGN(DNASequence);	
		static const char SEPARATION_CHAR;
		const static char DELETED_CHAR;
		const static Pos DELETED_POS;
		static const std::string complementary_;

		struct IteratorPtrHash
		{
			size_t operator()(SequencePosIterator * it) const
			{
				return reinterpret_cast<size_t>(&(**it));
			}
		};

		struct IteratorPtrCompare
		{
			bool operator()(SequencePosIterator * it, SequencePosIterator * jt) const
			{
				return (*it) == (*jt);
			}
		};

		typedef boost::unordered_set<SequencePosIterator*, IteratorPtrHash, IteratorPtrCompare> IteratorMap;
		typedef IteratorMap::iterator IteratorPlace;
		typedef std::pair<IteratorPlace, IteratorPlace> IteratorRange;		

		void SubscribeIterator(SequencePosIterator & it);
		void UnsubscribeIterator(SequencePosIterator & it);
		void NotifyBefore(SequencePosIterator begin, SequencePosIterator end, NotifyFunction before);
		void NotifyAfter(SequencePosIterator begin, SequencePosIterator end, NotifyFunction after);
		
		Sequence sequence_;
		std::vector<SequencePosIterator> posBegin_;
		std::vector<SequencePosIterator> posEnd_;
		IteratorMap iteratorStore_;
		std::vector<IteratorPlace> toReplace_;
	};	
	
	inline bool ProperKMer(DNASequence::StrandIterator it, size_t k)
	{
		for(size_t i = 0; i < k; i++, ++it)
		{
			if(!it.AtValidPosition())
			{
				return false;
			}
		}

		return true;
	}
}

#endif
