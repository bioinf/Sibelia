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
			Pos pos;
			DNACharacter() {}
			DNACharacter(char actual): actual(actual), pos(actual) {}
			DNACharacter(char actual, Pos pos): actual(actual), pos(pos) {}
			bool operator == (const DNACharacter & a) const
			{
				return pos == a.pos && actual == a.actual;
			}

			bool operator != (const DNACharacter & a) const
			{
				return !(*this == a);
			}
		};

		typedef unrolled_list<DNACharacter, 1000> Sequence;
		typedef Sequence::iterator SequencePosIterator;
		typedef Sequence::reverse_iterator SequenceNegIterator;

	private:
		const static char DELETED_CHAR;
		const static Pos DELETED_POS;

		class GenericIterator
		{
		public:
			virtual ~GenericIterator();
			virtual char Spell() const = 0;			
			virtual void MoveForward() = 0;
			virtual void MoveBackward() = 0;
			virtual char TranslateChar(char ch) const = 0;
			virtual GenericIterator* Clone() const = 0;
			virtual Direction GetDirection() const = 0;
			virtual DNACharacter* GetNaked() const = 0;
			virtual GenericIterator* Invert() const = 0;
			virtual SequencePosIterator Base() const = 0;
			virtual bool Equal(const GenericIterator& toCompare) const = 0;
		};

		class ForwardIterator: public GenericIterator
		{
		public:
			char Spell() const;
			void MoveForward();
			void MoveBackward();
			char TranslateChar(char ch) const;
			Direction GetDirection() const;
			GenericIterator* Clone() const;
			GenericIterator* Invert() const;
			SequencePosIterator Base() const;
			DNACharacter* GetNaked() const;			
			bool Equal(const GenericIterator& toCompare) const;
			ForwardIterator();
			ForwardIterator(SequencePosIterator it);
		private:
			SequencePosIterator it_;
		};

		class BackwardIterator: public GenericIterator
		{
		public:
			char Spell() const;
			void MoveForward();
			void MoveBackward();
			char TranslateChar(char ch) const;
			Direction GetDirection() const;
			GenericIterator* Clone() const;
			GenericIterator* Invert() const;
			SequencePosIterator Base() const;
			SequenceNegIterator Natural() const;
			DNACharacter* GetNaked() const;
			bool Equal(const GenericIterator& toCompare) const;
			BackwardIterator();
			BackwardIterator(SequenceNegIterator it);
		private:
			SequenceNegIterator it_;
		};

	public:			
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
			SequencePosIterator Base() const;
			char TranslateChar(char ch) const;
			const DNACharacter* GetNaked() const;
			bool AtValidPosition() const;
			bool operator < (const StrandIterator & comp) const;
			bool operator == (const StrandIterator & comp) const;
			bool operator != (const StrandIterator & comp) const;
			StrandIterator();
			StrandIterator(GenericIterator * it);
			StrandIterator(const StrandIterator & toCopy);
			StrandIterator(SequencePosIterator base, Direction direction);
			StrandIterator& operator = (const StrandIterator & toCopy);
		private:
			friend class DNASequence;
			std::auto_ptr<GenericIterator> it_;			
		};
		
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
			Sequence::notify_func before = 0,
			Sequence::notify_func after = 0);
		explicit DNASequence(const std::vector<FASTARecord> & record);
		DNASequence(const std::vector<FASTARecord> & record, const std::vector<std::vector<Pos> > & original);
		std::pair<size_t, size_t> SpellOriginal(StrandIterator it1, StrandIterator it2) const;
		
		static const char UNKNOWN_BASE;
		static const std::string alphabet;		
	private:
		DISALLOW_COPY_AND_ASSIGN(DNASequence);	
		static const char SEPARATION_CHAR;
		static const std::string complementary_;
		
		Sequence sequence_;
		std::vector<SequencePosIterator> posBegin_;
		std::vector<SequencePosIterator> posEnd_;
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
