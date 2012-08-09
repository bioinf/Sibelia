#ifndef _DNA_SEQUENCE_H_
#define _DNA_SEQUENCE_H_

#include "common.h"

#pragma warning(disable:4355)

namespace SyntenyBuilder
{
	class DNASequence
	{
	public:
		enum Direction
		{
			positive,
			negative
		};

	private:
		static const size_t NO_POS;

		struct DNACharacter
		{
			char actual;
			size_t pos;
			DNACharacter() {}
			DNACharacter(char actual): actual(actual), pos(actual) {}
			DNACharacter(char actual, size_t pos): actual(actual), pos(pos) {}
		};

		typedef std::list<DNACharacter> Sequence;
		typedef Sequence::const_iterator SequencePosIterator;
		typedef Sequence::const_reverse_iterator SequenceNegIterator;

		class GenericIterator
		{
		public:
			virtual ~GenericIterator();
			virtual char Spell() const = 0;			
			virtual void MoveForward() = 0;
			virtual void MoveBackward() = 0;
			virtual GenericIterator* Clone() const = 0;
			virtual Direction GetDirection() const = 0;
			virtual const DNACharacter* GetNaked() const = 0;			
			virtual bool Equal(const GenericIterator& toCompare) const = 0;
		};

		class ForwardIterator: public GenericIterator
		{
		public:
			char Spell() const;
			void MoveForward();
			void MoveBackward();
			Direction GetDirection() const;
			GenericIterator* Clone() const;
			const DNACharacter* GetNaked() const;
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
			Direction GetDirection() const;
			GenericIterator* Clone() const;
			const DNACharacter* GetNaked() const;
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
			bool operator < (const StrandIterator & comp) const;
			bool operator == (const StrandIterator & comp) const;
			bool operator != (const StrandIterator & comp) const;
			StrandIterator();
			StrandIterator(GenericIterator * it);
			StrandIterator(const StrandIterator & toCopy);
			StrandIterator& operator = (const StrandIterator & toCopy);
		private:
			std::auto_ptr<GenericIterator> it_;
		};

		DNASequence(const std::string & sequence);
		StrandIterator PositiveBegin() const;
		StrandIterator PositiveEnd() const;
		StrandIterator NegativeBegin() const;
		StrandIterator NegativeEnd() const;
		
		template<class Iterator>
			std::pair<size_t, size_t> SpellOriginal(StrandIterator it1, StrandIterator it2, Iterator out) const
			{
				size_t pos1 = NO_POS;
				size_t pos2 = NO_POS;
				for(;it1.GetOriginalPosition() == NO_POS && it1 != it2; ++it1);
				for(--it2; it2.GetOriginalPosition() == NO_POS && it1 != it2; --it2);
				size_t start = std::min(pos1, pos2);
				size_t end = std::max(pos1, pos2) + 1;
				if(it1.GetDirection() == positive)
				{
					std::copy(original_.begin() + start, original_.end() + end, out);
				}
				else
				{
					std::string rcomp;
					std::copy(original_.begin() + start, original_.end() + end, std::back_inserter(rcomp));
					for(size_t i = 0; i < rcomp.size(); i++)
					{
						rcomp[i] = Translate(rcomp[i]);
					}

					std::copy(rcomp.rbegin(), rcomp.rend(), out);
				}

				return std::make_pair(start, end);
			}

		template<class Iterator>
			void CopyN(Iterator start, size_t count, StrandIterator out)
			{				
			}

		void EraseN(StrandIterator out, size_t count)
		{
		}

		template<class Iterator>
			void SpellRaw(Iterator out) const
			{
			}

		size_t Size() const;
		static const std::string alphabet;
	private:
		DISALLOW_COPY_AND_ASSIGN(DNASequence);	
		static char Translate(char ch);
		static const std::string complementary_;		
		Sequence sequence_;
		std::string original_;
	};	
}

#endif
