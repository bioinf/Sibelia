#ifndef _DNA_SEQUENCE_H_
#define _DNA_SEQUENCE_H_

#include "common.h"
#include "indexiterator.h"

#pragma warning(disable:4355)

namespace SyntenyBuilder
{
	class DNASequence
	{
	private:
		struct WritingStrategy;
		template<class T> struct ReadingStrategy;
	public:			
		class StrandIterator;
		class StrandConstIterator;

		enum Direction
		{
			positive,
			negative
		};

		class StrandIterator: public std::iterator<std::bidirectional_iterator_tag, char, int>
		{
		public:
			typedef WritingStrategy WStrategy;
			typedef ReadingStrategy<IndexIterator> RStrategy;

			StrandIterator() {}
			StrandIterator(IndexIterator it, const RStrategy * rStrategy, const WStrategy * wStrategy):
				it_(it), rStrategy_(rStrategy), wStrategy_(wStrategy) {}

			Direction GetDirection()
			{
				return rStrategy_->GetDirection();
			}

			char operator * () const
			{
				return rStrategy_->GetBase(it_);
			}			

			StrandIterator& operator ++ ()
			{
				rStrategy_->AdvanceForward(it_);
				return *this;
			}

			StrandIterator operator ++ (int)
			{
				StrandIterator ret(it_, rStrategy_, wStrategy_);
				rStrategy_->AdvanceForward(it_);
				return ret;
			}

			void Invalidate()
			{
				wStrategy_->Invalidate(it_);
				++(*this);
			}

			void AssignBase(char newBase)
			{
				wStrategy_->AssignBase(it_, newBase);
			}

			StrandIterator& operator -- ()
			{
				rStrategy_->AdvanceBackward(it_);
				return *this;
			}

			StrandIterator operator -- (int)
			{
				StrandIterator ret(it_, rStrategy_, wStrategy_);
				rStrategy_->AdvanceBackward(it_);
				return ret;
			}

			bool operator == (StrandIterator comp) const
			{
				return it_ == comp.it_; 
			}

			bool operator != (StrandIterator comp) const
			{
				return !(*this == comp);
			}
			
			bool Valid() const
			{
				return it_.Valid();
			}

			size_t GetHashCode(int strSize)
			{
				return rStrategy_->GetHash(it_, strSize);
			}

			int GetPosition() const
			{
				return it_.GetPosition();
			}

		private:
			IndexIterator it_;
			const RStrategy * rStrategy_;
			const WStrategy * wStrategy_;
			friend class StrandConstIterator;
		};

		class StrandConstIterator: public std::iterator<std::bidirectional_iterator_tag, char, int>
		{
		public:
			typedef ReadingStrategy<IndexConstIterator> RStrategy;

 			StrandConstIterator() {}
			StrandConstIterator(IndexConstIterator it, const RStrategy * rStrategy): it_(it), rStrategy_(rStrategy) {}
			StrandConstIterator(StrandIterator other): it_(other.it_), rStrategy_(other.rStrategy_->Convert()) {}

			char operator * () const
			{
				return rStrategy_->GetBase(it_);
			}			

			Direction GetDirection()
			{
				return rStrategy_->GetDirection();
			}

			size_t GetHashCode(int strSize)
			{
				return rStrategy_->GetHash(it_, strSize);
			}

			StrandConstIterator& operator ++ ()
			{
				rStrategy_->AdvanceForward(it_);
				return *this;
			}

			StrandConstIterator operator ++ (int)
			{
				StrandConstIterator ret(it_, rStrategy_);
				rStrategy_->AdvanceForward(it_);
				return ret;
			}

			StrandConstIterator& operator -- ()
			{
				rStrategy_->AdvanceBackward(it_);
				return *this;
			}

			StrandConstIterator operator -- (int)
			{
				StrandConstIterator ret(it_, rStrategy_);
				rStrategy_->AdvanceBackward(it_);
				return ret;
			}

			bool operator == (StrandConstIterator comp) const
			{
				return it_ == comp.it_; 
			}

			bool operator != (StrandConstIterator comp) const
			{
				return !(*this == comp);
			}
			
			bool Valid() const
			{
				return it_.Valid();
			}

			int GetPosition() const
			{
				return it_.GetPosition();
			}

		private:
			IndexConstIterator it_;
			const RStrategy * rStrategy_;
		};

		typedef boost::function<void (void) > RefreshHook;
		typedef boost::function<void (int pos, bool)> ModifyHook;
		
		DNASequence(const std::string sequence, RefreshHook refresh = RefreshHook(), 
			ModifyHook before = ModifyHook(), ModifyHook after = ModifyHook()): 
			substrSize_(INVALID_HASH),
			sequence_(sequence),
			original_(sequence),
			refresh_(refresh),
			before_(before),
			after_(after),
			positiveReading_(this),
			negativeReading_(this),
			positiveConstReading_(this),
			negativeConstReading_(this),
			positiveWriting_(this),
			negativeWriting_(this)
		{
			position_.resize(sequence_.size());
			for(int i = 0; i < static_cast<int>(sequence_.size()); i++)
			{
				position_[i] = i;
			}
		}

		StrandIterator PositiveByIndex(int pos)
		{
			return StrandIterator(IndexIterator(sequence_, pos, DELETED_CHARACTER), &positiveReading_, &positiveWriting_);
		}

		StrandConstIterator PositiveByIndex(int pos) const
		{
			return StrandConstIterator(IndexConstIterator(sequence_, pos, DELETED_CHARACTER), &positiveConstReading_);
		}				

		StrandIterator NegativeByIndex(int pos)
		{
			return StrandIterator(IndexIterator(sequence_, pos, DELETED_CHARACTER, reverse), &negativeReading_, &negativeWriting_);
		}

		StrandConstIterator NegativeByIndex(int pos) const
		{
			return StrandConstIterator(IndexConstIterator(sequence_, pos, DELETED_CHARACTER, reverse), &negativeConstReading_);
		}

		StrandIterator PositiveBegin()
		{
			return StrandIterator(IndexIterator(sequence_, 0, DELETED_CHARACTER), &positiveReading_, &positiveWriting_);
		}

		StrandConstIterator PositiveBegin() const
		{
			return StrandConstIterator(IndexConstIterator(sequence_, 0, DELETED_CHARACTER), &positiveConstReading_);
		}

		StrandIterator PositiveRightEnd() 
		{
			return StrandIterator(MakeRightEnd(sequence_, DELETED_CHARACTER), &positiveReading_, &positiveWriting_);
		}

		StrandConstIterator PositiveRightEnd() const
		{
			return StrandConstIterator(MakeRightEnd(sequence_, DELETED_CHARACTER), &positiveConstReading_);
		}

		StrandIterator PositiveLeftEnd() 
		{
			return StrandIterator(MakeLeftEnd(sequence_, DELETED_CHARACTER), &positiveReading_, &positiveWriting_);
		}

		StrandConstIterator PositiveLeftEnd() const
		{
			return StrandConstIterator(MakeLeftEnd(sequence_, DELETED_CHARACTER), &positiveConstReading_);
		}

		StrandIterator NegativeBegin()
		{
			return StrandIterator(IndexIterator(sequence_, static_cast<int>(sequence_.size() - 1), DELETED_CHARACTER, reverse), 
				&negativeReading_, &negativeWriting_);
		}

		StrandConstIterator NegativeBegin() const
		{
			return StrandConstIterator(IndexConstIterator(sequence_, static_cast<int>(sequence_.size() - 1), DELETED_CHARACTER, reverse),
				&negativeConstReading_);
		}

		StrandIterator NegativeRightEnd()
		{
			return StrandIterator(MakeLeftEnd(sequence_, DELETED_CHARACTER), &negativeReading_, &negativeWriting_);
		}

		StrandConstIterator NegativeRightEnd() const
		{
			return StrandConstIterator(MakeLeftEnd(sequence_, DELETED_CHARACTER), &negativeConstReading_);
		}

		StrandIterator NegativeLeftEnd()
		{
			return StrandIterator(MakeLeftEnd(sequence_, DELETED_CHARACTER), &negativeReading_, &negativeWriting_);
		}

		StrandConstIterator NegativeLeftEnd() const
		{
			return StrandConstIterator(MakeLeftEnd(sequence_, DELETED_CHARACTER), &negativeConstReading_);
		}

		//This is a bad method (workaround-like). It must be replaced with a smarter representation later
		template<class Iterator>
			std::pair<int, int> SpellOriginal(StrandConstIterator it1, StrandConstIterator it2, Iterator out)
			{
				std::pair<int, int> ret(0, 0);
				int pos1 = it1.GetPosition();
				int pos2 = it2.GetPosition();
				if(pos1 > pos2)
				{
					pos1 = position_[pos1];
					pos2 = position_[pos2 + 1];
					ret = std::make_pair(pos2, pos1 + 1);
					for(; pos1 >= pos2; pos1--)
					{
						*out++ = complementary_[original_[pos1]];
					}
				}
				else
				{
					pos1 = position_[pos1];
					pos2 = position_[pos2 - 1];
					ret = std::make_pair(pos1, pos2 + 1);
					for(; pos1 <= pos2; pos1++)
					{
						*out++ = original_[pos1];
					}
				}

				return ret;
			}

		void KeepHash(int substrSize)
		{
			substrSize_ = substrSize;
			positiveHash_.resize(Size(), INVALID_HASH);
			negativeHash_.resize(Size(), INVALID_HASH);
		}

		void Optimize()
		{
			substrSize_ = INVALID_HASH;
			position_.erase(std::remove(position_.begin(), position_.end(), DELETED_CHARACTER), position_.end());
			sequence_.erase(std::remove(sequence_.begin(), sequence_.end(), DELETED_CHARACTER), sequence_.end());
			if(!refresh_.empty())
			{
				refresh_();
			}
		}

		template<class Iterator>
			void SpellRaw(Iterator out)
			{
				std::copy(sequence_.begin(), sequence_.end(), out);
			}

		int Size() const
		{
			return static_cast<int>(sequence_.size());
		}

		static const std::string alphabet;
	private:
		DISALLOW_COPY_AND_ASSIGN(DNASequence);	
		typedef long long hash_t;
		static const std::string complementary_;
		static const size_t MOD;
		static const size_t HASH_BASE;
		static const char INVALID_HASH;
		static const char DELETED_CHARACTER;

		template<class Iterator>
			struct ReadingStrategy
			{
			public:
				ReadingStrategy(DNASequence * sequence): sequence_(sequence) {}
				virtual char GetBase(Iterator it) const = 0;
				virtual Direction GetDirection() const = 0;
				virtual void AdvanceForward(Iterator & it) const = 0;
				virtual void AdvanceBackward(Iterator & it) const = 0;
				virtual size_t GetHash(Iterator it, int strSize) const = 0;				
				virtual const ReadingStrategy<IndexConstIterator>* Convert() const = 0;				
			protected:
				DNASequence * sequence_;
			};			

		struct WritingStrategy
		{
		public:
			WritingStrategy(DNASequence * sequence): sequence_(sequence) {}
			virtual void AssignBase(IndexIterator it, char newBase) const = 0;
			virtual void Invalidate(IndexIterator it) const
			{
				if(!sequence_->before_.empty())
				{
					sequence_->before_(it.GetPosition(), true);
				}

				*it = DELETED_CHARACTER;
				sequence_->position_[it.GetPosition()] = DELETED_CHARACTER;
				InvalidateHash(it.GetPosition());

				if(!sequence_->after_.empty())
				{
					sequence_->after_(it.GetPosition(), true);
				}
			}

			virtual void InvalidateHash(int pos) const
			{
				if(sequence_->substrSize_ != INVALID_HASH)
				{
					for(int i = 0; i < sequence_->substrSize_; i++)
					{
						if(pos - i >= 0)
						{
							sequence_->positiveHash_[pos - i] = INVALID_HASH;
						}

						if(pos + i < sequence_->Size())
						{
							sequence_->negativeHash_[pos + i] = INVALID_HASH;
						}
					}
				}
			}

		protected:
			DNASequence * sequence_;
		};		

		template<class Iterator>
			struct PositiveReadingStrategy: public ReadingStrategy<Iterator>
			{
				PositiveReadingStrategy(DNASequence * sequence): ReadingStrategy<Iterator>(sequence) {}
				Direction GetDirection() const
				{
					return positive;
				}

				char GetBase(Iterator it) const
				{
					return *it;
				}

				void AdvanceForward(Iterator & it) const
				{
					++it;
				}

				void AdvanceBackward(Iterator & it) const
				{
					--it;
				}

				const ReadingStrategy<IndexConstIterator>* Convert() const
				{
					return &this->sequence_->positiveConstReading_;
				}

				size_t GetHash(Iterator it, int strSize) const
				{
					hash_t ret;
					if(this->sequence_->substrSize_ == strSize)
					{
						if(this->sequence_->positiveHash_[it.GetPosition()] == INVALID_HASH)
						{
							this->sequence_->positiveHash_[it.GetPosition()] = this->sequence_->CalcHash(it, strSize);
						}

						ret = this->sequence_->positiveHash_[it.GetPosition()];
					}
					else
					{
						ret = this->sequence_->CalcHash(it, strSize);
					}

					return static_cast<size_t>(ret);
				}
			};

		struct PositiveWritingStrategy: public WritingStrategy
		{			
			PositiveWritingStrategy(DNASequence * sequence): WritingStrategy(sequence) {}
			void AssignBase(IndexIterator it, char newBase) const
			{
				if(!sequence_->before_.empty())
				{
					sequence_->before_(it.GetPosition(), false);
				}

				*it = newBase;
				InvalidateHash(it.GetPosition());

				if(!sequence_->after_.empty())
				{
					sequence_->after_(it.GetPosition(), false);
				}
			}
		};

		template<class Iterator>
			struct NegativeReadingStrategy: public ReadingStrategy<Iterator>
			{
				NegativeReadingStrategy(DNASequence * sequence): ReadingStrategy<Iterator>(sequence) {}
				Direction GetDirection() const
				{
					return negative;
				}

				char GetBase(Iterator it) const
				{
					return complementary_[*it];
				}

				void AdvanceForward(Iterator & it) const
				{
					--it;
				}

				void AdvanceBackward(Iterator & it) const
				{
					++it;
				}

				const ReadingStrategy<IndexConstIterator>* Convert() const
				{
					return &this->sequence_->negativeConstReading_;
				}

				size_t GetHash(Iterator it, int strSize) const
				{
					hash_t ret;
					StrandConstIterator src = this->sequence_->NegativeByIndex(it.GetPosition());
					if(this->sequence_->substrSize_ == strSize)
					{
						if(this->sequence_->negativeHash_[it.GetPosition()] == INVALID_HASH)
						{							
							this->sequence_->negativeHash_[it.GetPosition()] = this->sequence_->CalcHash(src, strSize);
						}

						ret = this->sequence_->negativeHash_[it.GetPosition()];
					}
					else
					{
						ret = this->sequence_->CalcHash(src, strSize);
					}

					return static_cast<size_t>(ret);
				}
			};
		
		struct NegativeWritingStrategy: public WritingStrategy
		{
			NegativeWritingStrategy(DNASequence * sequence): WritingStrategy(sequence) {}
			void AssignBase(IndexIterator it, char newBase) const
			{
				if(!sequence_->before_.empty())
				{
					sequence_->before_(it.GetPosition(), false);
				}

				*it = complementary_[newBase];
				InvalidateHash(it.GetPosition());

				if(!sequence_->after_.empty())
				{
					sequence_->after_(it.GetPosition(), false);
				}
			}
		};						

		template<class Iterator>
			static hash_t CalcHash(Iterator it, int k)
			{
				size_t base = 1;
				size_t hash = 0;
				for(int i = 0; i < k; i++)
				{			
					hash += *it++ * base;
					base *= HASH_BASE;
				}		

				return static_cast<hash_t>(hash % MOD);
			}

		RefreshHook refresh_;
		PositiveWritingStrategy positiveWriting_;
		NegativeWritingStrategy negativeWriting_;
		PositiveReadingStrategy<IndexIterator> positiveReading_;		
		NegativeReadingStrategy<IndexIterator> negativeReading_;		
		PositiveReadingStrategy<IndexConstIterator> positiveConstReading_;
		NegativeReadingStrategy<IndexConstIterator> negativeConstReading_;

		//Current version of the sequence (after possible simplification)		
		std::string sequence_;
		//Original version of the sequence (before doing any modifications)
		std::string original_;
		//Map from each position in the current sequence to the original sequence
		int substrSize_;
		ModifyHook before_;
		ModifyHook after_;
		std::vector<hash_t> positiveHash_;
		std::vector<hash_t> negativeHash_;
		std::vector<int> position_;
	};	
}

#endif
