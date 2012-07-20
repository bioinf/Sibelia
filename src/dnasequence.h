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
		struct ReadingStrategy;
	public:			
		class StrandIterator;

		enum Direction
		{
			positive,
			negative
		};

		class StrandIterator: public std::iterator<std::bidirectional_iterator_tag, char, size_t>
		{
		public:
			typedef ReadingStrategy RStrategy;

 			StrandIterator() {}
			StrandIterator(IndexConstIterator it, const RStrategy * rStrategy): it_(it), rStrategy_(rStrategy) {}

			char operator * () const
			{
				return rStrategy_->GetBase(it_);
			}			

			void Jump(size_t count)
			{
				rStrategy_->Jump(it_, count);
			}

			bool ProperKMer(size_t k) const
			{
				StrandIterator temp(*this);
				temp.Jump(k - 1);
				return temp.Valid();
			}

			StrandIterator Invert() const
			{
				return StrandIterator(it_, rStrategy_->Invert());
			}

			Direction GetDirection()
			{
				return rStrategy_->GetDirection();
			}

			size_t GetHashCode(size_t strSize)
			{
				return rStrategy_->GetHash(it_, strSize);
			}

			StrandIterator& operator ++ ()
			{
				rStrategy_->AdvanceForward(it_);
				return *this;
			}

			StrandIterator operator ++ (int)
			{
				StrandIterator ret(it_, rStrategy_);
				rStrategy_->AdvanceForward(it_);
				return ret;
			}

			StrandIterator& operator -- ()
			{
				rStrategy_->AdvanceBackward(it_);
				return *this;
			}

			StrandIterator operator -- (int)
			{
				StrandIterator ret(it_, rStrategy_);
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

			size_t GetPosition() const
			{
				return it_.GetPosition();
			}

			const ReadingStrategy* GetStrategy() const
			{
				return rStrategy_;
			}

		private:
			IndexConstIterator it_;
			const RStrategy * rStrategy_;
		};

		DNASequence(const std::string sequence):
			substrSize_(NOT_KEEP),
			sequence_(sequence),
			original_(sequence),
			positiveReading_(this),
			negativeReading_(this),
			deletions_(0)
		{
			position_.resize(sequence_.size());
			for(size_t i = 0; i < sequence_.size(); i++)
			{
				position_[i] = i;
			}
		}

		StrandIterator PositiveByIndex(size_t pos) const
		{
			return StrandIterator(IndexConstIterator(sequence_, pos, DELETED_CHARACTER), &positiveReading_);
		}				

		StrandIterator NegativeByIndex(size_t pos) const
		{
			return StrandIterator(IndexConstIterator(sequence_, pos, DELETED_CHARACTER, reverse), &negativeReading_);
		}

		StrandIterator PositiveBegin() const
		{
			return StrandIterator(IndexConstIterator(sequence_, 0, DELETED_CHARACTER), &positiveReading_);
		}

		StrandIterator PositiveRightEnd() const
		{
			return StrandIterator(MakeRightEnd(sequence_, DELETED_CHARACTER), &positiveReading_);
		}

		StrandIterator PositiveLeftEnd() const
		{
			return StrandIterator(MakeLeftEnd(sequence_, DELETED_CHARACTER), &positiveReading_);
		}

		StrandIterator NegativeBegin() const
		{
			return StrandIterator(IndexConstIterator(sequence_, sequence_.size() - 1, DELETED_CHARACTER, reverse),
				&negativeReading_);
		}

		StrandIterator NegativeRightEnd() const
		{
			return StrandIterator(MakeLeftEnd(sequence_, DELETED_CHARACTER), &negativeReading_);
		}

		StrandIterator NegativeLeftEnd() const
		{
			return StrandIterator(MakeLeftEnd(sequence_, DELETED_CHARACTER), &negativeReading_);
		}

		template<class Iterator>
			std::pair<size_t, size_t> SpellOriginal(StrandIterator it1, StrandIterator it2, Iterator out) const
			{
				StrandIterator out1(IndexConstIterator(original_, position_[it1.GetPosition()], DELETED_CHARACTER), it1.GetStrategy());
				StrandIterator out2(IndexConstIterator(original_, position_[(--it2).GetPosition()], DELETED_CHARACTER), it2.GetStrategy());
				std::copy(out1, ++out2, out);
				size_t pos1 = out1.GetPosition();
				size_t pos2 = (--out2).GetPosition();
				return std::make_pair(std::min(pos1, pos2), std::max(pos1, pos2) + 1);
			}

		template<class Iterator>
			void Copy(StrandIterator out, Iterator start, Iterator end)
			{

			}

		void Erase(StrandIterator start, StrandIterator end)
		{

		}

		void KeepHash(size_t substrSize)
		{
			substrSize_ = substrSize;
			positiveHash_.resize(Size(), INVALID_HASH);
			negativeHash_.resize(Size(), INVALID_HASH);
		}

		void DropHash()
		{
			substrSize_ = NOT_KEEP;
			positiveHash_.clear();
			negativeHash_.clear();
		}

		void Optimize()
		{
			DropHash();
			position_.erase(std::remove(position_.begin(), position_.end(), DELETED_CHARACTER), position_.end());
			sequence_.erase(std::remove(sequence_.begin(), sequence_.end(), DELETED_CHARACTER), sequence_.end());
		}

		template<class Iterator>
			void SpellRaw(Iterator out)
			{
				std::copy(sequence_.begin(), sequence_.end(), out);
			}

		size_t Size() const
		{
			return sequence_.size();
		}

		static const std::string alphabet;
	private:
		DISALLOW_COPY_AND_ASSIGN(DNASequence);	
		typedef long long hash_t;
		static const std::string complementary_;
		static const size_t MOD;
		static const size_t NOT_KEEP;
		static const size_t HASH_BASE;
		static const char INVALID_HASH;
		static const char DELETED_CHARACTER;

		struct ReadingStrategy
		{
		public:
			ReadingStrategy(DNASequence * sequence): sequence_(sequence) {}
			char GetBase(IndexConstIterator it) const
			{
				return Translate(*it);
			}
							
			virtual char Translate(char ch) const = 0;
			virtual Direction GetDirection() const = 0;
			virtual const ReadingStrategy* Invert() const = 0;
			virtual void AdvanceForward(IndexConstIterator & it) const = 0;
			virtual void AdvanceBackward(IndexConstIterator & it) const = 0;
			virtual void Jump(IndexConstIterator & it, size_t count) const = 0;
			virtual size_t GetHash(IndexConstIterator it, size_t strSize) const = 0;				
		protected:
			DNASequence * sequence_;
		};			

		struct PositiveReadingStrategy: public ReadingStrategy
		{
			PositiveReadingStrategy(DNASequence * sequence): ReadingStrategy(sequence) {}
			Direction GetDirection() const
			{
				return positive;
			}

			char Translate(char ch) const
			{
				return ch;
			}

			void AdvanceForward(IndexConstIterator & it) const
			{
				++it;
			}

			void AdvanceBackward(IndexConstIterator & it) const
			{
				--it;
			}

			const ReadingStrategy* Invert() const
			{
				return &sequence_->negativeReading_;
			}

			void Jump(IndexConstIterator & it, size_t count) const
			{
				if(sequence_->deletions_ == 0)
				{
					it = IndexConstIterator(sequence_->sequence_, 
						std::min(it.GetPosition() + count, sequence_->sequence_.size()),
						sequence_->DELETED_CHARACTER);
				}
				else
				{
					for(size_t i = 0; i < count; i++)
					{
						AdvanceForward(it);
					}
				}
			}

			size_t GetHash(IndexConstIterator it, size_t strSize) const
			{
				hash_t ret;
				if(this->sequence_->substrSize_ == strSize)
				{
					ret = this->sequence_->positiveHash_[it.GetPosition()];
					assert(ret == this->sequence_->CalcHash(it, strSize));
				}
				else
				{
					ret = this->sequence_->CalcHash(it, strSize);
				}

				return static_cast<size_t>(ret);
			}
		};

		struct NegativeReadingStrategy: public ReadingStrategy
		{
			NegativeReadingStrategy(DNASequence * sequence): ReadingStrategy(sequence) {}
			Direction GetDirection() const
			{
				return negative;
			}

			char Translate(char ch) const
			{
				return complementary_[ch];
			}

			void AdvanceForward(IndexConstIterator & it) const
			{
				--it;
			}

			void AdvanceBackward(IndexConstIterator & it) const
			{
				++it;
			}

			void Jump(IndexConstIterator & it, size_t count) const
			{
				if(sequence_->deletions_ == 0)
				{
					it = IndexConstIterator(sequence_->sequence_, 
						it.GetPosition() < count ? IndexConstIterator::NPOS : it.GetPosition() - count,
						sequence_->DELETED_CHARACTER,
						reverse);

				}
				else
				{
					for(size_t i = 0; i < count; i++)
					{
						AdvanceForward(it);
					}
				}
			}

			const ReadingStrategy* Invert() const
			{
				return &sequence_->positiveReading_;
			}

			size_t GetHash(IndexConstIterator it, size_t strSize) const
			{
				hash_t ret;
				StrandIterator src = this->sequence_->NegativeByIndex(it.GetPosition());
				if(this->sequence_->substrSize_ == strSize)
				{
					ret = this->sequence_->negativeHash_[it.GetPosition()];
					assert(ret == this->sequence_->CalcHash(src, strSize));
				}
				else
				{
					ret = this->sequence_->CalcHash(src, strSize);
				}

				return static_cast<size_t>(ret);
			}
		};
		
		template<class Iterator>
			static hash_t CalcHash(Iterator it, size_t k)
			{
				size_t base = 1;
				size_t hash = 0;
				for(size_t i = 0; i < k; i++)
				{			
					hash += *it++ * base;
					base *= HASH_BASE;
				}		

				return static_cast<hash_t>(hash % MOD);
			}

		PositiveReadingStrategy positiveReading_;
		NegativeReadingStrategy negativeReading_;

		//Current version of the sequence (after possible simplification)		
		std::string sequence_;
		//Original version of the sequence (before doing any modifications)
		std::string original_;
		//Map from each position in the current sequence to the original sequence
		size_t substrSize_;
		std::vector<hash_t> positiveHash_;
		std::vector<hash_t> negativeHash_;
		std::vector<size_t> position_;
		size_t deletions_;
	};	
}

#endif
