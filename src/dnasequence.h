#ifndef _DNA_SEQUENCE_H_
#define _DNA_SEQUENCE_H_

#include "common.h"
#include "indexiterator.h"

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

		class StrandIterator: public std::iterator<std::bidirectional_iterator_tag, char, void>
		{
		public:
			typedef WritingStrategy WStrategy;
			typedef ReadingStrategy<IndexIterator> RStrategy;

			StrandIterator() {}
			StrandIterator(IndexIterator it, const RStrategy * rStrategy, const WStrategy * wStrategy):
				it_(it), rStrategy_(rStrategy), wStrategy_(wStrategy) {}

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

		class StrandConstIterator: public std::iterator<std::bidirectional_iterator_tag, char, void>
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
		typedef boost::function<void (int pos)> ModifyHook;
		
		DNASequence(const std::string sequence, RefreshHook refresh = RefreshHook(), 
			ModifyHook before = ModifyHook(), ModifyHook after = ModifyHook()): 
			sequence_(sequence),
			original_(sequence),
			refresh_(refresh),
			positiveWriting_(before, after),
			negativeWriting_(before, after)
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
			return StrandIterator(IndexIterator(sequence_, sequence_.size() - 1, DELETED_CHARACTER, reverse), 
				&negativeReading_, &negativeWriting_);
		}

		StrandConstIterator NegativeBegin() const
		{
			return StrandConstIterator(IndexConstIterator(sequence_, sequence_.size() - 1, DELETED_CHARACTER, reverse),
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

		void Optimize()
		{
			for(size_t i = 0; i < sequence_.size(); i++)
			{
				if(sequence_[i] == DELETED_CHARACTER)
				{
					position_[i] = DELETED_CHARACTER;
				}
			}

			position_.erase(std::remove(position_.begin(), position_.end(), DELETED_CHARACTER), position_.end());
			sequence_.erase(std::remove(sequence_.begin(), sequence_.end(), DELETED_CHARACTER), sequence_.end());
			if(!refresh_.empty())
			{
				refresh_();
			}
		}

		int Size() const
		{
			return static_cast<int>(sequence_.size());
		}

		static const std::string alphabet;
	private:
		DISALLOW_COPY_AND_ASSIGN(DNASequence);	
		static const std::string complementary_;
		static const char DELETED_CHARACTER = -1;

		template<class Iterator>
			struct ReadingStrategy
			{
				virtual char GetBase(Iterator it) const = 0;
				virtual void AdvanceForward(Iterator & it) const = 0;
				virtual void AdvanceBackward(Iterator & it) const = 0;
				virtual const ReadingStrategy<IndexConstIterator>* Convert() const = 0;
			};			

		struct WritingStrategy
		{
		public:
			WritingStrategy(ModifyHook before, ModifyHook after): before_(before), after_(after) {}
			virtual void Invalidate(IndexIterator it) const = 0;
			virtual void AssignBase(IndexIterator it, char newBase) const = 0;
		protected:
			ModifyHook before_;
			ModifyHook after_;
		};		

		template<class Iterator>
			struct PositiveReadingStrategy: public ReadingStrategy<Iterator>
			{
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
					return &DNASequence::positiveConstReading_;
				}
			};

		struct PositiveWritingStrategy: public WritingStrategy
		{
			void AssignBase(IndexIterator it, char newBase) const
			{
				if(!before_.empty())
				{
					before_(it.GetPosition());
				}

				*it = newBase;

				if(!after_.empty())
				{
					after_(it.GetPosition());
				}
			}

			void Invalidate(IndexIterator it) const
			{
				if(!before_.empty())
				{
					before_(it.GetPosition());
				}

				*it = DELETED_CHARACTER;

				if(!after_.empty())
				{
					after_(it.GetPosition());
				}
			}

			PositiveWritingStrategy(ModifyHook before, ModifyHook after): WritingStrategy(before, after) {}
		};

		template<class Iterator>
			struct NegativeReadingStrategy: public ReadingStrategy<Iterator>
			{
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
					return &DNASequence::negativeConstReading_;
				}
			};
		
		struct NegativeWritingStrategy: public WritingStrategy
		{
			void AssignBase(IndexIterator it, char newBase) const
			{
				if(!before_.empty())
				{
					before_(it.GetPosition());
				}

				*it = complementary_[newBase];

				if(!after_.empty())
				{
					after_(it.GetPosition());
				}
			}

			void Invalidate(IndexIterator it) const
			{
				if(!before_.empty())
				{
					before_(it.GetPosition());
				}

				*it = DELETED_CHARACTER;

				if(!after_.empty())
				{
					after_(it.GetPosition());
				}
			}

			NegativeWritingStrategy(ModifyHook before, ModifyHook after): WritingStrategy(before, after) {}
		};						

		RefreshHook refresh_;
		PositiveWritingStrategy positiveWriting_;
		NegativeWritingStrategy negativeWriting_;
		static PositiveReadingStrategy<IndexIterator> positiveReading_;		
		static NegativeReadingStrategy<IndexIterator> negativeReading_;		
		static PositiveReadingStrategy<IndexConstIterator> positiveConstReading_;
		static NegativeReadingStrategy<IndexConstIterator> negativeConstReading_;

		//Current version of the sequence (after possible simplification)		
		std::string sequence_;
		//Original version of the sequence (before doing any modifications)
		std::string original_;
		//Map from each position in the current sequence to the original sequence
		std::vector<int> position_;
	};	
}

#endif
