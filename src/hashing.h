#ifndef _HASHING_H_
#define _HASHING_H_

#include "dnasequence.h"

namespace SyntenyBuilder
{
	template<class Iterator>
		class SlidingWindow
		{
		public:
			static const size_t HASH_BASE;	
			SlidingWindow() {}
			SlidingWindow(Iterator kMerStart, Iterator end, size_t k): k_(k), highPow_(1), 
				kMerStart_(kMerStart), kMerEnd_(AdvanceForward(kMerStart, k - 1)), end_(end)
			{
				for(size_t i = 1; i < k; i++)
				{
					highPow_ *= HASH_BASE;
				}

				value_ = CalcKMerHash(kMerStart, k);
			}

			size_t GetValue() const
			{
				return value_;
			}

			size_t GetK() const
			{
				return k_;
			}

			Iterator GetBegin() const
			{
				return kMerStart_;
			}
			
			Iterator GetEnd() const
			{
				return ++StrandIterator(kMerEnd_);
			}

			bool Move()
			{
				value_ = (value_ - *kMerStart_ * highPow_) * HASH_BASE;
				++kMerStart_;
				++kMerEnd_;
				if(Valid())
				{
					value_ += *kMerEnd_;
					assert(value_ == CalcKMerHash(kMerStart_, k_));
					return true;
				}

				return false;
			}

			bool Valid() const
			{
				return kMerEnd_ != end_;
			}

			static size_t CalcKMerHash(Iterator it, size_t k)
			{
				size_t base = 1;
				size_t hash = 0;
				std::advance(it, k - 1);
				for(size_t i = 0; i < k; i++)
				{			
					hash += *it * base;
					base *= HASH_BASE;
					if(i != k - 1)
					{
						--it;
					}
				}		

				return hash;
			}

		private:
			size_t k_;
			size_t highPow_;
			Iterator kMerStart_;
			Iterator kMerEnd_;
			Iterator end_;
			size_t value_;
		};

	template<class T>
		const size_t SlidingWindow<T>::HASH_BASE = 57;
	
	class KMerHashFunction
	{
	public:
		size_t operator ()(DNASequence::StrandIterator it) const
		{
			return SlidingWindow<DNASequence::StrandIterator>::CalcKMerHash(it, k_);
		}

		KMerHashFunction(size_t k): k_(k) {}
	private:
		size_t k_;
	};		
		
	class KMerEqualTo
	{
	public:
		bool operator()(DNASequence::StrandIterator it1, DNASequence::StrandIterator it2) const
		{	
			DNASequence::StrandIterator end1(AdvanceForward(it1, k_));
			return std::mismatch(it1, end1, it2).first == end1;
		}

		KMerEqualTo(size_t k): k_(k) {}
	private:
		size_t k_;
	};

	class WindowHashFunction
	{
	public:
		WindowHashFunction(const SlidingWindow<DNASequence::StrandIterator> & window): window_(window) {}

		size_t operator () (DNASequence::StrandIterator it) const
		{
			if(window_.GetBegin() == it)
			{
				return window_.GetValue();
			}

			return SlidingWindow<DNASequence::StrandIterator>::CalcKMerHash(it, window_.GetK());
		}

	private:
		const SlidingWindow<DNASequence::StrandIterator> & window_;
	};
}

#endif