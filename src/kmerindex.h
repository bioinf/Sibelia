#ifndef _KMER_INDEX_H_
#define _KMER_INDEX_H_

#include "hashing.h"
#include "indexmultiset.h"

namespace SyntenyBuilder
{
	class KMerIndex
	{
	public:				
		typedef DNASequence::StrandIterator StrandIterator;
		void SetupIndex(size_t k);		
		size_t GetK() const;
		size_t CountEquivalentKMers(StrandIterator it) const;
		size_t ListEquivalentKmers(StrandIterator it, std::vector<StrandIterator> & ret) const;
		KMerIndex(const DNASequence * sequence);				

		class KMerHashFunction
		{
		public:
			size_t operator ()(StrandIterator it) const
			{
				return SlidingWindow<StrandIterator>::CalcKMerHash(it, k_);
			}

			KMerHashFunction(size_t k): k_(k) {}
		private:
			size_t k_;
		};		
		
		class KMerEqualTo
		{
		public:
			bool operator()(StrandIterator it1, StrandIterator it2) const
			{	
				StrandIterator end1(AdvanceForward(it1, k_));
				return std::mismatch(it1, end1, it2).first == end1;
			}

			KMerEqualTo(size_t k): k_(k) {}
		private:
			size_t k_;
		};

		class WindowHashFunction
		{
		public:
			WindowHashFunction(const SlidingWindow<StrandIterator> & window): window_(window) {}

			size_t operator () (StrandIterator it) const
			{
				if(window_.GetBegin() == it)
				{
					return window_.GetValue();
				}

				return SlidingWindow<StrandIterator>::CalcKMerHash(it, window_.GetK());
			}

		private:
			const SlidingWindow<StrandIterator> & window_;
		};

	private:
		DISALLOW_COPY_AND_ASSIGN(KMerIndex);				
		
		class IndexTransformer
		{
		public:
			StrandIterator operator()(size_t index) const
			{
				return sequence_->PositiveByIndex(index);
			}

			IndexTransformer(const DNASequence * sequence): sequence_(sequence){}
		private:
			const DNASequence * sequence_;
		};		
		
		typedef IndexMultiSet<StrandIterator, IndexTransformer, WindowHashFunction, KMerEqualTo> KMerMultiSet;		
		//Container in which we store all kmers
		size_t k_;
		KMerMultiSet * kmer_;
		const DNASequence * sequence_;
		SlidingWindow<StrandIterator> window_;

		void IndexKMers(StrandIterator start, StrandIterator end);
	};

}

#endif