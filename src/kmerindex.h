#ifndef _KMER_INDEX_H_
#define _KMER_INDEX_H_

#include "dnasequence.h"
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
				return it.GetHashCode(k_);
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
				if(it1.GetHashCode(k_) != it2.GetHashCode(k_))
				{
					return false;
				}

				StrandIterator end1(AdvanceForward(it1, k_));
				return std::mismatch(it1, end1, it2).first == end1;
			}

			KMerEqualTo(size_t k): k_(k) {}
		private:
			size_t k_;
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
		
		typedef IndexMultiSet<StrandIterator, IndexTransformer, KMerHashFunction, KMerEqualTo> KMerMultiSet;		
		//Container in which we store all kmers
		size_t k_;
		KMerMultiSet * kmer_;
		const DNASequence * sequence_;

		void IndexKMers(StrandIterator start, StrandIterator end);
	};

}

#endif