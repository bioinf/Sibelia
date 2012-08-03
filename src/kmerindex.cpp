#include "kmerindex.h"

namespace SyntenyBuilder
{
	KMerIndex::KMerIndex(const DNASequence * sequence): kmer_(0), sequence_(sequence)
	{
	}

	KMerIndex::~KMerIndex()
	{
		delete kmer_;
	}

	void KMerIndex::SetupIndex(size_t k)
	{
		k_ = k;
		delete kmer_;
		window_ = SlidingWindow<StrandIterator>(sequence_->PositiveBegin(), 
			sequence_->PositiveRightEnd(), k);

		kmer_ = new KMerMultiSet(sequence_->Size(),
			IndexTransformer(sequence_), 
			WindowHashFunction(window_),
			KMerEqualTo(k));
		
		const size_t MOD = 1000000;
		for(StrandIterator it = sequence_->PositiveBegin(); it != sequence_->PositiveRightEnd(); it++)
		{
			if(it.GetPosition() % MOD == 0)
			{
				std::cerr << "Pos = " << it.GetPosition() << std::endl;
			}

			if(it.ProperKMer(k_))
			{
				kmer_->Insert(it.GetPosition());
				window_.Move();
			}
		}
	}

	void KMerIndex::IndexKMers(StrandIterator start, StrandIterator end)
	{
		for(; start != end && start.ProperKMer(k_); ++start)
		{
			kmer_->Insert(start.GetPosition());
		}
	}

	size_t KMerIndex::GetK() const
	{
		return k_;
	}

	size_t KMerIndex::CountEquivalentKMers(StrandIterator it) const
	{
		size_t ret = kmer_->Count(it);
		it.Jump(k_ - 1);
		return ret + kmer_->Count(it.Invert());
	}

	size_t KMerIndex::ListEquivalentKmers(StrandIterator it, std::vector<StrandIterator> & ret) const
	{
		ret.clear();
		std::vector<size_t> positive;
		std::vector<size_t> negative;
		kmer_->Find(it, std::back_inserter(positive));
		it.Jump(k_ - 1);
		kmer_->Find(it.Invert(), std::back_inserter(negative));		
		boost::function<StrandIterator (size_t)> createPos =
			boost::bind(&DNASequence::PositiveByIndex, sequence_, _1);

		std::transform(negative.begin(), negative.end(), std::back_inserter(ret), createPos);
		std::for_each(ret.begin(), ret.end(), boost::bind(std::advance<StrandIterator, size_t>, _1, k_ - 1));
		std::for_each(ret.begin(), ret.end(), boost::bind(&StrandIterator::MakeInverted, _1));
		std::transform(positive.begin(), positive.end(), std::back_inserter(ret), createPos);
		return ret.size();
	}
}