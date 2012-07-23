#include "kmerindex.h"

namespace SyntenyBuilder
{
	namespace
	{
		void MakeRangePositive(KMerIndex::StrandIterator & start, KMerIndex::StrandIterator & end)
		{
			if(start.GetDirection() == DNASequence::negative)
			{
				KMerIndex::StrandIterator newStart = (--end).Invert();
				KMerIndex::StrandIterator newEnd = ++(start.Invert());
				start = newStart;
				end = newEnd;
			}
		}
	}

	KMerIndex::KMerIndex(DNASequence * sequence): kmer_(0), sequence_(sequence)
	{
		sequence->SetHandlers(boost::bind(&KMerIndex::Invalidate, this, _1, _2),
			boost::bind(&KMerIndex::UpdateAfterCopy, this, _1, _2),
			boost::bind(&KMerIndex::UpdateAfterDelete, this, _1, _2));
	}

	void KMerIndex::SetupIndex(size_t k)
	{
		k_ = k;
		delete kmer_;
		kmer_ = new KMerMultiSet(sequence_->Size(),
			IndexTransformer(sequence_), 
			KMerHashFunction(k),
			KMerEqualTo(k));
		for(StrandIterator it = sequence_->PositiveBegin(); it != sequence_->PositiveRightEnd(); it++)
		{
			if(it.ProperKMer(k_))
			{
				kmer_->Insert(it.GetPosition());
			}
		}
	}

	void KMerIndex::IndexKMers(StrandIterator start, StrandIterator end)
	{
		for(; start != end; ++start)
		{
			kmer_->Insert(start.GetPosition());
		}
	}

	//This must be refactored later!
	void KMerIndex::Invalidate(StrandIterator start, StrandIterator end)
	{		
		MakeRangePositive(start, end);		
		for(start = AdvanceBackward(start, sequence_->PositiveBegin(), k_ - 1); start != end; ++start)
		{
			kmer_->Erase(start.GetPosition());
		}
	}

	void KMerIndex::UpdateAfterCopy(StrandIterator start, StrandIterator end)
	{
		MakeRangePositive(start, end);
		IndexKMers(AdvanceBackward(start, sequence_->PositiveBegin(), k_ - 1), end);
		/*
		std::vector<size_t> temp;
		kmer_->DumpIndex(std::back_inserter(temp));
		std::sort(temp.begin(), temp.end());		
		std::copy(temp.begin(), temp.end(), std::ostream_iterator<size_t>(std::cerr, "\n"));
		std::vector<StrandIterator> temp2;
		kmer_->Dump(std::back_inserter(temp2));
		for(size_t i = 0; i < temp2.size(); i++)
		{
			CopyN(temp2[i], k_, std::ostream_iterator<char>(std::cerr));
			std::cerr << std::endl;
		}*/
	}

	void KMerIndex::UpdateAfterDelete(StrandIterator start, StrandIterator end)
	{
		if(start.GetDirection() == DNASequence::positive)
		{
			--start;
		}
		else
		{
			start = end.Invert();
		}

		end = start;
		++end;
		start = AdvanceBackward(start, k_ - 2);
		IndexKMers(start, end);
		/*
		std::vector<size_t> temp;
		kmer_->DumpIndex(std::back_inserter(temp));
		std::sort(temp.begin(), temp.end());		
		std::copy(temp.begin(), temp.end(), std::ostream_iterator<size_t>(std::cerr, "\n"));
		std::vector<StrandIterator> temp2;
		kmer_->Dump(std::back_inserter(temp2));
		for(size_t i = 0; i < temp2.size(); i++)
		{
			CopyN(temp2[i], k_, std::ostream_iterator<char>(std::cerr));
			std::cerr << std::endl;
		}*/
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
		return static_cast<int>(ret.size());
	}
}