#include "graphalgorithm.h"

namespace SyntenyBuilder
{
	namespace
	{
		bool Valid(StrandIterator it, const DNASequence & sequence)
		{
			return it != sequence.PositiveEnd() && it != sequence.NegativeEnd();
		}

		size_t Extend(const DNASequence & sequence, size_t k, std::vector<StrandIterator> kmer)
		{
			size_t ret = 1;
			boost::function<StrandIterator& (StrandIterator&)> extender = boost::bind(&StrandIterator::operator++, _1);
			for(size_t i = 0; i < k; i++)
			{
				std::for_each(kmer.begin(), kmer.end(), extender);			
			}

			for(bool fail = false; !fail; ret++)
			{
				fail = !Valid(kmer[0], sequence);
				if(!fail)
				{
					char consensus = *kmer[0];
					for(size_t i = 1; i < kmer.size() && !fail; i++)
					{
						fail = !Valid(kmer[i], sequence) || *kmer[i] != consensus;
					}
				}

				if(!fail)
				{
					std::for_each(kmer.begin(), kmer.end(), extender);
				}
			}

 			return ret;
		}
	}

	void GraphAlgorithm::ListNonBranchingPaths(const DNASequence & sequence,
		const BifurcationStorage & bifStorage,
		size_t k,
		std::ostream & out,
		std::ostream & indexOut)
	{
		size_t count = 0;
		std::vector<std::pair<size_t, size_t> > multiKmer;
		for(size_t i = 0; i < bifStorage.GetMaxId(); i++)
		{
			size_t count = bifStorage.CountBifurcations(i);
			if(count > 1)
			{
				multiKmer.push_back(std::make_pair(count, i));
			}
		}
		
		std::string buf;
		std::vector<StrandIterator> kmer;
		std::sort(multiKmer.begin(), multiKmer.end(), std::greater<std::pair<size_t, size_t> >());	
		for(size_t i = 0; i < multiKmer.size(); i++)
		{
			bifStorage.ListPositions(multiKmer[i].second, std::back_inserter(kmer));
			if(kmer.size() > 1)
			{
				size_t forward = Extend(sequence, k, kmer);
				if(forward < 1)
				{
					continue;
				}

				out << "Consensus: " << std::endl;
				StrandIterator start = kmer[0];
				StrandIterator end = AdvanceForward(kmer[0], forward);				
				std::copy(start, end, std::ostream_iterator<char>(out));
				out << std::endl;	

				for(size_t j = 0; j < kmer.size(); j++)
				{
					buf.clear();
					start = kmer[j];
					end = AdvanceForward(kmer[j], forward);
					std::pair<size_t, size_t> coord = sequence.SpellOriginal(start, end, std::back_inserter(buf));
					out << (kmer[j].GetDirection() == DNASequence::positive ? '+' : '-') << "s, ";
					out << coord.first << ':' << coord.second << " " << buf << std::endl;
					indexOut << coord.second - coord.first << ' ' << coord.first << ' ' << coord.second << std::endl;
				}

				indexOut << DELIMITER << std::endl;
			}
		}
	}	
}
