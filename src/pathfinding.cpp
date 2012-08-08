#include "graphalgorithm.h"

namespace SyntenyBuilder
{
	namespace
	{
		typedef boost::unordered_set<size_t> ElementVisitSet;
		bool Invalid(const ElementVisitSet & visit, const StrandIterator & g)
		{
			return visit.find(g.GetElementId()) != visit.end();
		}

		void Invalidate(ElementVisitSet & visit, const StrandIterator & g)
		{
			visit.insert(g.GetElementId());
		}

		size_t Extend(ElementVisitSet & visit,
			std::vector<StrandIterator> kmer,
			boost::function<void (StrandIterator&)> extender,
			boost::function<void (const StrandIterator&)> invalidator)
		{
			size_t ret = 0;			
		/*	std::for_each(kmer.begin(), kmer.end(), extender);			
			for(bool fail = false; !fail; ret++)
			{
				fail = !kmer[0].Valid() || visit[kmer[0].GetPosition()];
				if(!fail)
				{
					char consensus = *kmer[0];
					for(size_t i = 1; i < kmer.size() && !fail; i++)
					{
						fail = !kmer[i].Valid() || visit[kmer[i].GetPosition()] || *kmer[i] != consensus;
					}
				}

				if(!fail)
				{
					std::for_each(kmer.begin(), kmer.end(), invalidator);
					std::for_each(kmer.begin(), kmer.end(), extender);
				}
			}

 			return ret - 1;*/
		}
	}

	void GraphAlgorithm::ListNonBranchingPaths(const DNASequence & sequence,
		BifurcationStorage & bifStorage,
		size_t k,
		std::ostream & out,
		std::ostream & indexOut)
	{
		size_t count = 0;
		std::vector<std::pair<size_t, size_t> > multiKmer;
		for(size_t i = 0; i < bifStorage.GetMaxId(); i++)
		{
			size_t count = bifStorage.CountBifurcations(i);
			multiKmer.push_back(std::make_pair(count, i));
		}
		
		std::string buf;
		std::vector<StrandIterator> kmer;
		std::sort(multiKmer.begin(), multiKmer.end(), std::greater<std::pair<size_t, size_t> >());	
		ElementVisitSet visit;

		boost::function<StrandIterator& (StrandIterator&)> moveForward = boost::bind(&StrandIterator::operator++, _1);
		boost::function<StrandIterator& (StrandIterator&)> moveBackward = boost::bind(&StrandIterator::operator--, _1);		
		boost::function<bool (const StrandIterator&)> invalid = boost::bind(Invalid, boost::ref(visit), _1);
		boost::function<void (const StrandIterator&)> invalidator = boost::bind(Invalidate, boost::ref(visit), _1);

		for(size_t i = 0; i < multiKmer.size(); i++)
		{
			bifStorage.ListPositions(multiKmer[i].second, std::back_inserter(kmer), sequence);
			kmer.erase(std::remove_if(kmer.begin(), kmer.end(), invalid), kmer.end());
			if(kmer.size() > 1)
			{
				size_t forward = Extend(visit, kmer, moveForward, invalidator) + 1;
				size_t backward = Extend(visit, kmer, moveBackward, invalidator);
				if(forward + backward < k)
				{
					continue;
				}

				std::for_each(kmer.begin(), kmer.end(), invalidator);
				out << "Consensus: " << std::endl;
				StrandIterator end = AdvanceForward(kmer[0], forward);
				StrandIterator start = AdvanceBackward(kmer[0], backward);
				std::copy(start, end, std::ostream_iterator<char>(out));
				out << std::endl;	

				for(size_t j = 0; j < kmer.size(); j++)
				{
					buf.clear();
					end = AdvanceForward(kmer[j], forward);
					start = AdvanceBackward(kmer[j], backward);
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
