#include "graphalgorithm.h"

namespace SyntenyBuilder
{
	namespace
	{		
		bool IsPositive(StrandIterator it)
		{
			return it.GetDirection() == DNASequence::positive;
		}

		size_t Extend(const DNASequence & sequence, size_t k, std::vector<StrandIterator> kmer)
		{
			size_t ret = k;
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

 			return ret - 1;
		}

		void FindMultiEdges(const DNASequence & sequence,
			const BifurcationStorage & bifStorage,
			size_t k,
			std::vector<std::vector<StrandIterator> > & ret)
		{
			ret.clear();
			std::vector<StrandIterator> kmer;
			std::vector<StrandIterator> edge[4];
			for(size_t i = 0; i < bifStorage.GetMaxId(); i++)
			{
				kmer.clear();
				bifStorage.ListPositions(i, std::back_inserter(kmer));
				for(size_t j = 0; j < kmer.size(); j++)
				{
					StrandIterator it = AdvanceForward(kmer[j], k);
					if(Valid(it, sequence))
					{
						edge[DNASequence::alphabet.find(*it)].push_back(kmer[j]);
					}
				}

				for(size_t j = 0; j < 4; j++)
				{
					if(edge[j].size() > 1 && std::find_if(edge[j].begin(), edge[j].end(), IsPositive) != edge[j].end())
					{
						ret.push_back(edge[j]);
					}

					edge[j].clear();
				}
			}
		}

	}

	void GraphAlgorithm::ListNonBranchingPaths(const DNASequence & sequence,
		const BifurcationStorage & bifStorage,
		size_t k,
		std::ostream & out,
		std::ostream & indexOut)
	{
		size_t count = 0;
		std::vector<std::vector<StrandIterator> > edge;
		FindMultiEdges(sequence, bifStorage, k, edge);
		std::vector<std::pair<size_t, size_t> > multiEdge;
		for(size_t i = 0; i < edge.size(); i++)
		{
			multiEdge.push_back(std::make_pair(edge[i].size(), i));
		}

		std::string buf;
		std::sort(multiEdge.begin(), multiEdge.end(), std::greater<std::pair<size_t, size_t> >());
		for(size_t i = 0; i < multiEdge.size(); i++)
		{
			bool repeat = false;
			size_t edgeId = multiEdge[i].second;
			size_t forward = Extend(sequence, k, edge[edgeId]);
			if(edge[edgeId].size())
			{
				out << "Consensus: " << std::endl;
				StrandIterator start = edge[edgeId][0];
				StrandIterator end = AdvanceForward(edge[edgeId][0], forward);				
				std::copy(start, end, std::ostream_iterator<char>(out));
				out << std::endl;	

				for(size_t j = 0; j < edge[edgeId].size(); j++)
				{
					buf.clear();
					start = edge[edgeId][j];
					end = AdvanceForward(edge[edgeId][j], forward);
					std::pair<size_t, size_t> coord = sequence.SpellOriginal(start, end, std::back_inserter(buf));
					out << (edge[edgeId][j].GetDirection() == DNASequence::positive ? '+' : '-') << "s, ";
					out << coord.first << ':' << coord.second << " " << buf << std::endl;
					indexOut << coord.second - coord.first << ' ' << coord.first << ' ' << coord.second << std::endl;
				}

				indexOut << DELIMITER << std::endl;
			}			
		}
	}	
}
