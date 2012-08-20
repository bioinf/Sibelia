#include "graphalgorithm.h"

namespace SyntenyBuilder
{
	namespace
	{	
		struct PathData
		{
		public:
			PathData() {}
			PathData(size_t bif1, size_t bif2, size_t multiplicity)
			{
				data_.push_back(bif1);
				data_.push_back(bif2);
				data_.push_back(multiplicity);
				std::sort(data_.begin(), data_.begin() + 2);
			}

			bool operator < (const PathData & toCompare) const
			{
				return data_ < toCompare.data_;
			}

		private:
			std::vector<size_t> data_;
			friend std::ostream& operator << (std::ostream & out, const PathData & data);
		};

		std::ostream& operator << (std::ostream & out, const PathData & o)
		{
			for(size_t i = 0; i < o.data_.size(); i++)
			{
				out << o.data_[i] << ' ';
			}

			return out;
		}

		bool IsPositive(StrandIterator it)
		{
			return it.GetDirection() == DNASequence::positive;
		}

		size_t Extend(const DNASequence & sequence, 
			size_t k,
			std::vector<StrandIterator> kmer,
			bool forward)
		{
			size_t ret = 0;
			boost::function<bool (const StrandIterator&)> atBegin;
			boost::function<StrandIterator& (StrandIterator&)> extender;
			if(forward)
			{
				extender = boost::bind(&StrandIterator::operator++, _1);
			}		
			else
			{				
				extender = boost::bind(&StrandIterator::operator--, _1);
				atBegin = boost::bind(AtBegin, _1, boost::cref(sequence));
			}

			bool fail = false;
			do
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
					++ret;
					if(!forward && std::find_if(kmer.begin(), kmer.end(), atBegin) != kmer.end())
					{
						fail = true;
					}
					else
					{
						std::for_each(kmer.begin(), kmer.end(), extender);
					}
				}
			}
			while(!fail);

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
						size_t pos = DNASequence::alphabet.find(*it);
						if(pos != DNASequence::alphabet.npos)
						{
							edge[pos].push_back(kmer[j]);
						}
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
		std::set<PathData> visit;
		std::sort(multiEdge.begin(), multiEdge.end(), std::greater<std::pair<size_t, size_t> >());		
		for(size_t i = 0; i < multiEdge.size(); i++)
		{
			size_t edgeId = multiEdge[i].second;
			size_t forward = Extend(sequence, k, edge[edgeId], true) + 1;
			size_t backward = Extend(sequence, k, edge[edgeId], false);
			StrandIterator start = AdvanceBackward(edge[edgeId][0], backward);
			StrandIterator end = AdvanceForward(edge[edgeId][0], forward - k);
			PathData nowData(bifStorage.GetBifurcation(start), bifStorage.GetBifurcation(end), edge[edgeId].size());
			if(visit.count(nowData) == 0)
			{
				visit.insert(nowData);
				std::stringstream outBuf;
				std::stringstream indexOutBuf;

				outBuf << "Consensus: " << std::endl;
				std::copy(start, AdvanceForward(start, forward), std::ostream_iterator<char>(outBuf));
				outBuf << std::endl;

				size_t count = 0;
				for(size_t j = 0; j < edge[edgeId].size(); j++)
				{
					buf.clear();
					start = AdvanceBackward(edge[edgeId][j], backward);
					end = AdvanceForward(edge[edgeId][j], forward);
					std::pair<size_t, size_t> coord = sequence.SpellOriginal(start, end, std::back_inserter(buf));
					if(buf.size() >= k)
					{
						count++;
						outBuf << (edge[edgeId][j].GetDirection() == DNASequence::positive ? '+' : '-') << "s, ";
						outBuf << coord.first << ':' << coord.second << " " << buf << std::endl;
						indexOutBuf << coord.second - coord.first << ' ' << coord.first << ' ' << coord.second << std::endl;
					}					
				}

				if(count > 1)
				{
					out << outBuf.str();
					indexOut << indexOutBuf.str() << DELIMITER << std::endl;
				}
			}
		}			
	}	
}
