#include "graphalgorithm.h"
namespace SyntenyBuilder
{
	namespace
	{	
		typedef boost::unordered_set<StrandIterator, KMerIndex::KMerHashFunction, KMerIndex::KMerEqualTo> KMerSimpleSet;

		void OutputEdge(const KMerIndex & index, StrandIterator it, std::ostream & out)
		{
			CopyN(it, index.GetK(), std::ostream_iterator<char>(out));
			out << " -> ";
			CopyN(++StrandIterator(it), index.GetK(), std::ostream_iterator<char>(out));

			char buf[1 << 8];
			if(it.GetDirection() == DNASequence::positive)
			{
				sprintf(&buf[0], "[color=\"%s\", label=\"%lu\"];", "blue", static_cast<long long unsigned>(it.GetPosition()));
			}
			else
			{
				sprintf(&buf[0], "[color=\"%s\", label=\"%lu\"];", "red", static_cast<long long unsigned>(it.GetPosition()));
			}

			out << " " << buf;
		}
		
		void ProcessIterator(KMerSimpleSet & visit, const KMerIndex & index, DNASequence::StrandIterator it, std::ostream & out)
		{
			std::vector<StrandIterator> kmer;		
			if(it.ProperKMer(index.GetK()) && visit.find(it) == visit.end())
			{
				visit.insert(it);				
				index.ListEquivalentKmers(it, kmer);
				for(size_t i = 0; i < kmer.size(); i++)
				{
					if(kmer[i].ProperKMer(index.GetK() + 1))
					{
						OutputEdge(index, kmer[i], out);
						out << std::endl;
					}
				}
			}
		}

		void SerializeCondensed(const BifurcationStorage & bifStorage,
			StrandIterator start, StrandIterator end, std::ostream & out)
		{
			for(; start != end && bifStorage.GetBifurcation(start) == BifurcationStorage::NO_BIFURCATION; ++start);
			size_t prev = start != end ? bifStorage.GetBifurcation(start) : -1;
			for(; start != end; )
			{
				size_t step = 0;
				for(++start; start != end && bifStorage.GetBifurcation(start) == BifurcationStorage::NO_BIFURCATION; 
					++start, ++step);
				{
					size_t bifId = bifStorage.GetBifurcation(start);
					if(bifId != BifurcationStorage::NO_BIFURCATION)
					{
						char buf[1 << 8];
						out << prev << " -> " << bifId;
						if(start.GetDirection() == DNASequence::positive)
						{
							sprintf(&buf[0], "[color=\"%s\", label=\"%lu\"];", "blue", static_cast<long long unsigned>(step + 1));
						}
						else
						{
							sprintf(&buf[0], "[color=\"%s\", label=\"%lu\"];", "red", static_cast<long long unsigned>(step + 1));
						}

						out << " " << buf << std::endl;
					}

					prev = bifId;
				}
			}
		}
	}

	void GraphAlgorithm::PrintRaw(const DNASequence & s, std::ostream & out)
	{
		std::string rcomp;
		s.SpellRaw(std::ostream_iterator<char>(out));
		out << std::endl;
		for(size_t i = 0; i < s.Size(); i++)
		{
			out << i % 10;
		}
		
		out << std::endl;
		std::copy(s.PositiveBegin(), s.PositiveRightEnd(), std::ostream_iterator<char>(out));
		out << std::endl;
		std::copy(s.NegativeBegin(), s.NegativeRightEnd(), std::back_inserter(rcomp));
		std::copy(rcomp.rbegin(), rcomp.rend(), std::ostream_iterator<char>(out));
		out << std::endl;
	}

	void GraphAlgorithm::PrintPath(StrandIterator e, size_t k, size_t distance, std::ostream & out)
	{
		out << e.GetPosition() << " ";
		out << (e.GetDirection() == DNASequence::positive ? "s+ " : "s- ");
		CopyN(e, distance + k, std::ostream_iterator<char>(out));
		std::cerr << std::endl;
	}			

	void GraphAlgorithm::SerializeCondensedGraph(const DNASequence & sequence, size_t k, std::ostream & out)
	{
		BifurcationStorage bifStorage;
		size_t bifurcationCount = GraphAlgorithm::EnumerateBifurcations(sequence, k, bifStorage);
		std::cerr << "Total bifurcations: " << bifurcationCount << std::endl;
		out << "digraph G" << std::endl << "{" << std::endl;
		out << "rankdir=LR" << std::endl;
		SerializeCondensed(bifStorage, sequence.PositiveBegin(), sequence.PositiveRightEnd(), out);
		SerializeCondensed(bifStorage, sequence.NegativeBegin(), sequence.NegativeRightEnd(), out);

		out << "}" << std::endl;
	}

	void GraphAlgorithm::SerializeGraph(const DNASequence & sequence, size_t k, std::ostream & out)
	{
		out << "digraph G" << std::endl << "{" << std::endl;
		out << "rankdir=LR" << std::endl;
		KMerIndex index(&sequence);
		index.SetupIndex(k);
		KMerSimpleSet visit(sequence.Size(), KMerIndex::KMerHashFunction(index.GetK()),
			KMerIndex::KMerEqualTo(index.GetK()));
		for(StrandIterator it = sequence.PositiveBegin(); it.Valid(); ++it)
		{
			ProcessIterator(visit, index, it, out);
		}
		
		for(StrandIterator it = sequence.NegativeBegin(); it.Valid(); ++it)
		{
			ProcessIterator(visit, index, it, out);
		}

		out << "}";
	}
}
