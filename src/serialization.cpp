#include "graphalgorithm.h"

namespace SyntenyBuilder
{
	namespace
	{	
		void OutputEdge(size_t k, StrandIterator it, size_t pos, std::ostream & out)
		{
			CopyN(it, k, std::ostream_iterator<char>(out));
			out << " -> ";
			CopyN(++StrandIterator(it), k, std::ostream_iterator<char>(out));

			char buf[1 << 8];
			if(it.GetDirection() == DNASequence::positive)
			{
				sprintf(&buf[0], "[color=\"%s\", label=\"%lu\"];", "blue", static_cast<long long unsigned>(pos));
			}
			else
			{
				sprintf(&buf[0], "[color=\"%s\", label=\"%lu\"];", "red", static_cast<long long unsigned>(pos));
			}

			out << " " << buf << std::endl;
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

				if(start != end)
				{
					bool foo = start == end;
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
		for(size_t i = 0; i < s.Size(); i++)
		{
			out << i % 10;
		}
		
		out << std::endl;
		std::copy(s.PositiveBegin(), s.PositiveEnd(), std::ostream_iterator<char>(out));
		out << std::endl;
		std::copy(s.NegativeBegin(), s.NegativeEnd(), std::back_inserter(rcomp));
		std::copy(rcomp.rbegin(), rcomp.rend(), std::ostream_iterator<char>(out));
		out << std::endl;
	}

	void GraphAlgorithm::PrintPath(StrandIterator e, size_t k, size_t distance, std::ostream & out)
	{
		out << (e.GetDirection() == DNASequence::positive ? "s+ " : "s- ");
		CopyN(e, distance + k, std::ostream_iterator<char>(out));
		std::cerr << std::endl;
	}			
	
	void GraphAlgorithm::SerializeCondensedGraph(const DNASequence & sequence,
		const BifurcationStorage & bifStorage, size_t k, std::ostream & out)
	{
		out << "digraph G" << std::endl << "{" << std::endl;
		out << "rankdir=LR" << std::endl;
		SerializeCondensed(bifStorage, sequence.PositiveBegin(), sequence.PositiveEnd(), out);
		SerializeCondensed(bifStorage, sequence.NegativeBegin(), sequence.NegativeEnd(), out);
		out << "}" << std::endl;
	}

	void GraphAlgorithm::SerializeGraph(const DNASequence & sequence, size_t k, std::ostream & out)
	{
		out << "digraph G" << std::endl << "{" << std::endl;
		out << "rankdir=LR" << std::endl;
		StrandIterator jt = AdvanceForward(sequence.PositiveBegin(), sequence.PositiveEnd(), k); 
		size_t pos = 1;
		for(StrandIterator it = sequence.PositiveBegin(); jt != sequence.PositiveEnd(); ++it, ++jt, ++pos)
		{
			OutputEdge(k, it, pos, out);
		}
		
		pos = 1;
		jt = AdvanceForward(sequence.NegativeBegin(), sequence.NegativeEnd(), k); 
		for(StrandIterator it = sequence.NegativeBegin(); jt != sequence.NegativeEnd(); ++it, ++jt, ++pos)
		{
			OutputEdge(k, it, sequence.Size() - pos + 1, out);
		}

		out << "}";
	}
}
