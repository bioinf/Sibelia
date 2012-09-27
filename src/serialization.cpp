#include "graphalgorithm.h"

namespace SyntenyBuilder
{
	namespace
	{	
		typedef unsigned long long ull;

		void OutputEdge(size_t k, StrandIterator it, size_t chr, size_t pos, std::ostream & out)
		{
			CopyN(it, k, std::ostream_iterator<char>(out));
			out << " -> ";
			CopyN(++StrandIterator(it), k, std::ostream_iterator<char>(out));
			char buf[1 << 8];
			std::string color = it.GetDirection() == DNASequence::positive ? "blue" : "red";			
			sprintf(&buf[0], "[color=\"%s\", label=\"(%lu, %lu)\"];", color.c_str(), static_cast<ull>(chr), static_cast<ull>(pos));
			out << " " << buf << std::endl;
		}
		
		void SerializeCondensed(const BifurcationStorage & bifStorage, size_t k,
			StrandIterator start, StrandIterator end, size_t chr, std::ostream & out)
		{
			ull upos = 0;
			ull uchr = static_cast<ull>(chr);
			for(; start != end && bifStorage.GetBifurcation(start) == BifurcationStorage::NO_BIFURCATION; ++start, ++upos);
			size_t prev = start != end ? bifStorage.GetBifurcation(start) : -1;			
			for(; start != end; )
			{
				ull step = 1;				
				StrandIterator origin = start;
				for(++start; start != end && bifStorage.GetBifurcation(start) == BifurcationStorage::NO_BIFURCATION; ++start, ++step);
				if(start != end)
				{
					size_t bifId = bifStorage.GetBifurcation(start);
					if(bifId != BifurcationStorage::NO_BIFURCATION)
					{
						char buf[1 << 8];
						out << prev << " -> " << bifId;
						char ch = *AdvanceForward(origin, k);
						std::string color = start.GetDirection() == DNASequence::positive ? "blue" : "red";
						sprintf(&buf[0], "[color=\"%s\", label=\"chr=%lu pos=%lu len=%lu ch=%c\"];", color.c_str(), uchr, upos, step, ch);
						out << " " << buf << std::endl;
					}

					prev = bifId;
					upos += step;
				}
			}
		}
	}
	
	void GraphAlgorithm::PrintRaw(const DNASequence & sequence, std::ostream & out)
	{
		for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
		{
			out << "Sequence #" << chr << std::endl;
			std::string rcomp;
			StrandIterator it = sequence.PositiveBegin(chr);
			for(size_t i = 0; it != sequence.PositiveEnd(chr); i++, ++it)
			{
				out << i % 10;
			}

			out << std::endl;
			std::copy(sequence.PositiveBegin(chr), sequence.PositiveEnd(chr), std::ostream_iterator<char>(out));
			out << std::endl;
			std::copy(sequence.NegativeBegin(chr), sequence.NegativeEnd(chr), std::back_inserter(rcomp));
			std::copy(rcomp.rbegin(), rcomp.rend(), std::ostream_iterator<char>(out));
			out << std::endl;
		}		
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
		for(size_t strand = 0; strand < 2; strand++)
		{
			for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
			{
				StrandIterator begin = sequence.Begin((DNASequence::Direction)strand, chr);
				StrandIterator end = sequence.End((DNASequence::Direction)strand, chr);
				SerializeCondensed(bifStorage, k, begin, end, chr, out);
			}
		}
		
		out << "}" << std::endl;
	}

	void GraphAlgorithm::SerializeGraph(const DNASequence & sequence, size_t k, std::ostream & out)
	{
		out << "digraph G" << std::endl << "{" << std::endl;
		out << "rankdir=LR" << std::endl;
		for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
		{
			StrandIterator bound[4] = 
			{
				sequence.PositiveBegin(chr),
				sequence.NegativeBegin(chr),
				sequence.PositiveEnd(chr),
				sequence.NegativeEnd(chr)
			};

			for(size_t strand = 0; strand < 2; strand++)
			{
				size_t pos = 0;
				for(SlidingWindow<StrandIterator> window(bound[strand], bound[strand + 2], k + 1); window.Valid(); ++pos, window.Move())
				{
					OutputEdge(k, window.GetBegin(), chr, pos, out);
				}
			}
		}

		out << "}" << std::endl;
	}
}
