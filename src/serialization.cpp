//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "blockfinder.h"

namespace SyntenyFinder
{/*
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
			sprintf(&buf[0], "[color=\"%s\", label=\"(%i, %i)\"];", color.c_str(), static_cast<int>(chr), static_cast<int>(pos));
			out << " " << buf << std::endl;
		}
	}
	
	void BlockFinder::PrintRaw(const DNASequence & sequence, std::ostream & out)
	{
		for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
		{
			out << "Sequence #" << chr << std::endl;
			std::string rcomp;
			StrandIterator it = sequence.PositiveBegin(chr);
			StrandIterator end = sequence.PositiveEnd(chr);
			for(size_t i = 0; it != end; i++, ++it)
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

	void BlockFinder::PrintPath(const DNASequence & s, StrandIterator e, size_t k, size_t distance, std::ostream & out)
	{
		out << (e.GetDirection() == DNASequence::positive ? "+" : "-") << s.GlobalIndex(e) << ' ';
		CopyN(e, distance + k, std::ostream_iterator<char>(out));
		std::cerr << std::endl;
	}

	void BlockFinder::ListEdges(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k, std::vector<Edge> & edge) const
	{
		edge.clear();
		for(size_t strand = 0; strand < 2; strand++)
		{
			for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
			{
				size_t pos = 0;				
				StrandIterator start = sequence.Begin((DNASequence::Direction)strand, chr);
				StrandIterator end = sequence.End((DNASequence::Direction)strand, chr);				
				size_t prevVertex = bifStorage.GetBifurcation(start);
				size_t length = std::distance(start, end);
				for(; start != end; )
				{
					size_t step = 1;				
					StrandIterator origin = start;
					for(++start; start != end && bifStorage.GetBifurcation(start) == BifurcationStorage::NO_BIFURCATION; ++start, ++step);
					if(start != end)
					{
						char firstChar = *AdvanceForward(origin, k);
						size_t nowVertex = bifStorage.GetBifurcation(start);
						std::pair<size_t, size_t> coord = sequence.SpellOriginal(origin, AdvanceForward(start, k));
						size_t actualPos = strand == 0 ? pos : length - (pos + step + k);
						edge.push_back(Edge(chr, start.GetDirection(), prevVertex, nowVertex, actualPos, step + k, coord.first, coord.second - coord.first, firstChar));
						prevVertex = nowVertex;
						pos += step;
					}
				}
			}
		}
	}

	void BlockFinder::SerializeCondensedGraph(size_t k, std::ostream & out, ProgressCallBack f)
	{
		IndexedSequence iseq(rawSeq_, originalPos_, k, tempDir_);
		out << "digraph G" << std::endl << "{" << std::endl;
		out << "rankdir=LR" << std::endl;
		std::vector<Edge> edge;
		ListEdges(iseq.Sequence(), iseq.BifStorage(), k, edge);
		for(size_t i = 0; i < edge.size(); i++)
		{
			char buf[1 << 8];
			std::string color = edge[i].GetDirection() == DNASequence::positive ? "blue" : "red";
			int uchr = static_cast<int>(edge[i].GetChr());
			int uorpos = static_cast<int>(edge[i].GetOriginalPosition());
			int uorlength = static_cast<int>(edge[i].GetOriginalLength());
			int upos = static_cast<int>(edge[i].GetActualPosition());
			int ulength = static_cast<int>(edge[i].GetActualLength());
			out << edge[i].GetStartVertex() << " -> " << edge[i].GetEndVertex();
			sprintf(&buf[0], "[color=\"%s\", label=\"chr=%i pos=%i len=%i orpos=%i orlen=%i  ch='%c'\"];", color.c_str(), uchr, upos, ulength, uorpos, uorlength, edge[i].GetFirstChar());
			out << " " << buf << std::endl;
		}

		out << "}" << std::endl;
	}

	void BlockFinder::SerializeGraph(size_t k, std::ostream & out)
	{
		DNASequence sequence(rawSeq_, originalPos_);
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
	}*/
}
