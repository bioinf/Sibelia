//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "outputgenerator.h"

namespace SyntenyFinder
{
	namespace
	{
		const char COVERED = 1;
		struct FooIt: public std::iterator<std::forward_iterator_tag, char>
		{
			FooIt& operator++()
			{
				return *this;
			}

			FooIt operator++(int)
			{
				return *this;
			}

			char operator * ()
			{
				return COVERED;
			}
		};

		typedef std::pair<size_t, std::vector<BlockInstance> > GroupedBlock;
		typedef std::vector<GroupedBlock> GroupedBlockList;

		bool ByFirstElement(const GroupedBlock & a, const GroupedBlock & b)
		{
			return a.first < b.first;
		}

		std::string OutputIndex(const BlockInstance & block)
		{
			std::stringstream out;			
			out << block.GetChr() + 1 << '\t' << (block.GetSignedBlockId() < 0 ? '-' : '+') << '\t';
			out << block.GetStart() << '\t' << block.GetEnd() << '\t' << block.GetEnd() - block.GetStart();
			return out.str();
		}

		template<class Iterator>
			void OutputLines(Iterator start, size_t length, std::ostream & out)
			{
				for(size_t i = 1; i <= length; i++, ++start)
				{
					out << *start;
					if(i % 80 == 0 && i != length)
					{
						out << std::endl;
					}
				}
			}

		std::vector<double> CalculateCoverage(const ChrList & chrList, GroupedBlockList::const_iterator start, GroupedBlockList::const_iterator end)
		{
			std::vector<double> ret;
			std::vector<char> cover;
			double totalBp = 0;
			double totalCoveredBp = 0;
			for(size_t chr = 0; chr < chrList.size(); chr++)
			{
				totalBp += chrList[chr].sequence.size();
				cover.assign(chrList[chr].sequence.size(), 0);
				for(GroupedBlockList::const_iterator it = start; it != end; ++it)
				{
					for(size_t i = 0; i < it->second.size(); i++)
					{
						if(it->second[i].GetChr() == chr)
						{
							CopyN(FooIt(), it->second[i].GetLength(), cover.begin() + it->second[i].GetStart());
						}
					}
				}

				double nowCoveredBp = static_cast<double>(std::count(cover.begin(), cover.end(), COVERED));
				ret.push_back(nowCoveredBp / cover.size() * 100);
				totalCoveredBp += nowCoveredBp;
			}
			
			ret.insert(ret.begin(), totalCoveredBp / totalBp * 100);
			return ret;
		}		
	}

	void OutputGenerator::ListChrs(std::ostream & out) const
	{
		out << "Chr_id\tSize\tDescription" << std::endl;
		for(size_t i = 0; i < chrList_.size(); i++)
		{
			out << i + 1 << '\t' << chrList_[i].sequence.size() << '\t' << chrList_[i].description << std::endl;
		}

		out << DELIMITER << std::endl;
	}

	void OutputGenerator::GenerateReport(const std::string & fileName) const
	{
		std::ofstream out;
		TryOpenFile(fileName, out);
		GroupedBlockList sepBlock;
		std::vector<IndexPair> group;
		GroupBy(blockList_, compareById, std::back_inserter(group));
		for(std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
		{
			sepBlock.push_back(std::make_pair(it->second - it->first, std::vector<BlockInstance>(blockList_.begin() + it->first, blockList_.begin() + it->second)));
		}

		ListChrs(out);
		out << "Degree\tCount\tTotal";
		for(size_t i = 0; i < chrList_.size(); i++)
		{
			out << "\tChr " << i + 1;
		}

		out << std::endl;
		group.clear();
		GroupBy(sepBlock, ByFirstElement, std::back_inserter(group));
		group.push_back(IndexPair(0, sepBlock.size()));
		for(std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
		{			
			if(it != group.end() - 1)
			{
				out << sepBlock[it->first].first << '\t' << it->second - it->first << '\t';
			}
			else
			{
				out << "All\t" << it->second - it->first << "\t";
			}

			out.precision(2);
			out.setf(std::ostream::fixed);
			std::vector<double> coverage = CalculateCoverage(chrList_, sepBlock.begin() + it->first, sepBlock.begin() + it->second);
			std::copy(coverage.begin(), coverage.end(), std::ostream_iterator<double>(out, "%\t"));
			out << std::endl;
		}

		out << DELIMITER << std::endl;
	}

	void OutputGenerator::ListChromosomesAsPermutations(const std::string & fileName) const
	{
		std::ofstream out;
		TryOpenFile(fileName, out);
		std::vector<IndexPair> group;
		GroupBy(blockList_, compareByChr, std::back_inserter(group));
 		for(std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
		{			
			out.setf(std::ios_base::showpos);	
			size_t length = it->second - it->first;
			size_t chr = blockList_[it->first].GetChr();
			out << '>' << chrList_[chr].description << std::endl;
			std::sort(blockList_.begin() + it->first, blockList_.begin() + it->second);
			CopyN(CFancyIterator(blockList_.begin() + it->first, boost::bind(&BlockInstance::GetSignedBlockId, _1), 0), length, std::ostream_iterator<int>(out, " "));
			out << "$" << std::endl;
		}
	}

	void OutputGenerator::ListBlocksIndices(const std::string & fileName) const
	{
		std::ofstream out;
		TryOpenFile(fileName, out);
		ListChrs(out);		
		std::vector<IndexPair> group;
		GroupBy(blockList_, compareById, std::back_inserter(group));
		for(std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
		{
			size_t length = it->second - it->first;
			out << "Block #" << blockList_[it->first].GetBlockId() << std::endl;
			out << "Chr_id\tStrand\tStart\tEnd\tLength" << std::endl;
			CopyN(CFancyIterator(blockList_.begin() + it->first, OutputIndex, std::string()), length, std::ostream_iterator<std::string>(out, "\n"));
			out << DELIMITER << std::endl;
		}
	}

	void OutputGenerator::ListBlocksSequences(const std::string & fileName) const
	{
		std::ofstream out;
		TryOpenFile(fileName, out);
		std::vector<IndexPair> group;
		GroupBy(blockList_, compareById, std::back_inserter(group));
		for(std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
		{
			for(size_t block = it->first; block < it->second; block++)
			{
				size_t length = blockList_[block].GetLength();
				char strand = blockList_[block].GetSignedBlockId() > 0 ? '+' : '-';
				const std::string & chr = chrList_[blockList_[block].GetChr()].sequence;
				out << ">Seq=\"" << chrList_[blockList_[block].GetChr()].description << "\",Strand='" << strand << "',";
				out << "Block_id=" << blockList_[block].GetBlockId() << ",Start=" ;
				out << blockList_[block].GetStart() << ",End=" << blockList_[block].GetEnd() << std::endl;

				if(blockList_[block].GetSignedBlockId() > 0)
				{
					OutputLines(chr.begin() + blockList_[block].GetStart(), length, out);
				}
				else
				{
					std::string::const_reverse_iterator it(chr.begin() + blockList_[block].GetEnd());
					OutputLines(CFancyIterator(it, DNASequence::Translate, ' '), length, out);
				}

				out << std::endl;
			}
		}
	}

	void OutputGenerator::GenerateCircosOutput(const std::string & outFile, const std::string & outDir) const
	{
		// TODO: create directory outDir

		//copy template file
		std::ofstream out;
		TryOpenFile(outFile, out);
		out << circosTemplate_;

		//blocks must be sorted by id
		BlockList sortedBlocks = blockList_;
		std::sort(sortedBlocks.begin(), sortedBlocks.end(), compareById);

		//write link and highlights file
		int idLength = static_cast<int>(log10(static_cast<double>(sortedBlocks.size()))) + 1;
		int lastId = 0;
		int linkCount = 0;
		BlockList blocksToLink;
		std::ofstream linksFile;				
		std::ofstream highlightFile;
		TryOpenFile(outDir + "/circos.segdup.txt", linksFile);
		TryOpenFile(outDir + "/circos.highlight.txt", highlightFile);

		for(BlockList::iterator itBlock = sortedBlocks.begin(); itBlock != sortedBlocks.end(); ++itBlock)
		{
			highlightFile << "hs" << itBlock->GetChr() + 1 << " ";
			highlightFile << itBlock->GetStart() << " " << itBlock->GetEnd() << std::endl;

			if (itBlock->GetBlockId() != lastId)
			{
				blocksToLink.clear();
				lastId = itBlock->GetBlockId();
			}
			for (BlockList::iterator itPair = blocksToLink.begin(); itPair != blocksToLink.end(); ++itPair)
			{
				//link start
				linksFile << "block_" << std::setw(idLength) << std::setfill('0') << linkCount << " ";
				linksFile << "hs" << itBlock->GetChr() + 1 << " ";
				linksFile << itBlock->GetStart() << " " << itBlock->GetEnd() << std::endl;
				//link end
				linksFile << "block_" << std::setw(idLength) << std::setfill('0') << linkCount << " ";
				linksFile << "hs" << itPair->GetChr() + 1 << " ";
				linksFile << itPair->GetStart() << " " << itPair->GetEnd() << std::endl;

				++linkCount;
			}
			blocksToLink.push_back(*itBlock);
		}

		//write kariotype file
		std::ofstream karFile;
		TryOpenFile(outDir + "/circos.sequences.txt", karFile);		
		
		for (size_t i = 0; i < chrList_.size(); ++i)
		{
			karFile << "chr - hs" << i + 1 << " " << i + 1 << " 0 " << chrList_[i].sequence.length();
			karFile	<< " chr" << i + 1 << std::endl;
		}
	}

	void OutputGenerator::GenerateD3Output(const std::string & outFile) const
	{
        //open output file
        std::ofstream out;
        TryOpenFile(outFile, out);
        out << "chart_data = [" << std::endl;

        //blocks must be sorted by id
        BlockList sortedBlocks = blockList_;
        std::sort(sortedBlocks.begin(), sortedBlocks.end(), compareByStart);

        // write to output file
        int lastId = 0;
        bool first_line = true;
		for(BlockList::iterator itBlock = sortedBlocks.begin(); itBlock != sortedBlocks.end(); ++itBlock) // O(N^2) by number of blocks, can be optimized
		{
            if (!first_line)
                out << ",";
            else
                first_line = false;
            out << "    {";
            out << "\"name\":\"chr" << itBlock->GetChr() + 1 << "." << std::setfill('0') << std::setw(8) << itBlock->GetStart() << "-" << std::setfill('0') << std::setw(10) << itBlock->GetEnd() << "\",";
            out << "\"size\":" << (itBlock->GetEnd() - itBlock->GetStart()) << ",";
            out << "\"imports\":[";
            bool first = true;
            for (BlockList::iterator itPair = sortedBlocks.begin(); itPair != sortedBlocks.end(); ++itPair)
            {
                if (itPair->GetBlockId() == itBlock->GetBlockId() && itPair != itBlock)
                {
                    if (!first)
                        out << ",";
                    else
                        first = false;
                    out << "\"chr" << itPair->GetChr() + 1 << "." << std::setfill('0') << std::setw(8) << itPair->GetStart() << "-" << std::setfill('0') << std::setw(10) << itPair->GetEnd() << "\"";
                }
			}
            out << "]";
            out << "}" << std::endl;
		}
        out << "];" << std::endl;
	}

	void OutputGenerator::TryOpenFile(const std::string & fileName, std::ofstream & stream) const
	{
		stream.open(fileName.c_str());
		if(!stream)
		{
			throw std::runtime_error(("Cannot open file " + fileName).c_str());
		}
	}

	void OutputGenerator::OutputBuffer(const std::string & fileName, const std::string & buffer) const
	{
		std::ofstream out;
		TryOpenFile(fileName, out);
		out << buffer;
	}
}
