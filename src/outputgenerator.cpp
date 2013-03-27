//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "outputgenerator.h"
#include "platform.h"
#include "variantcalling/variant.h"

namespace SyntenyFinder
{
	namespace
	{
		const char COVERED = 1;
		typedef std::pair<size_t, std::vector<BlockInstance> > GroupedBlock;
		typedef std::vector<GroupedBlock> GroupedBlockList;

		bool ByFirstElement(const GroupedBlock & a, const GroupedBlock & b)
		{
			return a.first < b.first;
		}

		std::string OutputIndex(const BlockInstance & block)
		{
			std::stringstream out;
			out << block.GetChrInstance().GetConventionalId() << '\t' << (block.GetSignedBlockId() < 0 ? '-' : '+') << '\t';
			out << block.GetConventionalStart() << '\t' << block.GetConventionalEnd() << '\t' << block.GetEnd() - block.GetStart();
			return out.str();
		}

		// function to sort blocks in one chromosome by starting position (it's sorted lexicographically by D3)
		std::string OutputD3BlockID(const BlockInstance & block)
		{
			std::stringstream out;
			// label should (1) have no dots, and (2) be valid DOM id/class identifier
			// (exception: it may contains spaces, I'll handle it in javascript later)
//			std::string label = block.GetChrInstance().GetDescription();
//			std::replace(label.begin(), label.end(), '|', ' ');
//			std::replace(label.begin(), label.end(), ':', ' ');
//			std::replace(label.begin(), label.end(), '.', ' ');
			out << "seq" << block.GetChrInstance().GetConventionalId() << ".";
			out << "seq " << block.GetChrInstance().GetConventionalId() << " - ";
			out << std::setfill(' ') << std::setw(8) << block.GetConventionalStart() << " - ";
			out << std::setfill(' ') << std::setw(8) << block.GetConventionalEnd();
			return out.str();
		}

		void OutputLink(std::vector<BlockInstance>::iterator block, int color, int fillLength,
						int linkId, std::ostream& stream)
		{
			size_t start = block->GetConventionalStart();
			size_t end = block->GetConventionalEnd();
			if (start > end) std::swap(start, end);

			stream << "block_" << std::setw(fillLength) << std::setfill('0') << linkId << " ";
			stream << "seq" << block->GetChrId() + 1 << " ";
			stream << start << " " << end;
			stream << " color=chr" << color << "_a2" << std::endl;
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
				totalBp += chrList[chr].GetSequence().size();
				cover.assign(chrList[chr].GetSequence().size(), 0);
				for(GroupedBlockList::const_iterator it = start; it != end; ++it)
				{
					for(size_t i = 0; i < it->second.size(); i++)
					{
						if(it->second[i].GetChrInstance().GetId() == chr)
						{
							std::fill(cover.begin() + it->second[i].GetStart(), cover.begin() + it->second[i].GetEnd(), COVERED);
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

	const int OutputGenerator::CIRCOS_MAX_COLOR = 25;
	const int OutputGenerator::CIRCOS_DEFAULT_RADIUS = 1500;
	const int OutputGenerator::CIRCOS_RESERVED_FOR_LABEL = 500;
	const int OutputGenerator::CIRCOS_HIGHLIGHT_THICKNESS = 50;

	void OutputGenerator::ListChrs(std::ostream & out) const
	{
		out << "Seq_id\tSize\tDescription" << std::endl;
		for(size_t i = 0; i < chrList_.size(); i++)
		{
			out << i + 1 << '\t' << chrList_[i].GetSequence().size() << '\t' << chrList_[i].GetDescription() << std::endl;
		}

		out << DELIMITER << std::endl;
	}	

	void OutputGenerator::GenerateReport(const BlockList & block, const std::string & fileName) const
	{
		std::ofstream out;
		TryOpenFile(fileName, out);
		GroupedBlockList sepBlock;
		std::vector<IndexPair> group;
		BlockList blockList = block;
		GroupBy(blockList, compareById, std::back_inserter(group));
		for(std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
		{
			sepBlock.push_back(std::make_pair(it->second - it->first, std::vector<BlockInstance>(blockList.begin() + it->first, blockList.begin() + it->second)));
		}

		ListChrs(out);
		out << "Degree\tCount\tTotal";
		for(size_t i = 0; i < chrList_.size(); i++)
		{
			out << "\tSeq " << i + 1;
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

	void OutputGenerator::ListChromosomesAsPermutations(const BlockList & block, const std::string & fileName) const
	{
		std::ofstream out;
		TryOpenFile(fileName, out);
		std::vector<IndexPair> group;
		BlockList blockList = block;
		GroupBy(blockList, compareByChrId, std::back_inserter(group));
 		for(std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
		{
			out.setf(std::ios_base::showpos);
			size_t length = it->second - it->first;
			size_t chr = blockList[it->first].GetChrInstance().GetId();
			out << '>' << chrList_[chr].GetDescription() << std::endl;
			std::sort(blockList.begin() + it->first, blockList.begin() + it->second);
			CopyN(CFancyIterator(blockList.begin() + it->first, boost::bind(&BlockInstance::GetSignedBlockId, _1), 0), length, std::ostream_iterator<int>(out, " "));
			out << "$" << std::endl;
		}
	}

   	void OutputGenerator::RearrangementScenario(const std::vector<std::string> & steps, const std::string & fileName) const
	{
		std::ofstream out;
		TryOpenFile(fileName, out);
		for(size_t i = 0; i < steps.size(); ++i)
		{
			out << steps[i] << std::endl;
		}
	}

	void OutputGenerator::ListBlocksIndices(const BlockList & block, const std::string & fileName) const
	{
		std::ofstream out;
		TryOpenFile(fileName, out);
		ListChrs(out);
		std::vector<IndexPair> group;
		BlockList blockList = block;
		GroupBy(blockList, compareById, std::back_inserter(group));
		for(std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
		{
			size_t length = it->second - it->first;
			std::sort(blockList.begin() + it->first, blockList.begin() + it->second, compareByChrId);
			out << "Block #" << blockList[it->first].GetBlockId() << std::endl;
			out << "Seq_id\tStrand\tStart\tEnd\tLength" << std::endl;
			CopyN(CFancyIterator(blockList.begin() + it->first, OutputIndex, std::string()), length, std::ostream_iterator<std::string>(out, "\n"));
			out << DELIMITER << std::endl;
		}
	}

	void OutputGenerator::ListBlocksSequences(const BlockList & block, const std::string & fileName) const
	{
		std::ofstream out;
		TryOpenFile(fileName, out);
		std::vector<IndexPair> group;
		BlockList blockList = block;
		GroupBy(blockList, compareById, std::back_inserter(group));
		for(std::vector<IndexPair>::iterator it = group.begin(); it != group.end(); ++it)
		{
			for(size_t block = it->first; block < it->second; block++)
			{
				size_t length = blockList[block].GetLength();
				char strand = blockList[block].GetSignedBlockId() > 0 ? '+' : '-';
				const FASTARecord & chr = blockList[block].GetChrInstance();
				out << ">Seq=\"" << chr.GetDescription() << "\",Strand='" << strand << "',";
				out << "Block_id=" << blockList[block].GetBlockId() << ",Start=" ;
				out << blockList[block].GetConventionalStart() << ",End=" << blockList[block].GetConventionalEnd() << std::endl;

				if(blockList[block].GetSignedBlockId() > 0)
				{
					OutputLines(chr.GetSequence().begin() + blockList[block].GetStart(), length, out);
				}
				else
				{
					std::string::const_reverse_iterator it(chr.GetSequence().begin() + blockList[block].GetEnd());
					OutputLines(CFancyIterator(it, DNASequence::Translate, ' '), length, out);
				}

				out << std::endl;
			}
		}
	}

	void OutputGenerator::WriteCircosImageConfig(const std::string & outDir, const std::string & fileName, int r) const
	{
		std::ofstream imageConfig;
		TryOpenFile(outDir + "/" + fileName, imageConfig);
		imageConfig << circosImageConfig;
		imageConfig << "radius = " << r << "p" << std::endl;
	}

	void OutputGenerator::GenerateHierarchyCircosOutput(const std::vector<BlockList> & history, const std::string & outFile, const std::string & outDir) const
	{
		int r = 100;
		std::ofstream config;
		CreateDirectory(outDir);
		TryOpenFile(outFile, config);
		config << circosTemplate;		
		WriteCircosLinks(outDir, "circos.segdup.txt", history.back());
		WriteCircosKaryoType(outDir, "circos.sequences.txt");
		config << "<highlights>\n\tfill_color = green" << std::endl;		
		WriteCircosHighlight(outDir, "circos.highlight.txt", history.back(), 0, 0, true, config);		
		for(std::vector<BlockList>::const_reverse_iterator it = ++history.rbegin(); it != history.rend(); ++it)
		{			
			std::stringstream ss;
			ss << "circos.highlight" << it - history.rbegin() << ".txt";
			WriteCircosHighlight(outDir, ss.str(), *it, r, r + CIRCOS_HIGHLIGHT_THICKNESS, false, config);			
			r += static_cast<int>(CIRCOS_HIGHLIGHT_THICKNESS * 1.5);
		}

		config << "</highlights>" << std::endl;	
		std::stringstream ss;
		ss << "<ideogram>\n\tlabel_radius = 1r + " << r << "p\n</ideogram>" << std::endl;
		config << ss.str();
		WriteCircosImageConfig(outDir, "circos.image.conf", CIRCOS_DEFAULT_RADIUS + CIRCOS_RESERVED_FOR_LABEL + r);
	}

	void OutputGenerator::GenerateCircosOutput(const BlockList & blockList, const std::string & outFile, const std::string & outDir) const
	{		
		std::ofstream config;
		CreateDirectory(outDir);
		TryOpenFile(outFile, config);
		config << circosTemplate;		
		WriteCircosLinks(outDir, "circos.segdup.txt", blockList);
		WriteCircosKaryoType(outDir, "circos.sequences.txt");
		config << "<highlights>\n\tfill_color = green" << std::endl;		
		WriteCircosHighlight(outDir, "circos.highlight.txt", blockList, 0, 0, true, config);
		config << "</highlights>" << std::endl;
		config << "<ideogram>\n\tlabel_radius = 1.08r\n</ideogram>" << std::endl;
		WriteCircosImageConfig(outDir, "circos.image.conf", CIRCOS_DEFAULT_RADIUS);
	}

	void OutputGenerator::WriteCircosLinks(const std::string & outDir, const std::string & fileName, const BlockList & block) const
	{
		//blocks must be sorted by id
		BlockList sortedBlocks = block;
		std::sort(sortedBlocks.begin(), sortedBlocks.end(), compareById);

		//write link and highlights file
		int idLength = static_cast<int>(log10(static_cast<double>(sortedBlocks.size()))) + 1;
		int lastId = 0;
		int linkCount = 0;
		BlockList blocksToLink;
		std::ofstream linksFile;
		TryOpenFile(outDir + "/" + fileName, linksFile);

		int color = 0;
		for(BlockList::iterator itBlock = sortedBlocks.begin(); itBlock != sortedBlocks.end(); ++itBlock)
		{
			if (itBlock->GetBlockId() != lastId)
			{
				blocksToLink.clear();
				lastId = itBlock->GetBlockId();
			}

			for (BlockList::iterator itPair = blocksToLink.begin(); itPair != blocksToLink.end(); ++itPair)
			{
				color = (color + 1) % CIRCOS_MAX_COLOR;
				//link start
				OutputLink(itBlock, color, idLength, linkCount, linksFile);
				//link end
				OutputLink(itPair, color, idLength, linkCount, linksFile);
				++linkCount;
			}

			blocksToLink.push_back(*itBlock);
		}
	}

	void OutputGenerator::WriteCircosHighlight(const std::string & outDir, const std::string & fileName, const BlockList & block, int r0, int r1, bool ideogram, std::ofstream & config) const
	{
		int color = 0;
		BlockList sortedBlocks = block;
		std::sort(sortedBlocks.begin(), sortedBlocks.end(), compareById);
		BlockList blocksToLink;		
		std::ofstream highlightFile;		
		TryOpenFile(outDir + "/" + fileName, highlightFile);
		for(BlockList::iterator itBlock = sortedBlocks.begin(); itBlock != sortedBlocks.end(); ++itBlock)
		{
			highlightFile << "seq" << itBlock->GetChrInstance().GetConventionalId() << " ";
			size_t blockStart = itBlock->GetConventionalStart();
			size_t blockEnd = itBlock->GetConventionalEnd();
			if (blockStart > blockEnd)
			{
				std::swap(blockStart, blockEnd);
			}

			if(itBlock != sortedBlocks.begin() && itBlock->GetBlockId() != (itBlock - 1)->GetBlockId())
			{
				color = (color + 1) % CIRCOS_MAX_COLOR;
			}

			highlightFile << blockStart << " " << blockEnd;
			if(!ideogram)
			{
				highlightFile << " fill_color=chr" << color << "_a0";
			}

			highlightFile << std::endl;
		}

		std::string prefix = "\t\t";
		config << "\t<highlight>" << std::endl;
		config << prefix << "file = " << fileName << std::endl;
		config << prefix << "ideogram = " << (ideogram ? "yes" : "no") << std::endl;
		config << prefix << "fill_color = blue_a3" << std::endl;
		config << prefix << "stroke_color = black" << std::endl;
		config << prefix << "stroke_thickness = 4" << std::endl;
		if(!ideogram)
		{
			config << prefix << "r0 = 1r +" << r0 << "p" << std::endl;
			config << prefix << "r1 = 1r +" << r1 << "p" << std::endl;
		}

		config << "\t</highlight>" << std::endl;
	}

	void OutputGenerator::WriteCircosKaryoType(const std::string & outDir, const std::string & fileName) const
	{
		std::ofstream karFile;
		TryOpenFile(outDir + "/" + fileName, karFile);
		for (size_t i = 0; i < chrList_.size(); ++i)
		{
			int colorId = (i + 1) % CIRCOS_MAX_COLOR;
			karFile << "chr - seq" << i + 1 << " " << chrList_[i].GetDescription() << " 0 " << chrList_[i].GetSequence().length();
			karFile	<< " chr" << colorId << std::endl;
		}
	}

	void OutputGenerator::GenerateD3Output(const BlockList & blockList, const std::string & outFile) const
	{
		std::istringstream htmlTemplate(d3Template);

		//open output file
		std::ofstream out;
		TryOpenFile(outFile, out);

		std::string buffer;
		for(;;)
		{
			std::getline(htmlTemplate, buffer);
			if (buffer != "//SIBELIA_MARK_INSERT")
			{
				out << buffer << std::endl;
			}
			else
			{
				break;
			}
		}

		out << "chart_data = [" << std::endl;

		//blocks must be sorted by start
		BlockList sortedBlocks = blockList;
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
			out << "\"name\":\"" << OutputD3BlockID(*itBlock) << "\",";
			out << "\"size\":" << itBlock->GetLength() << ",";
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
					out << "\"" << OutputD3BlockID(*itPair) << "\"";
				}
			}
			out << "]";
			out << "}" << std::endl;
		}
		out << "];" << std::endl;

		// making data for chart legend
		out << "chart_legend = [" << std::endl;
		first_line = true;
		for(size_t i = 0; i < chrList_.size(); i++)
		{
			if (!first_line)
				out << ",";
			else
				first_line = false;
			out << "    \"seq " << chrList_[i].GetId() + 1 << " : " <<  chrList_[i].GetDescription() << "\"" << std::endl;
		}
		out << "];" << std::endl;

		//write rest of html template
		while (!htmlTemplate.eof())
		{
			std::getline(htmlTemplate, buffer);
			out << buffer << std::endl;
		}
	}

	void OutputGenerator::GenerateVariantOutput(const std::vector <Variant> & variants, const std::string & outFile) const 
	{
		std::ofstream out;
		TryOpenFile(outFile, out);
		
		std::istringstream vcfTemplate(variantTemplate);
		std::string buffer;
		std::getline(vcfTemplate, buffer);
		while (!buffer.empty())
		{
			if (buffer == "##reference=") 
			{
				out << buffer << chrList_[0].GetDescription() << std::endl;
			}
			else if (buffer == "##assembly=") 
			{
				out << buffer << chrList_[1].GetDescription() << std::endl;
			}
			else 
			{
				out << buffer << std::endl;
			}
			std::getline(vcfTemplate, buffer);
		}

		std::string chrName = ".";
		size_t chrNameStart = chrList_[1].GetDescription().find_last_of("|", chrList_[1].GetDescription().size() - 2);
		if (chrNameStart == std::string::npos) 
		{
			chrName = chrList_[1].GetDescription();
		}
		else 
		{
			chrName = chrList_[1].GetDescription().substr(chrNameStart + 1);
			size_t chrNameEnd = chrName.find_last_of(".");
			if (chrNameEnd == std::string::npos) 
			{
				chrNameEnd = chrName.size() - 1; 
			}
			chrName.erase(chrNameEnd, chrName.size() - chrNameEnd);
		}

		for (size_t i = 0; i < variants.size(); ++i)
		{
			out << chrName << "\t"
				<< variants[i].GetReferencePos() << "\t"
				<< ".\t"
				<< (variants[i].GetReferenceAllele().empty() ? "." : variants[i].GetReferenceAllele()) << "\t"
				<< (variants[i].GetAlternativeAllele().empty() ? "." : variants[i].GetAlternativeAllele()) << "\t"
				<< ".\t.\t." << std::endl;
		}
	}

	void OutputGenerator::TryOpenFile(const std::string & fileName, std::ofstream & stream) const
	{
		stream.open(fileName.c_str());
		if(!stream)
		{
			throw std::runtime_error(("Cannot open file " + fileName).c_str());
		}
	}

	void OutputGenerator::TryOpenResourceFile(const std::string & fileName, std::ifstream & stream) const
	{
		std::vector<std::string> dirs = GetResourceDirs();
		for (std::vector<std::string>::iterator itDirs = dirs.begin(); itDirs != dirs.end(); ++itDirs)
		{
			stream.open((*itDirs + "/" + fileName).c_str());
			if (stream) break;
		}
		if (!stream)
		{
			throw std::runtime_error(("Cannot find resource file: " + fileName).c_str());
		}
	}

	void OutputGenerator::OutputBuffer(const std::string & fileName, const std::string & buffer) const
	{
		std::ofstream out;
		TryOpenFile(fileName, out);
		out << buffer;
	}

	void GlueBlock(const std::string block[], size_t blockSize, std::string & buf)
	{
		std::stringstream ss;
		std::copy(block, block + blockSize, std::ostream_iterator<std::string>(ss, "\n"));
	}
}
