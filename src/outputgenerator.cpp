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

		template<class Iterator, class F, class ReturnType>
			struct FancyIterator
			{
			public:
				FancyIterator& operator++()
				{
					++it;
					return *this;
				}

				FancyIterator operator++(int)
				{
					FancyIterator ret(*this);
					++(*this);
					return ret;
				}

				bool operator == (FancyIterator toCompare) const
				{
					return it == toCompare.it;
				}

				ReturnType operator * () 
				{
					return f(*it);
				}

				FancyIterator() {}
				FancyIterator(Iterator it, F f): it(it), f(f) {}

			private:
				F f;
				Iterator it;
			};

		template<class Iterator, class F, class ReturnType>
			FancyIterator<Iterator, F, ReturnType> CFancyIterator(Iterator it, F f, ReturnType)
			{
				return FancyIterator<Iterator, F, ReturnType>(it, f);
			}

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

	void OutputGenerator::GenerateReport(std::ostream & out) const
	{
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
				out << "All\t" << group.size() << "\t";
			}

			out.precision(2);
			out.setf(std::ostream::fixed);
			std::vector<double> coverage = CalculateCoverage(chrList_, sepBlock.begin() + it->first, sepBlock.begin() + it->second);
			std::copy(coverage.begin(), coverage.end(), std::ostream_iterator<double>(out, "%\t"));
			out << std::endl;
		}

		out << DELIMITER << std::endl;
	}

	void OutputGenerator::ListChromosomesAsPermutations(std::ostream & out) const
	{
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

	void OutputGenerator::ListBlocksIndices(std::ostream & out) const
	{
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

	void OutputGenerator::ListBlocksSequences(std::ostream & out) const
	{
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
}