#include "auxiliary.h"

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

		typedef std::vector<std::pair<size_t, std::vector<GraphAlgorithm::BlockInstance> > > GroupedBlockList;
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
			
			ret.push_back(totalCoveredBp / totalBp * 100);
			return ret;
		}		
	}

	void GenerateReport(const ChrList & chrList, std::vector<GraphAlgorithm::BlockInstance> & block, std::ostream & out)
	{
		std::sort(block.begin(), block.end(), CompareBlocksById);
		GroupedBlockList sepBlock;
		for(size_t now = 0; now < block.size(); )
		{
			size_t prev = now;
			for(; now < block.size() && block[now].GetBlockId() == block[prev].GetBlockId(); now++);
			sepBlock.push_back(std::make_pair(now - prev, std::vector<GraphAlgorithm::BlockInstance>(block.begin() + prev, block.begin() + now)));
		}

		out << "Degree\tCount\t";
		for(size_t chr = 0; chr < chrList.size(); chr++)
		{
			out << "Chr " << chr + 1 << '\t';
		}
		
		out << "Total" << std::endl;
		std::sort(sepBlock.begin(), sepBlock.end());
		for(size_t now = 0; now < sepBlock.size(); )
		{
			size_t prev = now;
			for(; now < sepBlock.size() && sepBlock[now].first == sepBlock[prev].first; now++);
			std::vector<double> coverage = CalculateCoverage(chrList, sepBlock.begin() + prev, sepBlock.begin() + now);
			out << sepBlock[prev].first << '\t' << now - prev << '\t';
			out.setf(std::ostream::fixed);
			out.precision(2);
			std::copy(coverage.begin(), coverage.end(), std::ostream_iterator<double>(out, "%\t"));
			out << std::endl;
		}

		std::vector<double> coverage = CalculateCoverage(chrList, sepBlock.begin(), sepBlock.end());
		out << "All" << '\t' << sepBlock.size() << '\t';
		out.setf(std::ostream::fixed);
		out.precision(2);
		std::copy(coverage.begin(), coverage.end(), std::ostream_iterator<double>(out, "%\t"));
		out << std::endl << DELIMITER << std::endl;;
	}
}