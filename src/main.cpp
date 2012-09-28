#include "graphalgorithm.h"

std::string IntToStr(size_t value)
{
	std::stringstream ss;
	ss << value;
	return ss.str();
}

std::vector<std::pair<int, int> > ReadStageFile(const std::string & fileName)
{
	int count;
	std::ifstream stageFile(fileName.c_str());
	if(!stageFile)
	{
		return std::vector<std::pair<int, int> >();
	}

	stageFile >> count;
	std::vector<std::pair<int, int> > ret(count);
	for(int i = 0; i < count; i++)
	{
		stageFile >> ret[i].first >> ret[i].second;
	}

	return ret;
}

int Abs(int x)
{
	return x > 0 ? x : -x;
}

bool CompareBlocksById(SyntenyBuilder::GraphAlgorithm::BlockInstance & a, SyntenyBuilder::GraphAlgorithm::BlockInstance & b)
{
	return Abs(a.GetId()) < Abs(b.GetId());
}

int main(int argc, char * argv[])
{	
	if(argc < 3)
	{
		std::cerr << "Program for building synteny blocks from a FASTA file" << std::endl;
		std::cerr << "Usage: SyntenyBuilder <input filename> <stage list filename> [-dot]" << std::endl;
	}
	else
	{
		std::string fileName(argv[1]);
		std::string stageFile(argv[2]);
		std::cout.sync_with_stdio(false);
		std::vector<std::pair<int, int> > stage = ReadStageFile(stageFile);
		if(stage.empty())
		{
			std::cerr << "Stage file is empty or cannot be open" << std::endl;
			return 1;
		}

		for(size_t i = 0; i < stage.size(); i++)
		{
			if(stage[i].first < 2 || stage[i].second < 0)
			{
				std::cerr << "Incorrect stage record " << std::endl;
				return 1;
			}
		}

		SyntenyBuilder::FASTAReader reader(fileName.c_str());
		if(!reader.IsOk())
		{
			std::cerr << "Can't open the input file" << std::endl;
			return -1;
		}

		bool dot = argc >= 4 && argv[3] == std::string("-dot");
		
		std::vector<SyntenyBuilder::FASTAReader::FASTARecord> record;
		reader.GetSequences(record);
		SyntenyBuilder::DNASequence dnaseq(record);
		SyntenyBuilder::BifurcationStorage bifStorage;
		for(size_t i = 0; i < stage.size(); i++)
		{
			std::cerr << "Building the graph, stage = " << i + 1 << std::endl;			
			if(stage[i].first < dnaseq.TotalSize())
			{
				if(dot)
				{				
					std::ofstream before((fileName + "_stage_" + IntToStr(i + 1) + "_before.dot").c_str());
					SyntenyBuilder::GraphAlgorithm::SerializeGraph(dnaseq, stage[i].first, before);
				}

				std::cerr << "Simplifying the graph, stage = " << i + 1 << std::endl;
				SyntenyBuilder::GraphAlgorithm::EnumerateBifurcations(dnaseq, bifStorage, stage[i].first);
				SyntenyBuilder::GraphAlgorithm::SimplifyGraph(dnaseq, bifStorage, stage[i].first, stage[i].second);

				if(dot)
				{
					std::ofstream after((fileName + "_stage_" + IntToStr(i + 1) + "_after.dot").c_str());
					SyntenyBuilder::GraphAlgorithm::SerializeGraph(dnaseq, stage[i].first, after);
				}
			}
			else
			{
				std::cerr << "ERROR: The sequence is too short!" << std::endl;
				exit(1);
			}
		}
			
		std::string header = fileName;
		std::ofstream chr((header + "_chr").c_str());
		std::ofstream report((header + "_report").c_str());
		std::ofstream blocks((header + "_blocks").c_str());
		std::ofstream indices((header + "_indices").c_str());

		indices << "Chr_id\tChr description" << std::endl;
		blocks << "Chr_id\tChr description" << std::endl;
		for(size_t i = 0; i < record.size(); i++)
		{
			indices << i + 1 << '\t' << record[i].description << std::endl;
			blocks << i + 1 << '\t' << record[i].description << std::endl;
		}

		indices << SyntenyBuilder::DELIMITER << std::endl;
		blocks << SyntenyBuilder::DELIMITER << std::endl; 
		std::cerr << SyntenyBuilder::DELIMITER << std::endl;
		if(stage.back().first < dnaseq.TotalSize())
		{
			size_t k = stage.back().first;
			std::ofstream condensed((fileName + "_condensed.dot").c_str());
			std::cerr << "Finding synteny blocks" << std::endl;
			std::vector<std::vector<SyntenyBuilder::GraphAlgorithm::BlockInstance> > chrList;
			SyntenyBuilder::GraphAlgorithm::EnumerateBifurcations(dnaseq, bifStorage, k);
			SyntenyBuilder::GraphAlgorithm::SerializeCondensedGraph(dnaseq, bifStorage, stage.back().first, condensed);
			SyntenyBuilder::GraphAlgorithm::GenerateSyntenyBlocks(dnaseq, bifStorage, stage.back().first, chrList);
			std::vector<SyntenyBuilder::GraphAlgorithm::BlockInstance> block;
			for(size_t i = 0; i < chrList.size(); i++)
			{
				block.insert(block.end(), chrList[i].begin(), chrList[i].end());
				chr << '>' << record[i].description << std::endl;
				chr.setf(std::ios_base::showpos);
				for(size_t j = 0; j < chrList[i].size(); j++)
				{
					chr << chrList[i][j].GetId() << ' ';
				}

				chr << '$' << std::endl;
			}

			std::sort(block.begin(), block.end(), CompareBlocksById);
			for(size_t now = 0; now < block.size(); )
			{
				size_t prev = now;
				indices << "Block #" << Abs(block[now].GetId()) << std::endl;
				indices << "Chr_id\tStrand\tStart\tEnd\tLength" << std::endl;
				blocks << "Block #" << Abs(block[now].GetId()) << std::endl;
				for(; now < block.size() && Abs(block[now].GetId()) == Abs(block[prev].GetId()); now++)
				{
					std::string & str = record[block[now].GetChr()].sequence;			
					size_t length = block[now].GetEnd() - block[now].GetStart();
					indices << block[now].GetChr() + 1 << '\t' << (block[now].GetId() < 0 ? '-' : '+') << '\t' << block[now].GetStart() << '\t' << block[now].GetEnd() << '\t' << length << std::endl;
					blocks << (block[now].GetId() < 0 ? '-' : '+') << block[now].GetChr() + 1 << ':' << block[now].GetStart() << ':' << block[now].GetEnd() << std::endl;
					if(block[now].GetId() > 0)
					{
						std::copy(str.begin() + block[now].GetStart(), str.begin() + block[now].GetEnd(), std::ostream_iterator<char>(blocks));
					}
					else
					{
						std::string buf(str.begin() + block[now].GetStart(), str.begin() + block[now].GetEnd());
						for(size_t k = 0; k < buf.size(); k++)
						{
							buf[k] = SyntenyBuilder::DNASequence::Translate(buf[k]);
						}

						std::copy(buf.begin(), buf.end(), std::ostream_iterator<char>(blocks));
					}

					blocks << std::endl;
				}

				indices << SyntenyBuilder::DELIMITER << std::endl;
				blocks << SyntenyBuilder::DELIMITER << std::endl;
			}
		}

		std::cerr.setf(std::cerr.fixed);
		std::cerr.precision(2);
		std::cerr << "Time elapsed: " << double(clock()) / CLOCKS_PER_SEC << std::endl;
	}

	return 0;
}