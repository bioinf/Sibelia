#include "auxiliary.h"

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

std::vector<std::pair<int, int> > DefaultStageFile()
{
	std::pair<int, int> stage[] = 
	{
		std::make_pair(30, 150),
		std::make_pair(100, 1000),
		std::make_pair(1000, 5000),		
		std::make_pair(5000, 15000)
	};

	return std::vector<std::pair<int, int> >(stage, stage + sizeof(stage) / sizeof(stage[0]));
}

void PutProgressChr(size_t)
{
	std::cout << '.';
	std::cout.flush();
}

const std::string DELIMITER(80, '-');

int main(int argc, char * argv[])
{	
	if(argc != 2)
	{
		std::cerr << "Program for building synteny blocks from a FASTA file" << std::endl;
		std::cerr << "Usage: SyntenyFinder <input filename>" << std::endl; // <stage list filename> [-dot]" << std::endl;
	}
	else
	{
		std::string fileName(argv[1]);
	//	std::string stageFile(argv[2]);
		std::cout.sync_with_stdio(false);
		//std::vector<std::pair<int, int> > stage = ReadStageFile(stageFile);
		//if(stage.empty())
	//	{
//			std::cerr << "Stage file is empty or cannot be open" << std::endl;
			//return 1;
		//}

		std::vector<std::pair<int, int> > stage = DefaultStageFile();
		for(size_t i = 0; i < stage.size(); i++)
		{
			if(stage[i].first < 2 || stage[i].second < 0)
			{
				std::cerr << "Incorrect stage record " << std::endl;
				return 1;
			}
		}

		SyntenyFinder::FASTAReader reader(fileName.c_str());
		if(!reader.IsOk())
		{
			std::cerr << "Can't open the input file" << std::endl;
			return -1;
		}

		bool dot = false; //argc >= 4 && argv[3] == std::string("-dot");
		
		std::vector<SyntenyFinder::FASTAReader::FASTARecord> record;
		reader.GetSequences(record);
		SyntenyFinder::DNASequence dnaseq(record);
		SyntenyFinder::BifurcationStorage bifStorage;
		for(size_t i = 0; i < stage.size(); i++)
		{
			if(stage[i].first < dnaseq.TotalSize())
			{
				if(dot)
				{				
					std::ofstream before((fileName + "_stage_" + IntToStr(i + 1) + "_before.dot").c_str());
					SyntenyFinder::GraphAlgorithm::SerializeGraph(dnaseq, stage[i].first, before);
				}

				std::cout << "Simplification stage " << i + 1 << " of " << stage.size() << std::endl;
				std::cout << "Enumerating vertices of the graph..." << std::endl << "[";			
				SyntenyFinder::GraphAlgorithm::EnumerateBifurcations(dnaseq, bifStorage, stage[i].first, boost::bind(PutProgressChr, _1));
				std::cout << "]" << std::endl << "Performing bulge removal..." << std::endl << "[";
				size_t bulges = SyntenyFinder::GraphAlgorithm::SimplifyGraph(dnaseq, bifStorage, stage[i].first, stage[i].second, 4, boost::bind(PutProgressChr, _1));
				std::cout << "]" << std::endl << std::endl;

				if(dot)
				{
					std::ofstream after((fileName + "_stage_" + IntToStr(i + 1) + "_after.dot").c_str());
					SyntenyFinder::GraphAlgorithm::SerializeGraph(dnaseq, stage[i].first, after);
				}
			}
			else
			{
				std::cerr << "ERROR: The sequence is too short!" << std::endl;
				std::cerr << "The sequence must be at least " << stage[i].first << " bp long " << std::endl;
				exit(1);
			}
		}
			
		std::string header = fileName;
		std::ofstream chr((header + "_chr").c_str());
		std::ofstream report((header + "_report").c_str());
		std::ofstream blocks((header + "_blocks").c_str());
		std::ofstream indices((header + "_indices").c_str());

		std::ostream * outFiles[] = {&report, &blocks, &indices};
		for(size_t i = 0; i < 3; i++)
		{
			*outFiles[i] << "Chr_id\tSize\tDescription" << std::endl;
		}

		for(size_t i = 0; i < record.size(); i++)
		{
			for(size_t j = 0; j < 3; j++)
			{
				*outFiles[j] << i + 1 << '\t' << record[i].sequence.size() << '\t' << record[i].description << std::endl;
			}
		}

		indices << DELIMITER << std::endl;
		blocks << DELIMITER << std::endl; 
		report << DELIMITER << std::endl; 
		if(stage.back().first < dnaseq.TotalSize())
		{
			size_t k = stage.back().first;
	//		std::ofstream condensed((fileName + "_condensed.dot").c_str());
			std::cout << "Finding synteny blocks and generating the output..." << std::endl;
			std::vector<std::vector<SyntenyFinder::GraphAlgorithm::BlockInstance> > chrList;
			std::cout << "Enumerating vertices of the graph..." << std::endl << "[";			
			SyntenyFinder::GraphAlgorithm::EnumerateBifurcations(dnaseq, bifStorage, stage.back().first, boost::bind(PutProgressChr, _1));
			std::cout << "]" << std::endl;
	//		SyntenyFinder::GraphAlgorithm::SerializeCondensedGraph(dnaseq, bifStorage, stage.back().first, condensed);
			SyntenyFinder::GraphAlgorithm::GenerateSyntenyBlocks(dnaseq, bifStorage, stage.back().first, chrList);
			std::vector<SyntenyFinder::GraphAlgorithm::BlockInstance> block;
			for(size_t i = 0; i < chrList.size(); i++)
			{
				block.insert(block.end(), chrList[i].begin(), chrList[i].end());
				chr << '>' << record[i].description << std::endl;
				chr.setf(std::ios_base::showpos);
				for(size_t j = 0; j < chrList[i].size(); j++)
				{
					chr << chrList[i][j].GetSignedBlockId() << ' ';
				}

				chr << '$' << std::endl;
			}

			std::sort(block.begin(), block.end(), SyntenyFinder::CompareBlocksById);
			for(size_t now = 0; now < block.size(); )
			{
				size_t prev = now;
				indices << "Block #" << block[now].GetBlockId() << std::endl;
				indices << "Chr_id\tStrand\tStart\tEnd\tLength" << std::endl;
				blocks << "Block #" << block[now].GetBlockId() << std::endl;
				for(; now < block.size() && block[now].GetBlockId() == block[prev].GetBlockId(); now++)
				{
					std::string & str = record[block[now].GetChr()].sequence;			
					size_t length = block[now].GetEnd() - block[now].GetStart();
					indices << block[now].GetChr() + 1 << '\t' << (block[now].GetSignedBlockId() < 0 ? '-' : '+') << '\t' << block[now].GetStart() << '\t' << block[now].GetEnd() << '\t' << length << std::endl;
					blocks << (block[now].GetSignedBlockId() < 0 ? '-' : '+') << block[now].GetChr() + 1 << ':' << block[now].GetStart() << ':' << block[now].GetEnd() << std::endl;
					if(block[now].GetSignedBlockId() > 0)
					{
						std::copy(str.begin() + block[now].GetStart(), str.begin() + block[now].GetEnd(), std::ostream_iterator<char>(blocks));
					}
					else
					{
						std::string buf(str.begin() + block[now].GetStart(), str.begin() + block[now].GetEnd());
						for(size_t k = 0; k < buf.size(); k++)
						{
							buf[k] = SyntenyFinder::DNASequence::Translate(buf[k]);
						}

						std::copy(buf.rbegin(), buf.rend(), std::ostream_iterator<char>(blocks));
					}

					blocks << std::endl;
				}

				indices << DELIMITER << std::endl;
				blocks << DELIMITER << std::endl;				
			}

			SyntenyFinder::GenerateReport(record, block, report);
		}

		std::cout.setf(std::cout.fixed);
		std::cout.precision(2);
		std::cout << "Time elapsed: " << double(clock()) / CLOCKS_PER_SEC << " seconds" << std::endl;
	}

	return 0;
}