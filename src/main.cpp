#include "outputgenerator.h"

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
	if(argc != 3)
	{
		std::cerr << "Program for building synteny blocks from a FASTA file" << std::endl;
		std::cerr << "Usage: SyntenyFinder <input filename>" << std::endl; // <stage list filename> [-dot]" << std::endl;
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

		SyntenyFinder::FASTAReader reader(fileName.c_str());
		if(!reader.IsOk())
		{
			std::cerr << "Can't open the input file" << std::endl;
			return -1;
		}

		bool dot = false; //argc >= 4 && argv[3] == std::string("-dot");
		
		std::vector<SyntenyFinder::FASTAReader::FASTARecord> chrList;
		reader.GetSequences(chrList);
		SyntenyFinder::DNASequence dnaseq(chrList);
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
		std::vector<SyntenyFinder::BlockInstance> blockList;

		if(stage.back().first < dnaseq.TotalSize())
		{
			size_t k = stage.back().first;
	//		std::ofstream condensed((fileName + "_condensed.dot").c_str());
			std::cout << "Finding synteny blocks and generating the output..." << std::endl;
			std::cout << "Enumerating vertices of the graph..." << std::endl << "[";			
			SyntenyFinder::GraphAlgorithm::EnumerateBifurcations(dnaseq, bifStorage, stage.back().first, boost::bind(PutProgressChr, _1));
			std::cout << "]" << std::endl;
	//		SyntenyFinder::GraphAlgorithm::SerializeCondensedGraph(dnaseq, bifStorage, stage.back().first, condensed);

			SyntenyFinder::GraphAlgorithm::GenerateSyntenyBlocks(dnaseq, bifStorage, stage.back().first, blockList);
			SyntenyFinder::OutputGenerator generator(chrList, blockList);
			generator.GenerateReport(report);
			generator.ListBlocksIndices(indices);
			generator.ListBlocksSequences(blocks);
			generator.ListChromosomesAsPermutations(chr);
		}

		std::cout.setf(std::cout.fixed);
		std::cout.precision(2);
		std::cout << "Time elapsed: " << double(clock()) / CLOCKS_PER_SEC << " seconds" << std::endl;
	}

	return 0;
}