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

void PutProgressChr(size_t progress, SyntenyFinder::BlockFinder::State state)
{
	static size_t prev = 0;
	switch(state)
	{
	case SyntenyFinder::BlockFinder::start:
		prev = 0;
		std::cout << '[';
		break;
	case SyntenyFinder::BlockFinder::run:
		if(progress != prev)
		{
			prev = progress;
			std::cout << '.';
			std::cout.flush();
		}

		break;
	case SyntenyFinder::BlockFinder::end:
		std::cout << "]" << std::endl;
		break;
	}
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

		std::vector<SyntenyFinder::FASTARecord> chrList;
		reader.GetSequences(chrList);
		SyntenyFinder::BlockFinder finder(chrList);
		for(size_t i = 0; i < stage.size(); i++)
		{
			std::cout << "Simplification stage " << i + 1 << " of " << stage.size() << std::endl;
			std::cout << "Enumerating vertices of the graph, then performing bulge removal..." << std::endl;
			finder.PerformGraphSimplifications(stage[i].first, stage[i].second, 4, PutProgressChr);
		}
		
		std::string header = fileName;
		std::ofstream chr((header + "_chr").c_str());
		std::ofstream report((header + "_report").c_str());
		std::ofstream blocks((header + "_blocks").c_str());
		std::ofstream indices((header + "_indices").c_str());
		std::vector<SyntenyFinder::BlockInstance> blockList;
		std::cout << "Finding synteny blocks and generating the output..." << std::endl;
		finder.GenerateSyntenyBlocks(stage.back().first, blockList, PutProgressChr);
		SyntenyFinder::OutputGenerator generator(chrList, blockList);
		generator.GenerateReport(report);
		generator.ListBlocksIndices(indices);
		generator.ListBlocksSequences(blocks);
		generator.ListChromosomesAsPermutations(chr);

		std::cout.setf(std::cout.fixed);
		std::cout.precision(2);
		std::cout << "Time elapsed: " << double(clock()) / CLOCKS_PER_SEC << " seconds" << std::endl;
	}

	return 0;
}