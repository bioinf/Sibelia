#include "fastareader.h"
#include "debruijngraph.h"

int main(int argc, char * argv[])
{
	if(argc < 4)
	{
		std::cerr << "Program form building synteny blocks from a FASTA file" << std::endl;
		std::cerr << "Usage: SyntenyBuilder <k-mer size> <minimum cycle size> <input filename> [-debug]" << std::endl;
	}
	else
	{
		int k = atoi(argv[1]);
		int d = atoi(argv[2]);
		std::string fileName(argv[3]);

		if(k > 1 && k > 1)
		{
			SyntenyBuilder::FASTAReader reader(fileName.c_str());
			if(!reader.IsOk())
			{
				std::cerr << "Can't open the input file" << std::endl;
				return -1;
			}

			std::string sequence;
			std::vector<std::string> path;
			reader.GetSequence(sequence);
			for(size_t i = 0; i < sequence.size(); i++)
			{
				sequence[i] = tolower(sequence[i]);
			}

			sequence.erase(std::remove(sequence.begin(), sequence.end(), 'n'), sequence.end());
			SyntenyBuilder::DeBruijnGraph graph(sequence, k);
			bool debug = argc > 4 && argv[4] == std::string("-debug");
			if(debug)
			{
				std::cerr << "Before simplification: " << std::endl;
				graph.DebugOutput(std::cerr);
				std::cerr << std::string(80, '-');
			}
			
			//graph.Simplify(d);

			if(debug)
			{
				std::cerr << "After simplification: " << std::endl;
				graph.DebugOutput(std::cerr);
				std::cerr << std::string(80, '-');
			}

			graph.ListNonBranchingPaths(std::cout);
		}
		else
		{
			if(k <= 1)
			{
				std::cerr << "K-mer size is incorrect" << std::endl;
			}

			if(d < 0)
			{
				std::cerr << "Cycle size is incorrect" << std::endl;
			}

			return -1;
		}
	}

	return 0;
}