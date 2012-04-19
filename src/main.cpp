#include "fastareader.h"
#include "graphalgorithm.h"

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
			std::cout.sync_with_stdio(false);
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
			std::cerr << "Building the graph...";
			SyntenyBuilder::DeBruijnGraph g(sequence, k);
			std::cerr << "Done." << std::endl;
			bool debug = argc > 4 && argv[4] == std::string("-debug");
			if(debug)
			{
				std::ofstream before((fileName + "_before.dot").c_str());
				before << g;
			}

			std::cerr << "Simplifying the graph...";
			SyntenyBuilder::GraphAlgorithm::Simplify(g, d);
			std::cerr << "Done." << std::endl;

			if(debug)
			{
				std::ofstream after((fileName + "_after.dot").c_str());
				after << g;
			}

			std::ofstream general((fileName + "_blocks").c_str());
			std::ofstream indexes((fileName + "_indexes").c_str());

			std::cerr << "Finding non-branching paths...";
			SyntenyBuilder::GraphAlgorithm::ListNonBranchingPaths(g, general, indexes);
			std::cerr << "Done." << std::endl;
			std::cerr.setf(std::cerr.fixed);
			std::cerr.precision(2);
			std::cerr << "Time elapsed: " << double(clock()) / 1000 << std::endl;
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