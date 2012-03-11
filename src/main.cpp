#include "fastareader.h"
#include "debruijngraph.h"

int main(int argc, char * argv[])
{
	if(argc != 3)
	{
		std::cerr << "Program form building synteny blocks from a FASTA file" << std::endl;
		std::cerr << "Usage: SyntenyBuilder <K-mer size> <input filename>" << std::endl;
	}
	else
	{
		int k = atoi(argv[1]);
		std::string fileName(argv[2]);

		if(k > 1)
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
			//graph.DebugOutput(std::cout);
		}
		else
		{
			std::cerr << "K-mer size is incorrect"; 
			return -1;
		}
	}

	return 0;
}