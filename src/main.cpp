#include "fastareader.h"
#include "graphalgorithm.h"
/*
void Test()
{
	int k = 10;
	srand(time(0));
	SyntenyBuilder::DNASequence seq(std::string(1000, 'g'));
	seq.KeepHash(k);
	for(SyntenyBuilder::DNASequence::StrandConstIterator it = seq.PositiveBegin(); it != seq.PositiveRightEnd(); ++it)
	{
		SyntenyBuilder::DNASequence::StrandConstIterator jt = it;
		std::advance(jt, k);
		if(jt.Valid())
		{
			it.GetHashCode(k);
		}
	}

	for(SyntenyBuilder::DNASequence::StrandConstIterator it = seq.NegativeBegin(); it != seq.NegativeRightEnd(); ++it)
	{
		SyntenyBuilder::DNASequence::StrandConstIterator jt = it;
		std::advance(jt, k);
		if(jt.Valid())
		{
			it.GetHashCode(k);
		}
	}

	std::string alphabet = "acgt";
	for(SyntenyBuilder::DNASequence::StrandIterator it = seq.PositiveBegin(); it != seq.PositiveRightEnd(); ++it)
	{
		if(rand() % 6 == 0)
		{
			it.AssignBase(alphabet[rand() % 4]);
		}
		else if(rand() % 7 == 0)
		{
			it.Invalidate();
		}
	}

	for(SyntenyBuilder::DNASequence::StrandConstIterator it = seq.PositiveBegin(); it != seq.PositiveRightEnd(); ++it)
	{
		SyntenyBuilder::DNASequence::StrandConstIterator jt = it;
		std::advance(jt, k);
		if(jt.Valid())
		{
			it.GetHashCode(k);
		}
	}

	for(SyntenyBuilder::DNASequence::StrandConstIterator it = seq.NegativeBegin(); it != seq.NegativeRightEnd(); ++it)
	{
		SyntenyBuilder::DNASequence::StrandConstIterator jt = it;
		std::advance(jt, k);
		if(jt.Valid())
		{
			it.GetHashCode(k);
		}
	}
}
*/
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
			std::ofstream indexes((fileName + "_indices").c_str());
			std::ofstream consensus((fileName + "_consensus").c_str());
			std::copy(g.sequence.PositiveBegin(), g.sequence.PositiveRightEnd(), 
				std::ostream_iterator<char>(consensus));

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