#include "graphalgorithm.h"
//#define _DEBUG
//#undef _DEBUG
namespace SyntenyBuilder
{
	namespace
	{
		void WhirlRemovalBifurcationUpdate(BifurcationStorage & bifStorage, 
			DNASequence & sequence, size_t k, StrandIterator it, size_t step)
		{
			StrandIterator jt = AdvanceForward(it, step);
			std::vector<std::pair<size_t, size_t> > save;
			for(size_t i = 1; i < k; i++)
			{
				size_t bifurcation = bifStorage.GetBifurcation(++jt);
				if(bifurcation != BifurcationStorage::NO_BIFURCATION)
				{
					save.push_back(std::make_pair(i, bifurcation));
				}
			}

			jt = AdvanceForward(it, k);
			for(size_t i = 0; i < step; i++, ++jt)
			{
				bifStorage.ErasePoint(jt);
				bifStorage.ErasePoint(jt.Invert());
			}

			jt = it;
			for(size_t i = 1; i < k; i++)
			{
				bifStorage.ErasePoint(++jt);
			}

			size_t near = 0;
			for(size_t i = 1; i < k && near < save.size(); i++)
			{
				++it;
				if(i == save[near].first)
				{
					bifStorage.AddPoint(it, save[near++].second);
				}
			}
		}
	}

	size_t GraphAlgorithm::RemoveWhirls(BifurcationStorage & bifStorage, DNASequence & sequence, 
		size_t k, size_t minBranchSize)
	{
		size_t ret = 0;
		for(StrandIterator it = sequence.PositiveBegin(); it != sequence.PositiveRightEnd(); ++it)
		{
			size_t bifId = bifStorage.GetBifurcation(it);
			std::vector<StrandIterator> kmer(1, it);
			if(bifId != BifurcationStorage::NO_BIFURCATION)
			{
				bool remove = true;
				while(remove)
				{
					remove = false;
					StrandIterator jt = it;
					for(size_t step = 1; step < minBranchSize; step++)
					{
						if(bifStorage.GetBifurcation(++jt) == bifId)
						{
						#ifdef _DEBUG
							std::cerr << "Whirl #" << ret++ << std::endl;
							std::cerr << "Before: " << std::endl;
							PrintRaw(sequence, std::cerr);
							std::cerr << "Source branch: " << std::endl;			
							PrintPath(it, k, step, std::cerr);
							bifStorage.Dump(std::cerr);
						#endif

							WhirlRemovalBifurcationUpdate(bifStorage, sequence, k, it, step);
							sequence.EraseN(AdvanceForward(it, k), step);
							remove = true;

						#ifdef _DEBUG
							std::cerr << "After: " << std::endl;
							PrintRaw(sequence, std::cerr);
							std::cerr << "Source branch: " << std::endl;			
							PrintPath(it, k, 0, std::cerr);
							bifStorage.Dump(std::cerr);
							std::cerr << DELIMITER << std::endl;
							GraphAlgorithm::Test(sequence, bifStorage, k);
						#endif

							break;
						}
					}
				}
			}
		}
			
		return ret;
	}
}
