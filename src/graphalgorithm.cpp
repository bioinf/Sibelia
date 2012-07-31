#include "graphalgorithm.h"
//#define _DEBUG
//#undef _DEBUG
namespace SyntenyBuilder
{
	size_t GraphAlgorithm::EnumerateBifurcations(DNASequence & sequence, size_t k, BifurcationStorage & bifStorage)
	{
		bifStorage.Clear();
		std::cerr << DELIMITER << std::endl;
		std::cerr << "Finding all bifurcations in the graph..." << std::endl;		
		KMerIndex index(&sequence);
		index.SetupIndex(k);
		size_t bifurcationCount = 0;
		std::vector<StrandIterator> kmer;
		KMerSet visit(sequence.Size(), KMerIndex::KMerHashFunction(k), KMerIndex::KMerEqualTo(k));
		StrandIterator posBegin = sequence.PositiveBegin();
		StrandIterator negBegin = sequence.NegativeBegin();
		
	#ifdef _DEBUG
		std::cerr << "Found bifurcations:" << std::endl;
	#endif

		for(StrandIterator it = sequence.PositiveBegin(); it.ProperKMer(k); ++it)
		{
			if(visit.find(it) != visit.end())
			{
				continue;
			}

			char forward = -1;
			char backward = -1;
			bool properBifurcation = false;
			index.ListEquivalentKmers(it, kmer);
			for(size_t i = 0; i < kmer.size(); i++)
			{
				char nowForward = -1;
				char nowBackward = -1;
				if(kmer[i] != posBegin && kmer[i] != negBegin)
				{
					nowBackward = *AdvanceBackward(kmer[i], 1);
				}

				if(kmer[i].ProperKMer(k + 1))
				{
					nowForward = *AdvanceForward(kmer[i], k);
				}

				if(nowForward != -1 && forward == -1)
				{
					forward = nowForward;
				}

				if(nowBackward != -1 && backward == -1)
				{
					backward = nowBackward;
				}

				if((nowForward != -1 && nowForward != forward) || (nowBackward != -1 && nowBackward != backward))
				{
					properBifurcation = true;
					break;
				}				
			}

			if(properBifurcation)
			{
				visit.insert(it);				
				boost::function<void (StrandIterator)> adder = boost::bind(
								&BifurcationStorage::AddPoint, 
								boost::ref(bifStorage),
								_1, 
								bifurcationCount++);
				std::for_each(kmer.begin(), kmer.end(), adder);

		#ifdef _DEBUG
				std::cerr << "Bifurcation No. " << bifurcationCount - 1 << std::endl;
				for(size_t i = 0; i < kmer.size(); i++)
				{
					CopyN(kmer[i], k, std::ostream_iterator<char>(std::cerr));
					std::cerr << ", " << (kmer[i].GetDirection() == DNASequence::positive ? '+' : '-') <<
						kmer[i].GetPosition() << std::endl;
				}

				std::cerr << DELIMITER << std::endl;
		#endif
			}
		}

		return bifurcationCount;
	}

	void GraphAlgorithm::SimplifyGraph(DNASequence & sequence, size_t k, size_t minBranchSize)
	{
		size_t totalBulges;
		size_t totalWhirls;
		bool anyChanges = true;
		BifurcationStorage bifStorage;
		do
		{
			totalBulges = 0;
			size_t counter = 0;
			size_t bifurcationCount = GraphAlgorithm::EnumerateBifurcations(sequence, k, bifStorage);
			std::cerr << "Total bifurcations: " << bifurcationCount << std::endl;
			std::cerr << "Removing whirls..." << std::endl;
			totalWhirls = RemoveWhirls(bifStorage, sequence, k, minBranchSize);

			std::cerr << "Removing bulges..." << std::endl;
			for(size_t id = 0; id < bifurcationCount; id++)
			{
				totalBulges += RemoveBulges(bifStorage, sequence, k, minBranchSize, id);
			}

			std::cerr << "Total whirls: " << totalWhirls << std::endl;
			std::cerr << "Total bulges: " << totalBulges << std::endl;		
			sequence.Optimize();
		}
		while(totalBulges > 0 || totalWhirls > 0);
	}
}
