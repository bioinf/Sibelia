#include "graphalgorithm.h"
//#define _DEBUG
//#undef _DEBUG
namespace SyntenyBuilder
{
	namespace
	{
		void ClearVisit(VertexVisitMap & visit,
			BifurcationStorage & bifStorage,
			const std::vector<StrandIterator> & startKMer,
			VisitData targetData)
		{	
			StrandIterator target = startKMer[targetData.kmerId];
			for(size_t i = 0; i < targetData.distance; i++)
			{
				size_t bifurcation = bifStorage.GetBifurcation(++target);
				if(bifurcation != BifurcationStorage::NO_BIFURCATION)
				{
					std::pair<VertexVisitMap::iterator, VertexVisitMap::iterator>
						range = visit.equal_range(bifurcation);
					for(VertexVisitMap::iterator it = range.first; it != range.second; )
					{
						if(targetData.kmerId == it->second.kmerId)
						{
							it = visit.erase(it);						
						}
						else
						{
							++it;
						}
					}
				}
			}
		}

		IteratorPair InvertRange(IteratorPair it, size_t k)
		{
			it.first.Jump(k - 2);
			it.second.Jump(k - 2);
			it.first.Invert();
			it.second.Invert();
			return it;
		}

		void EraseBifurcations(DNASequence & sequence,
			BifurcationStorage & bifStorage,
			size_t k,
			const std::vector<StrandIterator> & startKMer,
			VisitData targetData)
		{
			IteratorPair target = std::make_pair(AdvanceForward(startKMer[targetData.kmerId], 1),
				AdvanceForward(startKMer[targetData.kmerId], targetData.distance));
			IteratorPair inverseTarget = InvertRange(target, k);

			while(target.first != target.second)
			{
				bifStorage.ErasePoint(++target.first);
			}

			while(inverseTarget.first != inverseTarget.second)
			{
				bifStorage.ErasePoint(++inverseTarget.first);
			}
		}

		void UpdateBifurcations(DNASequence & sequence,
			BifurcationStorage & bifStorage,
			size_t k,
			const std::vector<StrandIterator> & startKMer,
			VisitData sourceData,
			VisitData targetData)
		{
			IteratorPair target = std::make_pair(AdvanceForward(startKMer[targetData.kmerId], 1),
				AdvanceForward(startKMer[targetData.kmerId], sourceData.distance));
			IteratorPair inverseTarget = InvertRange(target, k);
			IteratorPair source = std::make_pair(AdvanceForward(startKMer[sourceData.kmerId], 1),
				AdvanceForward(startKMer[sourceData.kmerId], sourceData.distance));
			IteratorPair inverseSource = InvertRange(source, k);
			for(;target.first != target.second; ++target.first, ++source.first)
			{
				size_t bifurcation = bifStorage.GetBifurcation(source.first);
				if(bifurcation != BifurcationStorage::NO_BIFURCATION)
				{
					bifStorage.AddPoint(source.first, bifurcation);
				}
			}

			for(; inverseTarget.first != inverseTarget.second; ++inverseTarget.first, ++inverseSource.first)
			{
				size_t bifurcation = bifStorage.GetBifurcation(inverseTarget.first);
				if(bifurcation != BifurcationStorage::NO_BIFURCATION)
				{
					bifStorage.AddPoint(inverseSource.first, bifurcation);
				}
			}
		}

		void CollapseBulge(DNASequence & sequence,
			BifurcationStorage & bifStorage,
			VertexVisitMap & visit,
			size_t k,
			const std::vector<StrandIterator> & startKMer,
			VisitData sourceData,
			VisitData targetData)
		{
		#ifdef _DEBUG
			static size_t bulge = 0;
			std::cerr << "Bulge #" << bulge++ << std::endl;
			std::cerr << "Before: " << std::endl;
			GraphAlgorithm::PrintRaw(sequence, std::cerr);
			std::cerr << "Source branch: " << std::endl;			
			GraphAlgorithm::PrintPath(startKMer[sourceData.kmerId], k, sourceData.distance, std::cerr);
			std::cerr << "Target branch: " << std::endl;			
			GraphAlgorithm::PrintPath(startKMer[targetData.kmerId], k, targetData.distance, std::cerr);
		#endif

			ClearVisit(visit, bifStorage, startKMer, targetData);
			EraseBifurcations(sequence, bifStorage, k, startKMer, sourceData);
			StrandIterator sourceIt = startKMer[sourceData.kmerId];
			StrandIterator targetIt = startKMer[targetData.kmerId];
			size_t diff = targetData.distance - sourceData.distance;
			sequence.CopyN(sourceIt, sourceData.distance, targetIt);
			targetIt.Jump(sourceData.distance);
			sequence.EraseN(targetIt, diff);
			UpdateBifurcations(sequence, bifStorage, k, startKMer, sourceData, targetData);

		#ifdef _DEBUG
			std::cerr << "After: " << std::endl;
			GraphAlgorithm::PrintRaw(sequence, std::cerr);
			std::cerr << "Source branch: " << std::endl;			
			GraphAlgorithm::PrintPath(startKMer[sourceData.kmerId], k, sourceData.distance, std::cerr);
			std::cerr << "Target branch: " << std::endl;			
			GraphAlgorithm::PrintPath(startKMer[targetData.kmerId], k, sourceData.distance, std::cerr);
			std::cerr << DELIMITER << std::endl;
		#endif
		}	
	}

	size_t GraphAlgorithm::RemoveBulges(BifurcationStorage & bifStorage, 
		DNASequence & sequence, size_t k, size_t minBranchSize, size_t bifId)
	{	
		size_t ret = 0;
		VertexVisitMap visit;
		std::vector<StrandIterator> nowVertex;
		std::vector<StrandIterator> startKMer;
		if(bifStorage.ListPositions(bifId, std::back_inserter(startKMer), sequence) < 2)
		{
			return ret;
		}

		nowVertex.resize(startKMer.size());
		std::vector<char> endChar(startKMer.size(), ' ');
		for(size_t i = 0; i < startKMer.size(); i++)
		{
			if(startKMer[i].ProperKMer(k + 1))
			{                    
				nowVertex[i] = AdvanceForward(startKMer[i], 1);
				endChar[i] = *AdvanceForward(startKMer[i], k);
			}
		}

		std::vector<size_t> travelRange(startKMer.size(), 1);
		for(size_t step = 0; step < minBranchSize; step++)
		{
			for(size_t kmerId = 0; kmerId < nowVertex.size(); kmerId++)
			{
				StrandIterator & kmer = nowVertex[kmerId];
				if(kmer.Valid())
				{
					bool collapsed = false;
					size_t bifurcation = bifStorage.GetBifurcation(kmer);
					VisitData nowData = VisitData(kmerId, travelRange[kmerId]);
					if(bifurcation != BifurcationStorage::NO_BIFURCATION)
					{
						if(bifurcation == bifId)
						{
							continue;
						}

						std::pair<VertexVisitMap::iterator, VertexVisitMap::iterator> range = visit.equal_range(bifurcation);						
						for(VertexVisitMap::iterator it = range.first; it != range.second; ++it)
						{			
							StrandIterator opposite = AdvanceForward(kmer, k - 1).Invert();
							if(endChar[kmerId] != endChar[it->second.kmerId] && opposite != startKMer[it->second.kmerId])
							{
								ret++;
								collapsed = true;
								CollapseBulge(sequence, bifStorage, visit, k, startKMer, it->second, nowData);
								travelRange[kmerId] = it->second.distance;
								endChar[kmerId] = endChar[it->second.kmerId];
								break;
							}
						}

						if(!collapsed)
						{
							visit.insert(std::make_pair(bifurcation, nowData));
						}
					}
						
					if(!collapsed)
					{
						++kmer;
						travelRange[kmerId]++;
					}
				}
			}
		}

		return ret;
	}
}