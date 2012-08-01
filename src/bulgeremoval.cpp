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
			return std::make_pair(it.second.Invert(), it.first.Invert());
		}

		void EraseBifurcations(DNASequence & sequence,
			BifurcationStorage & bifStorage,
			size_t k,
			const std::vector<StrandIterator> & startKMer,
			VisitData targetData,
			std::vector<std::pair<size_t, size_t> > & lookBack)
		{
			lookBack.clear();
			IteratorPair target = std::make_pair(startKMer[targetData.kmerId],
				AdvanceForward(startKMer[targetData.kmerId], targetData.distance + 1));
			IteratorPair inverseTarget = InvertRange(target, k);
			IteratorPair inverseKMer = std::make_pair(AdvanceForward(startKMer[targetData.kmerId], k - 1).Invert(),
				AdvanceBackward(startKMer[targetData.kmerId], 1).Invert());

			for(size_t step = 0; inverseKMer.first != inverseKMer.second; ++inverseKMer.first, step++)
			{
				size_t bifurcation = bifStorage.GetBifurcation(inverseKMer.first);
				if(bifurcation != BifurcationStorage::NO_BIFURCATION)
				{
					lookBack.push_back(std::make_pair(step, bifurcation));
				}
			}

			for(;target.first != target.second; ++target.first)
			{
				bifStorage.ErasePoint(target.first);
			}

			for(;inverseTarget.first != inverseTarget.second; ++inverseTarget.first)
			{
				bifStorage.ErasePoint(inverseTarget.first);
			}
		}

		void UpdateBifurcations(DNASequence & sequence,
			BifurcationStorage & bifStorage,
			size_t k,
			const std::vector<StrandIterator> & startKMer,
			VisitData sourceData,
			VisitData targetData,
			std::vector<std::pair<size_t, size_t> > & lookBack)
		{
			IteratorPair target = std::make_pair(startKMer[targetData.kmerId],
				AdvanceForward(startKMer[targetData.kmerId], sourceData.distance + 1));
			IteratorPair inverseTarget = InvertRange(target, k);
			IteratorPair source = std::make_pair(startKMer[sourceData.kmerId],
				AdvanceForward(startKMer[sourceData.kmerId], sourceData.distance + 1));
			IteratorPair inverseSource = InvertRange(source, k);
			IteratorPair inverseKMer = std::make_pair(AdvanceForward(startKMer[targetData.kmerId], k - 1).Invert(),
				AdvanceBackward(startKMer[targetData.kmerId], 1).Invert());

			size_t near = 0;
			for(size_t i = 0; i < k && near < lookBack.size(); i++, ++inverseKMer.first)
			{
				if(i == lookBack[near].first)
				{
					bifStorage.AddPoint(inverseKMer.first, lookBack[near++].second);
				}
			}

			for(;target.first != target.second; ++target.first, ++source.first)
			{
				size_t bifurcation = bifStorage.GetBifurcation(source.first);
				if(bifurcation != BifurcationStorage::NO_BIFURCATION)
				{
					bifStorage.AddPoint(target.first, bifurcation);
				}
			}

			for(; inverseTarget.first != inverseTarget.second; ++inverseTarget.first, ++inverseSource.first)
			{
				size_t bifurcation = bifStorage.GetBifurcation(inverseSource.first);
				if(bifurcation != BifurcationStorage::NO_BIFURCATION)
				{
					bifStorage.AddPoint(inverseTarget.first, bifurcation);
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
			bifStorage.Dump(std::cerr);
		#endif

			ClearVisit(visit, bifStorage, startKMer, targetData);
			std::vector<std::pair<size_t, size_t> > lookBack;
			EraseBifurcations(sequence, bifStorage, k, startKMer, targetData, lookBack);
			StrandIterator sourceIt = startKMer[sourceData.kmerId];
			StrandIterator targetIt = startKMer[targetData.kmerId];
			size_t diff = targetData.distance - sourceData.distance;
			sequence.CopyN(sourceIt, sourceData.distance, targetIt);
			targetIt.Jump(sourceData.distance);
			sequence.EraseN(targetIt, diff);
			UpdateBifurcations(sequence, bifStorage, k, startKMer, sourceData, targetData, lookBack);

		#ifdef _DEBUG
			std::cerr << "After: " << std::endl;
			GraphAlgorithm::PrintRaw(sequence, std::cerr);
			std::cerr << "Source branch: " << std::endl;			
			GraphAlgorithm::PrintPath(startKMer[sourceData.kmerId], k, sourceData.distance, std::cerr);
			std::cerr << "Target branch: " << std::endl;			
			GraphAlgorithm::PrintPath(startKMer[targetData.kmerId], k, sourceData.distance, std::cerr);
			bifStorage.Dump(std::cerr);
			std::cerr << DELIMITER << std::endl;
			GraphAlgorithm::Test(sequence, bifStorage, k);
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

		std::map<size_t, size_t> restrict;
		nowVertex.resize(startKMer.size());
		std::vector<char> endChar(startKMer.size(), ' ');
		for(size_t i = 0; i < startKMer.size(); i++)
		{
			if(startKMer[i].ProperKMer(k + 1))
			{                    
				nowVertex[i] = AdvanceForward(startKMer[i], 1);
				endChar[i] = *AdvanceForward(startKMer[i], k);
			}

			StrandIterator jt(startKMer[i]);
			for(size_t j = 0; j < k; j++)
			{
				restrict[jt.GetPosition()] = i;
				++jt;
			}
		}

		std::vector<size_t> travelRange(startKMer.size(), 1);
		for(size_t step = 1; step < minBranchSize; step++)
		{
			for(size_t kmerId = 0; kmerId < nowVertex.size(); kmerId++)
			{
				StrandIterator & kmer = nowVertex[kmerId];
				std::map<size_t, size_t>::iterator rt = restrict.find(kmer.GetPosition());
				if(kmer.Valid() && (rt == restrict.end() || rt->second == kmerId))
				{
					bool collapsed = false;
					size_t bifurcation = bifStorage.GetBifurcation(kmer);
					if(bifurcation != BifurcationStorage::NO_BIFURCATION)
					{
						if(bifurcation == bifId)
						{
							kmer = sequence.PositiveRightEnd();
							continue;
						}

						VisitData nowData(kmerId, travelRange[kmerId]);
						std::pair<VertexVisitMap::iterator, VertexVisitMap::iterator> range = visit.equal_range(bifurcation);						
						for(VertexVisitMap::iterator it = range.first; it != range.second; ++it)
						{			
							VisitData prevData = it->second;
							StrandIterator opposite = AdvanceForward(kmer, k - 1).Invert();
							if(endChar[nowData.kmerId] != endChar[prevData.kmerId] && opposite != startKMer[prevData.kmerId]
								&& nowData.distance >= prevData.distance)
							{
								ret++;
								collapsed = true;								
								CollapseBulge(sequence, bifStorage, visit, k, startKMer, prevData, nowData);
								travelRange[nowData.kmerId] = prevData.distance;
								endChar[nowData.kmerId] = endChar[prevData.kmerId];
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