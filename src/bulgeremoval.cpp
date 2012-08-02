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
			std::vector<std::pair<size_t, size_t> > & lookForward,
			std::vector<std::pair<size_t, size_t> > & lookBack)
		{
			lookBack.clear();
			lookForward.clear();
			StrandIterator amer = AdvanceForward(startKMer[targetData.kmerId], k - 1).Invert();
			StrandIterator bmer = AdvanceForward(startKMer[targetData.kmerId], targetData.distance);
			for(size_t i = 0; i < k; i++, ++amer, ++bmer)
			{
				size_t bifId = bifStorage.GetBifurcation(amer);
				if(bifId != BifurcationStorage::NO_BIFURCATION)
				{
					bifStorage.ErasePoint(amer);
					lookBack.push_back(std::make_pair(i, bifId));
				}

				bifId = bifStorage.GetBifurcation(bmer);
				if(bifId != BifurcationStorage::NO_BIFURCATION)
				{
					bifStorage.ErasePoint(bmer);
					lookForward.push_back(std::make_pair(i, bifId));
				}
			}

			amer = startKMer[targetData.kmerId];
			bmer = AdvanceForward(amer, k + targetData.distance - 1).Invert();
			for(size_t i = 0; i < k + targetData.distance; i++, ++amer, ++bmer)
			{
				bifStorage.ErasePoint(amer);
				bifStorage.ErasePoint(bmer);
			}
		}

		void UpdateBifurcations(DNASequence & sequence,
			BifurcationStorage & bifStorage,
			size_t k,
			const std::vector<StrandIterator> & startKMer,
			VisitData sourceData,
			VisitData targetData,
			const std::vector<std::pair<size_t, size_t> > & lookForward,
			const std::vector<std::pair<size_t, size_t> > & lookBack)
		{
			size_t anear = 0;
			size_t bnear = 0;
			StrandIterator amer = AdvanceForward(startKMer[targetData.kmerId], k - 1).Invert();
			StrandIterator bmer = AdvanceForward(startKMer[targetData.kmerId], sourceData.distance);
			for(size_t i = 0; i < k; i++, ++amer, ++bmer)
			{
				if(anear < lookBack.size() && i == lookBack[anear].first)
				{
					bifStorage.AddPoint(amer, lookBack[anear++].second);
				}

				if(bnear < lookForward.size() && i == lookForward[bnear].first)
				{
					bifStorage.AddPoint(bmer, lookForward[bnear++].second);
				}
			}

			amer = startKMer[targetData.kmerId];
			bmer = AdvanceForward(startKMer[targetData.kmerId], sourceData.distance + k - 1).Invert();
			StrandIterator srcAMer = startKMer[sourceData.kmerId];
			StrandIterator srcBMer = AdvanceForward(startKMer[sourceData.kmerId], sourceData.distance + k - 1).Invert();
			for(size_t i = 0; i < sourceData.distance + 1; i++, ++amer, ++bmer, ++srcAMer, ++srcBMer)
			{
				size_t bifId = bifStorage.GetBifurcation(srcAMer);
				if(bifId != BifurcationStorage::NO_BIFURCATION)
				{
					bifStorage.AddPoint(amer, bifId);
				}

				bifId = bifStorage.GetBifurcation(srcBMer);
				if(bifId != BifurcationStorage::NO_BIFURCATION)
				{
					bifStorage.AddPoint(bmer, bifId);
				}
			}
		}

		bool Overlap(size_t k,
			const std::vector<StrandIterator> & startKMer,
			VisitData sourceData,
			VisitData targetData)
		{
			std::vector<size_t> occur;
			StrandIterator it = startKMer[sourceData.kmerId];
			for(size_t i = 0; i < sourceData.distance + k; i++, ++it)
			{
				occur.push_back(it.GetPosition());
			}

			it = startKMer[targetData.kmerId];
			std::sort(occur.begin(), occur.end());
			for(size_t i = 0; i < targetData.distance + k; i++, ++it)
			{
				if(std::binary_search(occur.begin(), occur.end(), it.GetPosition()))
				{
					return true;
				}
			}

			return false;
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
			std::vector<std::pair<size_t, size_t> > lookForward;
			std::vector<std::pair<size_t, size_t> > lookBack;
			EraseBifurcations(sequence, bifStorage, k, startKMer, targetData, lookForward, lookBack);
			StrandIterator sourceIt = startKMer[sourceData.kmerId];
			StrandIterator targetIt = startKMer[targetData.kmerId];
			size_t diff = targetData.distance - sourceData.distance;
			sequence.CopyN(sourceIt, sourceData.distance, targetIt);
			targetIt.Jump(sourceData.distance);
			sequence.EraseN(targetIt, diff);
			UpdateBifurcations(sequence, bifStorage, k, startKMer, sourceData, targetData, lookForward, lookBack);

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
		}

		std::vector<size_t> travelRange(startKMer.size(), 1);
		for(size_t step = 1; step < minBranchSize; step++)
		{
			for(size_t kmerId = 0; kmerId < nowVertex.size(); kmerId++)
			{
				StrandIterator & kmer = nowVertex[kmerId];
				if(kmer.Valid())
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
								if(!Overlap(k, startKMer, prevData, nowData))
								{
									ret++;
									collapsed = true;								
									CollapseBulge(sequence, bifStorage, visit, k, startKMer, prevData, nowData);
									travelRange[nowData.kmerId] = prevData.distance;
									endChar[nowData.kmerId] = endChar[prevData.kmerId];
								}
								else
								{
									nowVertex[prevData.kmerId] = nowVertex[nowData.kmerId] = sequence.PositiveRightEnd();
								}

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