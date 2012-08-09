#include "graphalgorithm.h"

namespace SyntenyBuilder
{
	namespace
	{
		struct BifurcationMark
		{
			size_t bifId;
			size_t distance;
			BifurcationMark() {}
			BifurcationMark(size_t bifId, size_t distance): bifId(bifId),
				distance(distance) {}

			bool operator < (const BifurcationMark & compare)
			{
				if(bifId != compare.bifId)
				{
					return bifId < compare.bifId;
				}

				return distance < compare.distance;
			}
		};

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
			StrandIterator amer = AdvanceForward(startKMer[targetData.kmerId], k).Invert();
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
			bmer = AdvanceForward(amer, k + targetData.distance).Invert();
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
			StrandIterator amer = AdvanceForward(startKMer[targetData.kmerId], k).Invert();
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
			bmer = AdvanceForward(startKMer[targetData.kmerId], sourceData.distance + k).Invert();
			StrandIterator srcAMer = startKMer[sourceData.kmerId];
			StrandIterator srcBMer = AdvanceForward(startKMer[sourceData.kmerId], sourceData.distance + k).Invert();
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
				occur.push_back(it.GetElementId());
			}

			it = startKMer[targetData.kmerId];
			std::sort(occur.begin(), occur.end());
			for(size_t i = 0; i < targetData.distance + k; i++, ++it)
			{
				if(std::binary_search(occur.begin(), occur.end(), it.GetElementId()))
				{
					return true;
				}
			}

			return false;
		}
		
		void CollapseBulge(DNASequence & sequence,
			BifurcationStorage & bifStorage,
			size_t k,
			std::vector<StrandIterator> & startKMer,
			const std::multimap<size_t, size_t> & restricted,
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

			StrandIterator it = startKMer[targetData.kmerId];
			for(size_t step = 0; step < targetData.distance + k; step++, ++it)
			{
				typedef std::multimap<size_t, size_t>::const_iterator MMIterator;
				std::pair<MMIterator, MMIterator> kt = restricted.equal_range(it.GetElementId());
				for(; kt.first != kt.second; ++kt.first)
				{
					if(kt.first->second != targetData.kmerId)
					{
						startKMer[kt.first->second] = sequence.PositiveEnd();
					}
				}
			}

			std::vector<std::pair<size_t, size_t> > lookForward;
			std::vector<std::pair<size_t, size_t> > lookBack;
			EraseBifurcations(sequence, bifStorage, k, startKMer, targetData, lookForward, lookBack);
			StrandIterator sourceIt = startKMer[sourceData.kmerId];
			StrandIterator targetIt = startKMer[targetData.kmerId];
			size_t diff = targetData.distance - sourceData.distance;
			sequence.CopyN(sourceIt, sourceData.distance, targetIt);
			targetIt = AdvanceForward(targetIt, sourceData.distance);
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

		void FillVisit(const BifurcationStorage & bifStorage, 
			StrandIterator kmer,
			size_t minBranchSize,
			std::vector<BifurcationMark> & visit)			
		{
			visit.clear();
			for(size_t step = 1; step < minBranchSize; step++)
			{
				size_t bifId = bifStorage.GetBifurcation(++kmer);
				if(bifId != BifurcationStorage::NO_BIFURCATION)
				{
					visit.push_back(BifurcationMark(bifId, step));
				}
			}

			std::sort(visit.begin(), visit.end());
		}

		void SpellBulges(const DNASequence & sequence, size_t k,
			size_t bifStart,
			size_t bifEnd,
			const std::vector<StrandIterator> & startKMer,
			const std::vector<VisitData> & visitData)
		{
			static size_t bulge = 0;
			std::cerr << "Bulge #" << bulge++ << ", " << "(" << bifStart << ", " << bifEnd << ")" << std::endl;
			for(size_t i = 0; i < visitData.size(); i++)
			{
				std::cerr << "Branch #" << i << ", size = " << visitData[i].distance + k << ":" << std::endl;			
				GraphAlgorithm::PrintPath(startKMer[visitData[i].kmerId], k, visitData[i].distance, std::cerr);
			}

			std::cerr << DELIMITER << std::endl;
		}
	}

	/*
	size_t GraphAlgorithm::FindBulges(const DNASequence & sequence, const BifurcationStorage & bifStorage,
		size_t k, size_t bifId)
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
		std::set<std::vector<size_t> > bulgeVisit;
		std::vector<char> endChar(startKMer.size(), ' ');
		for(size_t i = 0; i < startKMer.size(); i++)
		{
			if(startKMer[i].ProperKMer(k + 1))
			{                    
				nowVertex[i] = AdvanceForward(startKMer[i], 1);
				endChar[i] = *AdvanceForward(startKMer[i], k);
			}
		}

		bool goOn = true;
		std::vector<size_t> travelRange(startKMer.size(), 1);
		for(size_t step = 0; goOn; step++)
		{
			goOn = false;
			for(size_t kmerId = 0; kmerId < nowVertex.size(); kmerId++)
			{
				StrandIterator & kmer = nowVertex[kmerId];
				if(kmer.Valid())
				{
					goOn = true;
					size_t bifurcation = bifStorage.GetBifurcation(kmer);
					VisitData nowData = VisitData(kmerId, travelRange[kmerId]);
					if(bifurcation != BifurcationStorage::NO_BIFURCATION)
					{
						if(bifurcation == bifId)
						{
							nowVertex[kmerId] = sequence.PositiveEnd();
							continue;
						}						
					}
						
					++kmer;
					travelRange[kmerId]++;
					if(bifurcation != BifurcationStorage::NO_BIFURCATION)
					{
						visit.insert(std::make_pair(bifurcation, nowData));
					}
				}
			}
		}

		for(VertexVisitMap::iterator it = visit.begin(); it != visit.end(); )
		{	
			std::set<char> ends;
			std::vector<size_t> bulge;
			std::vector<VisitData> branch;
			VertexVisitMap::iterator jt = it;
			for(; jt != visit.end() && jt->first == it->first; jt++)
			{				
				branch.push_back(jt->second);
				bulge.push_back(jt->second.kmerId);
				ends.insert(endChar[jt->second.kmerId]);
			}

			std::sort(bulge.begin(), bulge.end());
			if(bulge.size() > 1 && bulgeVisit.count(bulge) == 0 && ends.size() > 1)
			{
				bulgeVisit.insert(bulge);
				SpellBulges(sequence, k, bifId, it->first, startKMer, branch);
			}

			it = jt;
		}

		return ret;
	}*/

	size_t GraphAlgorithm::RemoveBulges(DNASequence & sequence,
		BifurcationStorage & bifStorage, size_t k, size_t minBranchSize, size_t bifId)
	{	
		size_t ret = 0;		
		std::vector<StrandIterator> startKMer;
		std::multimap<size_t, size_t> restricted;
		if(bifStorage.ListPositions(bifId, std::back_inserter(startKMer)) < 2)
		{
			return ret;
		}

		std::vector<char> endChar(startKMer.size(), ' ');
		for(size_t i = 0; i < startKMer.size(); i++)
		{
			if(ProperKMer(startKMer[i], sequence, k + 1))
			{                    
				endChar[i] = *AdvanceForward(startKMer[i], k);
			}

			StrandIterator it = startKMer[i];
			for(size_t j = 0; j < k; j++, ++it)
			{
				restricted.insert(std::make_pair(it.GetElementId(), i));
			}
		}

		std::vector<BifurcationMark> visit;
		for(size_t kmerI = 0; kmerI < startKMer.size(); kmerI++)
		{
			if(!Valid(startKMer[kmerI], sequence))
			{
				continue;
			}

			FillVisit(bifStorage, startKMer[kmerI], minBranchSize, visit);
			for(size_t kmerJ = kmerI + 1; kmerJ < startKMer.size(); kmerJ++)
			{
				if(!Valid(startKMer[kmerJ], sequence) || endChar[kmerI] == endChar[kmerJ])
				{
					continue;
				}

				StrandIterator kmer = ++StrandIterator(startKMer[kmerJ]);
				for(size_t step = 1; Valid(kmer, sequence) && step < minBranchSize; ++kmer, step++)
				{
					size_t nowBif = bifStorage.GetBifurcation(kmer);
					if(nowBif != BifurcationStorage::NO_BIFURCATION)
					{
						std::vector<BifurcationMark>::iterator vt = std::lower_bound(
							visit.begin(), visit.end(), BifurcationMark(nowBif, 0));
						if(vt != visit.end() && vt->bifId == nowBif)
						{
							VisitData jdata(kmerJ, step);
							VisitData idata(kmerI, vt->distance);
							if(Overlap(k, startKMer, idata, jdata) || nowBif == bifId)
							{
								break;
							}

							++ret;
							bool iless = idata.distance < jdata.distance || (idata.distance == jdata.distance &&
								idata.kmerId < jdata.kmerId);
							if(iless)
							{
								endChar[jdata.kmerId] = endChar[idata.kmerId];
								CollapseBulge(sequence, bifStorage, k, startKMer, restricted, idata, jdata);
							}
							else
							{
								endChar[idata.kmerId] = endChar[jdata.kmerId];
								CollapseBulge(sequence, bifStorage, k, startKMer, restricted, jdata, idata);
								FillVisit(bifStorage, startKMer[kmerI], minBranchSize, visit);
							}

							break;
						}
					}
				}
			}
		}

		return ret;
	}
}