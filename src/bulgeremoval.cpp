//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "blockfinder.h"

namespace SyntenyFinder
{
	namespace
	{
		const char EMPTY = ' ';

		bool CmpSizePair(const std::pair<size_t, size_t> & a, const std::pair<size_t, size_t> & b)
		{
			return a.first == b.first && a.second == b.second;
		}

		struct BifurcationMark
		{
			size_t bifId;
			size_t distance;
			BifurcationMark() {}
			BifurcationMark(size_t bifId, size_t distance): bifId(bifId),
				distance(distance) {}

			bool operator < (const BifurcationMark & compare) const
			{
				if(bifId != compare.bifId)
				{
					return bifId < compare.bifId;
				}

				return distance < compare.distance;
			}
		};

		size_t MaxBifurcationMultiplicity(const BifurcationStorage & bifStorage,
			StrandIterator it, size_t distance)
		{
			size_t ret = 0;
			for(size_t i = 0; i < distance - 1; i++)
			{
				size_t bifId = bifStorage.GetBifurcation(++it);
				if(bifId != BifurcationStorage::NO_BIFURCATION)
				{
					ret = std::max(ret, bifStorage.CountBifurcations(bifId));
				}
			}

			return ret;
		}

		void EraseBifurcations(DNASequence & sequence,
			BifurcationStorage & bifStorage,
			size_t k,
			const IteratorProxyVector  & startKMer,
			VisitData targetData,
			std::vector<std::pair<size_t, size_t> > & lookForward,
			std::vector<std::pair<size_t, size_t> > & lookBackward)
		{			
			lookForward.clear();
			lookBackward.clear();
			StrandIterator amer = AdvanceForward(*startKMer[targetData.kmerId], k).Invert();
			StrandIterator bmer = AdvanceForward(*startKMer[targetData.kmerId], targetData.distance);
			for(size_t i = 0; i < k; i++, ++amer, ++bmer)
			{
				size_t bifId = bifStorage.GetBifurcation(amer);
				if(bifId != BifurcationStorage::NO_BIFURCATION)
				{
					bifStorage.ErasePoint(amer);
					lookBackward.push_back(std::make_pair(i, bifId));
				}

				bifId = bifStorage.GetBifurcation(bmer);
				if(bifId != BifurcationStorage::NO_BIFURCATION)
				{
					bifStorage.ErasePoint(bmer);
					lookForward.push_back(std::make_pair(i, bifId));
				}
			}

			amer = *startKMer[targetData.kmerId];
			bmer = AdvanceForward(amer, k + targetData.distance).Invert();
			for(size_t i = 0; i < k + targetData.distance; i++, ++amer, ++bmer)
			{
				if(i > 0)
				{
					bifStorage.ErasePoint(amer);
				}

				bifStorage.ErasePoint(bmer);
			}
		}

		bool Overlap(size_t k,
			const IteratorProxyVector & startKMer,
			VisitData sourceData,
			VisitData targetData)
		{
			std::vector<size_t> occur;
			StrandIterator it = *startKMer[sourceData.kmerId];
			for(size_t i = 0; i < sourceData.distance + k; i++, ++it)
			{
				occur.push_back(it.GetElementId());
			}

			it = *startKMer[targetData.kmerId];
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

		void FillVisit(const DNASequence & sequence,
			const BifurcationStorage & bifStorage,
			StrandIterator kmer,
			size_t minBranchSize,
			std::vector<BifurcationMark> & visit)
		{
			visit.clear();
			size_t start = bifStorage.GetBifurcation(kmer);
			++kmer;
			for(size_t step = 1; step < minBranchSize && kmer.AtValidPosition(); ++kmer, step++)
			{
				size_t bifId = bifStorage.GetBifurcation(kmer);
				if(bifId == start)
				{
					break;
				}

				if(bifId != BifurcationStorage::NO_BIFURCATION)
				{
					visit.push_back(BifurcationMark(bifId, step));
				}
			}

			std::sort(visit.begin(), visit.end());
		}
		
		struct BranchData
		{
			BranchData() {}
			BranchData(char ch): endChar(ch) {}
			char endChar;
			std::vector<size_t> branchIds;
		};

		typedef std::vector< std::vector<size_t> > BulgedBranches;

		bool AnyBulges(DNASequence & sequence,
			BifurcationStorage & bifStorage,
			size_t k,
			const IteratorProxyVector & startKMer,
			const std::vector<char> & endChar,
			BulgedBranches & bulges,
			size_t minBranchSize)
		{
			//bool result = false;
			bulges.clear();
			boost::unordered_map<size_t, BranchData> visit;
			for(size_t i = 0; i < startKMer.size(); i++)
			{
				if(endChar[i] != EMPTY)
				{
					StrandIterator kmer = *startKMer[i];
					size_t start = bifStorage.GetBifurcation(kmer);
					++kmer;
					for(size_t step = 1; step < minBranchSize && kmer.AtValidPosition(); ++kmer, step++)
					{
						size_t bifId = bifStorage.GetBifurcation(kmer);
						if(bifId == start)
						{
							break;
						}

						if(bifId != BifurcationStorage::NO_BIFURCATION)
						{
							boost::unordered_map<size_t, BranchData>::iterator kt = visit.find(bifId);
							if(kt == visit.end())
							{
								BranchData bData(endChar[i]);
								bData.branchIds.push_back(i);
								visit[bifId] = bData;//endChar[i];
							}
							else if(kt->second.endChar != endChar[i])
							{
								kt->second.branchIds.push_back(i);
								break;
							}
						}
					}
				}
			}

			for (boost::unordered_map<size_t, BranchData>::iterator kt = visit.begin(); kt !=visit.end(); ++kt)
			{
				if (kt->second.branchIds.size() > 1)
				{
					bulges.push_back(kt->second.branchIds);
					/*std::cerr << "[";
					for (size_t i = 0; i < kt->second.branchIds.size(); ++i)
					{
						std::cerr << kt->second.branchIds[i] << ",";
					}
					std::cerr << "]";*/
				}
			}
			//std::cerr << std::endl;
			return !bulges.empty();
		}
	}

	void BlockFinder::ScanBifurcationsScale(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, const IteratorProxyVector & startKMer, VisitData sourceData,
		VisitData copyData, std::vector<std::pair<size_t, size_t> > & lookForward, std::vector<std::pair<size_t, size_t> > & lookBackward)
	{
		size_t bifId;
		lookForward.clear();
		lookBackward.clear();
		StrandIterator srcAMer = ++(*startKMer[sourceData.kmerId]);
		StrandIterator srcBMer = ++(AdvanceForward(*startKMer[sourceData.kmerId], sourceData.distance + k).Invert());
		for(size_t i = 1; i < sourceData.distance; i++, ++srcAMer, ++srcBMer)
		{
			bifId = bifStorage.GetBifurcation(srcAMer);
			if(bifId != BifurcationStorage::NO_BIFURCATION)
			{
				size_t pos = i * copyData.distance / sourceData.distance;
				lookForward.push_back(std::make_pair(i, bifId));
			}

			bifId = bifStorage.GetBifurcation(srcBMer);
			if(bifId != BifurcationStorage::NO_BIFURCATION)
			{
				size_t pos = i * copyData.distance / sourceData.distance;
				lookBackward.push_back(std::make_pair(pos, bifId));
			}
		}
	}

	void BlockFinder::ScanBifurcations(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, const IteratorProxyVector & startKMer, VisitData sourceData,
		std::vector<std::pair<size_t, size_t> > & lookForward, std::vector<std::pair<size_t, size_t> > & lookBackward)
	{
		size_t bifId;
		lookForward.clear();
		lookBackward.clear();
		StrandIterator srcAMer = *startKMer[sourceData.kmerId];
		StrandIterator srcBMer = AdvanceForward(*startKMer[sourceData.kmerId], sourceData.distance + k).Invert();
		for(size_t i = 0; i < sourceData.distance + 1; i++, ++srcAMer, ++srcBMer)
		{
			bifId = bifStorage.GetBifurcation(srcAMer);
			if(bifId != BifurcationStorage::NO_BIFURCATION)
			{
				lookForward.push_back(std::make_pair(i, bifId));
			}

			bifId = bifStorage.GetBifurcation(srcBMer);
			if(bifId != BifurcationStorage::NO_BIFURCATION)
			{
				lookBackward.push_back(std::make_pair(i, bifId));
			}
		}
	}

	void BlockFinder::SpellBulges(const DNASequence & sequence, size_t k,
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
			BlockFinder::PrintPath(sequence, startKMer[visitData[i].kmerId], k, visitData[i].distance, std::cerr);
		}

		std::cerr << DELIMITER << std::endl;
	}

	void BlockFinder::RestoreCornerBifurcations(DNASequence & sequence,
			BifurcationStorage & bifStorage,
			size_t k,
			const IteratorProxyVector & startKMer,
			VisitData sourceData,
			VisitData targetData,
			const std::vector<std::pair<size_t, size_t> > & lookForward,
			const std::vector<std::pair<size_t, size_t> > & lookBackward)
	{
		size_t anear = 0;
		size_t bnear = 0;
		StrandIterator amer = AdvanceForward(*startKMer[targetData.kmerId], k).Invert();
		StrandIterator bmer = AdvanceForward(*startKMer[targetData.kmerId], sourceData.distance);
		for(size_t i = 0; i < k; i++, ++amer, ++bmer)
		{
			if(anear < lookBackward.size() && i == lookBackward[anear].first)
			{
				bifStorage.AddPoint(amer, lookBackward[anear++].second);
			}

			if(bnear < lookForward.size() && i == lookForward[bnear].first)
			{
				bifStorage.AddPoint(bmer, lookForward[bnear++].second);
			}
		}
	}

	void BlockFinder::CollapseBulgeGreedily(DNASequence & sequence,
		BifurcationStorage & bifStorage,
		size_t k,
		IteratorProxyVector & startKMer,
		VisitData sourceData,
		VisitData targetData)
	{
	#ifdef _DEBUG
		static size_t bulge = 0;
		std::cerr << "Bulge #" << bulge++ << std::endl;
		std::cerr << "Before: " << std::endl;
//		BlockFinder::PrintRaw(sequence, std::cerr);
		std::cerr << "Source branch: " << std::endl;
		BlockFinder::PrintPath(sequence, *startKMer[sourceData.kmerId], k, sourceData.distance, std::cerr);
		std::cerr << "Target branch: " << std::endl;
		BlockFinder::PrintPath(sequence, *startKMer[targetData.kmerId], k, targetData.distance, std::cerr);
//		bifStorage.Dump(sequence, k, std::cerr);
		iseq_->Test();
	#endif
		std::vector<std::pair<size_t, size_t> > lookForward;
		std::vector<std::pair<size_t, size_t> > lookBackward;
		std::vector<std::pair<size_t, size_t> > lookForwardTarget;
		std::vector<std::pair<size_t, size_t> > lookBackwardTarget;
		EraseBifurcations(sequence, bifStorage, k, startKMer, targetData, lookForward, lookBackward);
		ScanBifurcationsScale(sequence, bifStorage, k, startKMer, targetData, sourceData, lookForwardTarget, lookBackwardTarget);
		StrandIterator sourceIt = *startKMer[sourceData.kmerId];
		StrandIterator targetIt = *startKMer[targetData.kmerId];
		sequence.Replace(AdvanceForward(sourceIt, k),
			sourceData.distance,
			AdvanceForward(targetIt, k),
			targetData.distance,
			boost::bind(&BifurcationStorage::NotifyBefore, boost::ref(bifStorage), _1, _2),
			boost::bind(&BifurcationStorage::NotifyAfter, boost::ref(bifStorage), _1, _2));
		RestoreCornerBifurcations(sequence, bifStorage, k, startKMer, sourceData, targetData, lookForward, lookBackward);
		ScanBifurcations(sequence, bifStorage, k, startKMer, sourceData, lookForward, lookBackward);
		VisitData newTargetData(targetData.kmerId, sourceData.distance);
		RestoreMainBifurcations(sequence, bifStorage, k, startKMer, sourceData, targetData, lookForward, lookBackward);
		RestoreMainBifurcations(sequence, bifStorage, k, startKMer, sourceData, targetData, lookForwardTarget, lookBackwardTarget);
		RestoreMainBifurcations(sequence, bifStorage, k, startKMer, sourceData, newTargetData, lookForwardTarget, lookBackwardTarget);
	#ifdef _DEBUG
		std::cerr << "After: " << std::endl;
//		BlockFinder::PrintRaw(sequence, std::cerr);
		std::cerr << "Source branch: " << std::endl;
		BlockFinder::PrintPath(sequence, *startKMer[sourceData.kmerId], k, sourceData.distance, std::cerr);
		std::cerr << "Target branch: " << std::endl;
		BlockFinder::PrintPath(sequence, *startKMer[targetData.kmerId], k, sourceData.distance, std::cerr);
//		bifStorage.Dump(sequence, k, std::cerr);
		iseq_->Test();
		std::cerr << DELIMITER << std::endl;
	#endif
	}


	size_t BlockFinder::RemoveBulges(DNASequence & sequence,
		BifurcationStorage & bifStorage, size_t k, size_t minBranchSize, size_t bifId)
	{
		size_t ret = 0;
		IteratorProxyVector startKMer;
		if(bifStorage.ListPositions(bifId, std::back_inserter(startKMer)) < 2)
		{
			return ret;
		}

		std::vector<char> endChar(startKMer.size(), EMPTY);
		for(size_t i = 0; i < startKMer.size(); i++)
		{
			if(ProperKMer(*startKMer[i], k + 1))
			{
				endChar[i] = *AdvanceForward(*startKMer[i], k);
			}
		}

		//std::vector<bool> isBulge(startKMer.size(), false);
		BulgedBranches bulges;
		if(!AnyBulges(sequence, bifStorage, k, startKMer, endChar, bulges, minBranchSize))
		{
			return ret;
		}

		std::vector<BifurcationMark> visit;
		for (size_t numBulge = 0; numBulge < bulges.size(); ++numBulge)
		{

			//for(size_t kmerI = 0; kmerI < startKMer.size(); kmerI++)
			for (size_t  idI = 0; idI < bulges[numBulge].size(); ++idI)
			{
				size_t kmerI = bulges[numBulge][idI];
				//if(!startKMer[kmerI].Valid() || !isBulge[kmerI])
				if(!startKMer[kmerI].Valid())
				{
					//if (!isBulge[kmerI]) std::cerr << "optimised!\n";
					continue;
				}

				FillVisit(sequence, bifStorage, *startKMer[kmerI], minBranchSize, visit);
				//for(size_t kmerJ = kmerI + 1; kmerJ < startKMer.size(); kmerJ++)
				for(size_t  idJ = idI + 1; idJ < bulges[numBulge].size(); ++idJ)
				{
					size_t kmerJ = bulges[numBulge][idJ];
					//if(!startKMer[kmerJ].Valid() || endChar[kmerI] == endChar[kmerJ] || !isBulge[kmerJ])
					if(!startKMer[kmerJ].Valid() || endChar[kmerI] == endChar[kmerJ])
					{
						//if (!isBulge[kmerJ]) std::cerr << "optimised!\n";
						continue;
					}

					StrandIterator kmer = ++StrandIterator(*startKMer[kmerJ]);
					for(size_t step = 1; kmer.AtValidPosition() && step < minBranchSize; ++kmer, step++)
					{
						size_t nowBif = bifStorage.GetBifurcation(kmer);
						if(nowBif != BifurcationStorage::NO_BIFURCATION)
						{
							if(nowBif == bifId)
							{
								break;
							}

							std::vector<BifurcationMark>::iterator vt = std::lower_bound(visit.begin(), visit.end(), BifurcationMark(nowBif, 0));
							if(vt != visit.end() && vt->bifId == nowBif)
							{
								VisitData jdata(kmerJ, step);
								VisitData idata(kmerI, vt->distance);
								if(Overlap(k, startKMer, idata, jdata) || nowBif == bifId)
								{
									break;
								}

								++ret;
								size_t imlp = MaxBifurcationMultiplicity(bifStorage, *startKMer[kmerI], idata.distance);
								size_t jmlp = MaxBifurcationMultiplicity(bifStorage, *startKMer[kmerJ], jdata.distance);
								bool iless = imlp > jmlp || (imlp == jmlp && idata.kmerId < jdata.kmerId);
								{
									if(iless)
									{
										endChar[jdata.kmerId] = endChar[idata.kmerId];
										CollapseBulgeGreedily(sequence, bifStorage, k, startKMer, idata, jdata);
									}
									else
									{
										endChar[idata.kmerId] = endChar[jdata.kmerId];
										CollapseBulgeGreedily(sequence, bifStorage, k, startKMer, jdata, idata);
										FillVisit(sequence, bifStorage, *startKMer[kmerI], minBranchSize, visit);
									}

									break;
								}
							}
						}
					}
				}
			}
		}

		bifStorage.Cleanup();
		return ret;
	}

	namespace
	{
		void CollectBifurcations(BifurcationStorage & bifStorage, StrandIterator it, size_t minBranchSize, size_t distLeft, std::map<size_t, size_t> & bif, size_t startBifId)
		{
			for(size_t i = 0; i < minBranchSize && i < distLeft && it.AtValidPosition(); i++, ++it)
			{
				size_t bifId = bifStorage.GetBifurcation(it);
				if(bifId != BifurcationStorage::NO_BIFURCATION && i > 0)
				{
					if(bifId == startBifId)
					{
						break;
					}

					bif[bifId] = i;
				}
			}
		}

		const size_t MLP_CONSTANT = 100;

		size_t GetMaxRange(size_t minBranchSize)
		{
			return minBranchSize * MLP_CONSTANT;
		}

		std::string GetString(StrandIterator it, size_t range)
		{
			return std::string(it, AdvanceForward(it, range));
		}

	}

	void FindSuperBulges(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, size_t minBranchSize, size_t minPathLength, size_t bifId, std::vector<SuperBulge> & ret)
	{
		IteratorProxyVector startKMer;
		if(bifStorage.ListPositions(bifId, std::back_inserter(startKMer)) < 2)
		{
			return;
		}

		std::vector<char> endChar(startKMer.size(), EMPTY);
		for(size_t i = 0; i < startKMer.size(); i++)
		{
			if(ProperKMer(*startKMer[i], k + 1))
			{
				endChar[i] = *AdvanceForward(*startKMer[i], k);
			}
		}

		size_t maxRange = GetMaxRange(minBranchSize);
		for(size_t i = 0; i < startKMer.size(); i++)
		{				
			for(size_t j = i + 1; j < startKMer.size(); j++)
			{
				bool fail = false;
				size_t counterJ = 0;
				size_t counterI = 0;
				StrandIterator it = *startKMer[i];
				StrandIterator jt = *startKMer[j];
				size_t nextCommonBif = BifurcationStorage::NO_BIFURCATION;
				while(counterI < maxRange && counterJ < maxRange && it.AtValidPosition() && jt.AtValidPosition() && !fail)
				{
					fail = true;
					std::map<size_t, size_t> nextIBif;
					std::map<size_t, size_t> nextJBif;
					CollectBifurcations(bifStorage, it, minBranchSize, maxRange - counterI, nextIBif, bifId);
					CollectBifurcations(bifStorage, jt, minBranchSize, maxRange - counterJ, nextJBif, bifId);
					for(std::map<size_t, size_t>::iterator bifIt = nextIBif.begin(); bifIt != nextIBif.end() && fail; ++bifIt)
					{
						if(nextJBif.count(bifIt->first) > 0)
						{
							fail = false;
							nextCommonBif = bifIt->first;
							counterI += nextIBif[nextCommonBif];
							counterJ += nextJBif[nextCommonBif];
							std::advance(it, nextIBif[nextCommonBif]);
							std::advance(jt, nextJBif[nextCommonBif]);
						}
					}
				}

				VisitData idata(i, counterI);
				VisitData jdata(j, counterJ);				
				if(nextCommonBif != BifurcationStorage::NO_BIFURCATION && !Overlap(k, startKMer, idata, jdata))
				{
					std::vector<std::string> branchSet;
					branchSet.push_back(GetString(*startKMer[i], counterI + k));
					branchSet.push_back(GetString(*startKMer[j], counterJ + k));					
					if(branchSet[0] != branchSet[1] && branchSet[0].size() + branchSet[1].size() >= 2 * minPathLength)
					{
						ret.push_back(SuperBulge(counterI + counterJ, bifId, nextCommonBif, idata, jdata, branchSet));
					}
				}
			}				
		}
	}	

	bool BlockFinder::SimplifySuperBulge(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, size_t minBranchSize, SuperBulge bulge, std::set<size_t> & deprecateId)
	{
		IteratorProxyVector startKMer;
		VisitData branch[] = {bulge.idata, bulge.jdata};
		std::vector<std::string> nowBranchSet;
		bifStorage.ListPositions(bulge.startId, std::back_inserter(startKMer));
		for(size_t i = 0; i < 2; i++)
		{
			if(branch[i].kmerId >= startKMer.size())
			{
				return false;
			}

			StrandIterator it = *startKMer[branch[i].kmerId];
			for(size_t j = 0; j < branch[i].distance; j++, ++it)
			{
				if(!it.AtValidPosition())
				{
					return false;
				}
			}

			nowBranchSet.push_back(GetString(*startKMer[branch[i].kmerId], branch[i].distance + k));
			if(bifStorage.GetBifurcation(it) != bulge.endId)
			{
				return false;
			}
		}

		if(Overlap(k, startKMer, bulge.idata, bulge.jdata) || bulge.branchSet != nowBranchSet)
		{						
			return false;
		}

		std::vector<size_t> degree(2);
		for(size_t i = 0; i < degree.size(); i++)
		{
			degree[i] = MaxBifurcationMultiplicity(bifStorage, *startKMer[branch[i].kmerId], branch[i].distance);
		}

		size_t sample = std::max_element(degree.begin(), degree.end()) - degree.begin();
		VisitData sampleData(branch[sample]);
		StrandIterator it = *startKMer[branch[sample].kmerId];
		std::string sampleString(it, AdvanceForward(it, branch[sample].distance + k));
		std::cerr << "Found a superbulge #" << bulgeId++ << std::endl;
		std::cerr << "Branch set:" << std::endl;
		std::copy(bulge.branchSet.begin(), bulge.branchSet.end(), std::ostream_iterator<std::string>(std::cerr, "\n"));
		for(size_t i = 0; i < 2; i++)
		{
			VisitData idata(branch[i]);
			StrandIterator it = *startKMer[branch[i].kmerId];
			if(i != sample)
			{
				CollapseBulgeGreedily(sequence, bifStorage, k, startKMer, sampleData, idata);
			}
		}

		return true;
	}

	void BlockFinder::RestoreMainBifurcations(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, const IteratorProxyVector & startKMer, VisitData sourceData, VisitData targetData,
			const std::vector<std::pair<size_t, size_t> > & lookForward, const std::vector<std::pair<size_t, size_t> > & lookBackward)
	{
		size_t anear = 0;
		size_t bnear = 0;
		StrandIterator amer = *startKMer[targetData.kmerId];
		StrandIterator bmer = AdvanceForward(*startKMer[targetData.kmerId], sourceData.distance + k).Invert();
		for(size_t i = 0; i < sourceData.distance + 1; i++, ++amer, ++bmer)
		{
			if(anear < lookForward.size() && i == lookForward[anear].first)
			{
				bifStorage.AddPoint(amer, lookForward[anear++].second);
			}

			if(bnear < lookBackward.size() && i == lookBackward[bnear].first)
			{
				bifStorage.AddPoint(bmer, lookBackward[bnear++].second);
			}
		}			
	}	

	size_t BlockFinder::SimplifyGraphEasily(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, size_t minBranchSize, size_t minPathSize, size_t maxIterations, ProgressCallBack callBack)
	{
		size_t count = 0;
		size_t iterations = 0;
		size_t totalProgress = 0;
		if(!callBack.empty())
		{
			callBack(totalProgress, start);
		}
				
		bulgeId = 0;
		bool simplify = false;
		size_t totalBulges = 0;
		size_t threshold = (bifStorage.GetMaxId() * maxIterations * 2) / PROGRESS_STRIDE;
		std::set<std::vector<std::string> > prevBulge;
		do
		{
			iterations++;
			std::vector<SuperBulge> superBulge;
			for(size_t id = 0; id <= bifStorage.GetMaxId(); id++)
			{
				FindSuperBulges(sequence, bifStorage, k, minBranchSize, minPathSize, id, superBulge);
				if(++count >= threshold && !callBack.empty())
				{
					count = 0;
					totalProgress = std::min(totalProgress + 1, PROGRESS_STRIDE);
					callBack(totalProgress, run);
				}
			}
			
			simplify = false;
			std::set<size_t> deprecateId;
			std::sort(superBulge.begin(), superBulge.end());
			for(size_t i = 0; i < superBulge.size(); i++)
			{
				if(SimplifySuperBulge(sequence, bifStorage, k, minBranchSize, superBulge[i], deprecateId))
				{
					totalBulges++;
					simplify = true;
					prevBulge.insert(superBulge[i].branchSet);
				}
			}
		}
		while(simplify && iterations < maxIterations);

		if(!callBack.empty())
		{
			callBack(PROGRESS_STRIDE, end);
		}
		
		return totalBulges;
	}
}
