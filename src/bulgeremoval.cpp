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
			std::vector<std::pair<size_t, size_t> > & lookBack)
		{
			lookBack.clear();
			lookForward.clear();
			StrandIterator amer = AdvanceForward(*startKMer[targetData.kmerId], k).Invert();
			StrandIterator bmer = AdvanceForward(*startKMer[targetData.kmerId], targetData.distance);
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

	void BlockFinder::UpdateBifurcations(DNASequence & sequence,
			BifurcationStorage & bifStorage,
			size_t k,
			const IteratorProxyVector & startKMer,
			VisitData sourceData,
			VisitData targetData,
			const std::vector<std::pair<size_t, size_t> > & lookForward,
			const std::vector<std::pair<size_t, size_t> > & lookBack)
	{
		size_t anear = 0;
		size_t bnear = 0;
		StrandIterator amer = AdvanceForward(*startKMer[targetData.kmerId], k).Invert();
		StrandIterator bmer = AdvanceForward(*startKMer[targetData.kmerId], sourceData.distance);
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

		amer = *startKMer[targetData.kmerId];
		bmer = AdvanceForward(*startKMer[targetData.kmerId], sourceData.distance + k).Invert();
		StrandIterator srcAMer = *startKMer[sourceData.kmerId];
		StrandIterator srcBMer = AdvanceForward(*startKMer[sourceData.kmerId], sourceData.distance + k).Invert();
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
		BlockFinder::PrintRaw(sequence, std::cerr);
		std::cerr << "Source branch: " << std::endl;
		BlockFinder::PrintPath(sequence, *startKMer[sourceData.kmerId], k, sourceData.distance, std::cerr);
		std::cerr << "Target branch: " << std::endl;
		BlockFinder::PrintPath(sequence, *startKMer[targetData.kmerId], k, targetData.distance, std::cerr);
		bifStorage.Dump(sequence, k, std::cerr);
		iseq_->Test();
	#endif
		std::vector<std::pair<size_t, size_t> > lookForward;
		std::vector<std::pair<size_t, size_t> > lookBack;
		EraseBifurcations(sequence, bifStorage, k, startKMer, targetData, lookForward, lookBack);
		StrandIterator sourceIt = *startKMer[sourceData.kmerId];
		StrandIterator targetIt = *startKMer[targetData.kmerId];
		sequence.Replace(AdvanceForward(sourceIt, k),
			sourceData.distance,
			AdvanceForward(targetIt, k),
			targetData.distance,
			boost::bind(&BifurcationStorage::NotifyBefore, boost::ref(bifStorage), _1, _2),
			boost::bind(&BifurcationStorage::NotifyAfter, boost::ref(bifStorage), _1, _2));
		UpdateBifurcations(sequence, bifStorage, k, startKMer, sourceData, targetData, lookForward, lookBack);

	#ifdef _DEBUG
		std::cerr << "After: " << std::endl;
		BlockFinder::PrintRaw(sequence, std::cerr);
		std::cerr << "Source branch: " << std::endl;
		BlockFinder::PrintPath(sequence, *startKMer[sourceData.kmerId], k, sourceData.distance, std::cerr);
		std::cerr << "Target branch: " << std::endl;
		BlockFinder::PrintPath(sequence, *startKMer[targetData.kmerId], k, sourceData.distance, std::cerr);
		bifStorage.Dump(sequence, k, std::cerr);
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
				if(bifId != BifurcationStorage::NO_BIFURCATION)
				{
					if(bifId == startBifId && i > 0)
					{
						break;
					}

					bif[bifId] = i;
				}
			}
		}

		struct BranchVisitData
		{							
			std::vector<size_t> firstBranchId;
			std::vector<size_t> secondBranchId;
			std::vector<size_t> firstBranchLength;
			std::vector<size_t> secondBranchLength;
		};

		const size_t MLP_CONSTANT = 6;

		size_t GetMaxRange(size_t minBranchSize)
		{
			return minBranchSize * MLP_CONSTANT;
		}

		std::string GetString(StrandIterator it, size_t range)
		{
			return std::string(it, AdvanceForward(it, range));
		}

	}

	void FindSuperBulges(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, size_t minBranchSize, size_t bifId, std::vector<SuperBulge> & ret)
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
		boost::unordered_map<size_t, BranchVisitData> visit;
		for(size_t i = 0; i < startKMer.size(); i++)
		{				
			for(size_t j = i + 1; j < startKMer.size(); j++)
			{
				size_t counterJ = 0;
				size_t counterI = 0;				
				StrandIterator it = *startKMer[i];
				StrandIterator jt = *startKMer[j];
				bool fail = false;
				while(counterI < maxRange && counterJ < maxRange && it.AtValidPosition() && jt.AtValidPosition() && !fail)
				{
					while(counterI < maxRange && counterJ < maxRange && it.AtValidPosition() && jt.AtValidPosition() && *it == *jt)
					{
						size_t ibif = bifStorage.GetBifurcation(it);
						size_t jbif = bifStorage.GetBifurcation(jt);
						if((ibif == bifId && counterI > 0) || (jbif == bifId && counterJ > 0))
						{
							fail = true;
							break;
						}

						if(ibif != BifurcationStorage::NO_BIFURCATION && ibif == jbif && counterI > 0 && counterJ > 0)
						{
							visit[ibif].firstBranchId.push_back(i);
							visit[ibif].secondBranchId.push_back(j);
							visit[ibif].firstBranchLength.push_back(counterI);
							visit[ibif].secondBranchLength.push_back(counterJ);
							std::cout << 1 << std::endl << GetString(*startKMer[i], counterI + k) << std::endl << GetString(*startKMer[j], counterJ + k) << std::endl;
						}

						++it;
						++jt;
						counterI++;
						counterJ++;
						char ch1 = *it;
						char ch2 = *jt;
					}

					if(fail)
					{
						break;
					}
						
					std::map<size_t, size_t> nextIBif;
					size_t nextCommonBif = BifurcationStorage::NO_BIFURCATION;
					CollectBifurcations(bifStorage, it, minBranchSize, maxRange - counterI, nextIBif, bifId);
					for(size_t jdist = 0; jdist < minBranchSize && jdist < maxRange - counterJ && jt.AtValidPosition(); ++jdist, ++jt)
					{
						size_t jbif = bifStorage.GetBifurcation(jt);
						if(jbif == bifId && jdist > 0)
						{
							break;
						}

						if(jbif != BifurcationStorage::NO_BIFURCATION && nextIBif.find(jbif) != nextIBif.end())
						{
							counterJ += jdist;
							counterI += nextIBif[bifId];
							it = AdvanceForward(it, nextIBif[bifId]);
							visit[jbif].firstBranchId.push_back(i);
							visit[jbif].secondBranchId.push_back(j);
							visit[jbif].firstBranchLength.push_back(counterI);
							visit[jbif].secondBranchLength.push_back(counterJ);
							std::cout << 2 << std::endl << GetString(*startKMer[i], counterI + k) << std::endl << GetString(*startKMer[j], counterJ + k) << std::endl;
						}
					}
				}
			}				
		}

		SuperBulge best;
		best.score = -1;
		for(boost::unordered_map<size_t, BranchVisitData>::iterator point = visit.begin(); point != visit.end(); ++point)
		{
			std::map<size_t, size_t> pathLength;
			for(size_t path = 0; path < point->second.firstBranchId.size(); ++path)
			{
				size_t branch[] = {point->second.firstBranchId[path], point->second.secondBranchId[path]};
				size_t length[] = {point->second.firstBranchLength[path], point->second.secondBranchLength[path]};
				VisitData idata(branch[0], length[0]);
				VisitData jdata(branch[1], length[1]);
				if(!Overlap(k, startKMer, idata, jdata))
				{
					for(size_t i = 0; i < 2; i++)
					{
						if(pathLength.count(branch[i]) == 0)
						{
							pathLength[branch[i]] = length[i];
						}
						else
						{
							pathLength[branch[i]] = std::max(length[i], pathLength[branch[i]]);
						}
					}
				}
			}

			size_t sum = 0;
			std::set<char> startCharSet;
			std::vector<size_t> branches;
			for(std::map<size_t, size_t>::iterator it = pathLength.begin(); it != pathLength.end(); ++it)
			{
				sum += it->second;
				startCharSet.insert(endChar[it->first]);
				branches.push_back(it->first);
			}

			sum *= pathLength.size() * pathLength.size();
			if((sum > best.score || best.score == -1) && startCharSet.size() > 1)
			{
				best = SuperBulge(sum, bifId, point->first, branches);
			}
		}

		if(best.score != -1)
		{
			ret.push_back(best);
		}
	}		

	void BlockFinder::SimplifySuperBulge(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, size_t minBranchSize, SuperBulge bulge, std::set<size_t> & deprecateId)
	{			
		IteratorProxyVector startKMer;
		if(bifStorage.ListPositions(bulge.startId, std::back_inserter(startKMer)) < bulge.branch.size())
		{
			return;
		}

		std::vector<size_t> toDeprecate;
		size_t maxRange = GetMaxRange(minBranchSize);
		std::vector<size_t> range(bulge.branch.size());
		std::vector<size_t> maxBifDegree(bulge.branch.size());
		for(size_t i = 0; i < bulge.branch.size(); i++)
		{
			bool found = false;
			StrandIterator it = *startKMer[bulge.branch[i]];
			for(size_t j = 0; j < maxRange; ++j, ++it)
			{
				if(!it.AtValidPosition())
				{
					return;
				}

				size_t bifId = bifStorage.GetBifurcation(it);
				if(bifId != BifurcationStorage::NO_BIFURCATION)
				{
					if(deprecateId.count(bifId))
					{
						return;
					}

					toDeprecate.push_back(bifId);
					maxBifDegree[i] = std::max(maxBifDegree[i], bifStorage.CountBifurcations(bifId));
				}

				if(bifId == bulge.endId)
				{
					found = true;
					range[i] = j;
					break;						
				}
			}
			 
			if(!found)
			{
				return;
			}
		}
			
		size_t sample = std::max_element(maxBifDegree.begin(), maxBifDegree.end()) - maxBifDegree.begin();
		VisitData sampleData(bulge.branch[sample], range[sample]);
		StrandIterator it = *startKMer[bulge.branch[sample]];
		std::string sampleString(it, AdvanceForward(it, range[sample]));
		for(size_t i = 0; i < bulge.branch.size(); i++)
		{
			VisitData idata(bulge.branch[i], range[i]);
			StrandIterator it = *startKMer[bulge.branch[i]];
			std::string nowString(it, AdvanceForward(it, range[i]));
			if(i != sample)
			{				
				CollapseBulgeGreedily(sequence, bifStorage, k, startKMer, sampleData, idata);
			}
		}

		deprecateId.insert(toDeprecate.begin(), toDeprecate.end());
	}

	size_t BlockFinder::SimplifyGraphEasily(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, size_t minBranchSize, size_t maxIterations, ProgressCallBack callBack)
	{
		size_t count = 0;
		size_t iterations = 0;
		size_t totalProgress = 0;
		bool anyChanges = true;
		if(!callBack.empty())
		{
			callBack(totalProgress, start);
		}

		size_t threshold = (bifStorage.GetMaxId() * maxIterations * 2) / PROGRESS_STRIDE;
		do
		{
			iterations++;
			std::vector<SuperBulge> superBulge;
			for(size_t id = 0; id <= bifStorage.GetMaxId(); id++)
			{
				FindSuperBulges(sequence, bifStorage, k, minBranchSize, id, superBulge);
				if(++count >= threshold && !callBack.empty())
				{
					count = 0;
					totalProgress = std::min(totalProgress + 1, PROGRESS_STRIDE);
					callBack(totalProgress, run);
				}
			}

			std::set<size_t> deprecateId;
			std::sort(superBulge.begin(), superBulge.end());
			for(size_t i = 0; i < superBulge.size(); i++)
			{
				if(deprecateId.count(superBulge[i].startId) == 0)
				{
					SimplifySuperBulge(sequence, bifStorage, k, minBranchSize, superBulge[i], deprecateId);
					break;
				}
			}
		}
		while(iterations < maxIterations);

		if(!callBack.empty())
		{
			callBack(PROGRESS_STRIDE, end);
		}
		
		return 0;
	}
}
