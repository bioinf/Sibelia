#include "graphalgorithm.h"
//#define _DEBUG
//#undef _DEBUG
namespace SyntenyBuilder
{
	namespace
	{
		typedef char Bool;	
		typedef long long int64;
		typedef DNASequence::StrandIterator StrandIterator;
		typedef google::sparse_hash_set<StrandIterator, KMerIndex::KMerHashFunction,
			KMerIndex::KMerEqualTo> KMerSet;	
		typedef std::pair<StrandIterator, StrandIterator> Vertex;
		typedef std::pair<StrandIterator, StrandIterator> IteratorPair;
		typedef std::vector<Vertex> VertexVector;

		struct VisitData
		{
			size_t kmerId;
			size_t distance;
			VisitData() {}
			VisitData(size_t kmerId, size_t distance): kmerId(kmerId), distance(distance) {}

		};

		class BifurcationStorage
		{
		public:
			static const size_t NO_BIFURCATION;

			void Clear()
			{
				for(size_t strand = 0; strand < 2; strand++)
				{
					posBifurcation[strand].clear();
					bifurcationPos[strand].clear();
				}
			}

			void AddPoint(StrandIterator it, size_t bifId)
			{
				size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
				BifMapIterator place = bifurcationPos[strand].insert(std::make_pair(bifId, it.GetPosition()));
				posBifurcation[strand].insert(place);
			}

			void ErasePoint(StrandIterator it)
			{
				size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
				BifMapIterator jt = temp.insert(std::make_pair(-1, it.GetPosition()));
				PosBifurcation::iterator kt = posBifurcation[strand].find(jt);
				if(kt != posBifurcation[strand].end())
				{
					bifurcationPos[strand].erase(*kt);
					posBifurcation[strand].erase(kt);
				}

				temp.clear();
			}

			size_t GetBifurcation(StrandIterator it) const
			{
				size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
				BifMapIterator jt = temp.insert(std::make_pair(-1, it.GetPosition()));
				PosBifurcation::const_iterator kt = posBifurcation[strand].find(jt);
				temp.clear();
				return kt == posBifurcation[strand].end() ? NO_BIFURCATION : (*kt)->first;
			}

			template<class Iterator>
				size_t ListPositions(size_t bifId, Iterator out, const DNASequence & seq) const
				{
					size_t ret = 0;
					typedef boost::function<StrandIterator (size_t)> Transformer;
					Transformer trans[2] = 
					{
						boost::bind(&DNASequence::PositiveByIndex, boost::cref(seq), _1),
						boost::bind(&DNASequence::NegativeByIndex, boost::cref(seq), _1)
					};

					for(size_t strand = 0; strand < 2; strand++)
					{
						std::pair<CBifMapIterator, CBifMapIterator> range = 
							bifurcationPos[strand].equal_range(bifId);
						for(CBifMapIterator it = range.first; it != range.second; ++it, ++ret)
						{
							*out++ = trans[strand](it->second);
						}
					}

					return ret;
				}

		private:
			typedef boost::unordered_multimap<size_t, size_t> BifurcationPos;
			typedef BifurcationPos::iterator BifMapIterator;
			typedef BifurcationPos::const_iterator CBifMapIterator;
			struct IteratorLess
			{
			public:
				bool operator () (BifMapIterator it1, BifMapIterator it2)
				{
					return it1->second < it2->second;
				}
			};

			typedef std::set<BifurcationPos::iterator, IteratorLess> PosBifurcation;

			mutable BifurcationPos temp;
			BifurcationPos bifurcationPos[2];
			PosBifurcation posBifurcation[2];
		};

		size_t deletedBulge;
		const size_t MOD = 100000;
		const size_t BifurcationStorage::NO_BIFURCATION = -1;

		typedef boost::unordered_multimap<size_t, VisitData> VertexVisitMap;
		
		void OutputEdge(const KMerIndex & index, StrandIterator it, std::ostream & out)
		{
			CopyN(it, index.GetK(), std::ostream_iterator<char>(out));
			out << " -> ";
			CopyN(++StrandIterator(it), index.GetK(), std::ostream_iterator<char>(out));

			char buf[1 << 8];
			if(it.GetDirection() == DNASequence::positive)
			{
				sprintf(&buf[0], "[color=\"%s\", label=\"%lu\"];", "blue", static_cast<long long unsigned>(it.GetPosition()));
			}
			else
			{
				sprintf(&buf[0], "[color=\"%s\", label=\"%lu\"];", "red", static_cast<long long unsigned>(it.GetPosition()));
			}

			out << " " << buf;
		}
		
		void ProcessIterator(KMerSet & visit, const KMerIndex & index, DNASequence::StrandIterator it, std::ostream & out)
		{
			std::vector<StrandIterator> kmer;		
			if(it.ProperKMer(index.GetK()) && visit.find(it) == visit.end())
			{
				visit.insert(it);				
				index.ListEquivalentKmers(it, kmer);
				for(size_t i = 0; i < kmer.size(); i++)
				{
					if(kmer[i].ProperKMer(index.GetK() + 1))
					{
						OutputEdge(index, kmer[i], out);
						out << std::endl;
					}
				}
			}
		}

		bool Less(const std::pair<size_t, StrandIterator> & p1, const std::pair<size_t, StrandIterator> & p2)
		{
			return p1.first < p2.first;
		}

		bool Invalid(const std::vector<Bool> & visit, const StrandIterator & g)
		{
			return visit[g.GetPosition()] == 1;
		}

		void Invalidate(std::vector<Bool> & visit, const StrandIterator & g)
		{
			visit[g.GetPosition()] = true;
		}
			
		size_t Extend(std::vector<Bool> & visit,
			std::vector<StrandIterator> kmer,
			boost::function<void (StrandIterator&)> extender,
			boost::function<void (const StrandIterator&)> invalidator)
		{
			size_t ret = 0;			
			std::for_each(kmer.begin(), kmer.end(), extender);			
			for(bool fail = false; !fail; ret++)
			{
				fail = !kmer[0].Valid() || visit[kmer[0].GetPosition()];
				if(!fail)
				{
					char consensus = *kmer[0];
					for(size_t i = 1; i < kmer.size() && !fail; i++)
					{
						fail = !kmer[i].Valid() || visit[kmer[i].GetPosition()] || *kmer[i] != consensus;
					}
				}

				if(!fail)
				{
					std::for_each(kmer.begin(), kmer.end(), invalidator);
					std::for_each(kmer.begin(), kmer.end(), extender);
				}
			}

 			return ret - 1;
		}

		void PrintRaw(DNASequence & s, std::ostream & out)
		{
			std::string rcomp;
			s.SpellRaw(std::ostream_iterator<char>(out));
			out << std::endl;
			for(size_t i = 0; i < s.Size(); i++)
			{
				out << i % 10;
			}
		
			out << std::endl;
			std::copy(s.PositiveBegin(), s.PositiveRightEnd(), std::ostream_iterator<char>(out));
			out << std::endl;
			std::copy(s.NegativeBegin(), s.NegativeRightEnd(), std::back_inserter(rcomp));
			std::copy(rcomp.rbegin(), rcomp.rend(), std::ostream_iterator<char>(out));
			out << std::endl;
		}

		void PrintPath(StrandIterator e, size_t k, size_t distance, std::ostream & out)
		{
			out << e.GetPosition() << " ";
			out << (e.GetDirection() == DNASequence::positive ? "s+ " : "s- ");
			CopyN(e, distance + k, std::ostream_iterator<char>(out));
			std::cerr << std::endl;
		}			

		void ClearVisit(VertexVisitMap & visit,
			BifurcationStorage & bifStorage,
			const std::vector<StrandIterator> & startVertex,
			VisitData targetData)
		{	
			StrandIterator target = startVertex[targetData.kmerId];
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
			const std::vector<StrandIterator> & startVertex,
			VisitData targetData
			)
		{
			IteratorPair target = std::make_pair(AdvanceForward(startVertex[targetData.kmerId], 1),
				AdvanceForward(startVertex[targetData.kmerId], targetData.distance));
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
			const std::vector<StrandIterator> & startVertex,
			VisitData sourceData,
			VisitData targetData
			)
		{
			IteratorPair target = std::make_pair(AdvanceForward(startVertex[targetData.kmerId], 1),
				AdvanceForward(startVertex[targetData.kmerId], sourceData.distance - 1));
			IteratorPair inverseTarget = InvertRange(target, k);
			IteratorPair source = std::make_pair(AdvanceForward(startVertex[sourceData.kmerId], 1),
				AdvanceForward(startVertex[sourceData.kmerId], sourceData.distance - 1));
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
			const std::vector<StrandIterator> & startVertex,
			VisitData sourceData,
			VisitData targetData
			)
		{
		#ifdef _DEBUG
			static size_t bulge = 0;
			std::cerr << "Bulge #" << bulge++ << std::endl;
			std::cerr << "Before: " << std::endl;
			PrintRaw(sequence, std::cerr);
			std::cerr << "Source branch: " << std::endl;			
			PrintPath(startVertex[sourceData.kmerId], k, sourceData.distance, std::cerr);
			std::cerr << "Target branch: " << std::endl;			
			PrintPath(startVertex[targetData.kmerId], k, targetData.distance, std::cerr);
		#endif

			ClearVisit(visit, bifStorage, startVertex, targetData);
			EraseBifurcations(sequence, bifStorage, k, startVertex, sourceData);
			StrandIterator sourceIt = startVertex[sourceData.kmerId];
			StrandIterator targetIt = startVertex[targetData.kmerId];
			sourceIt.Jump(k);
			targetIt.Jump(k);
			size_t diff = targetData.distance - sourceData.distance;
			deletedBulge += diff;
			sequence.CopyN(sourceIt, sourceData.distance, targetIt);
			targetIt.Jump(sourceData.distance);
			sequence.EraseN(targetIt, diff);
			UpdateBifurcations(sequence, bifStorage, k, startVertex, sourceData, targetData);

		#ifdef _DEBUG
			std::cerr << "After: " << std::endl;
			PrintRaw(sequence, std::cerr);
			std::cerr << "Source branch: " << std::endl;			
			PrintPath(startVertex[sourceData.kmerId], k, sourceData.distance, std::cerr);
			std::cerr << "Target branch: " << std::endl;			
			PrintPath(startVertex[targetData.kmerId], k, sourceData.distance, std::cerr);
			std::cerr << std::string(80, '-') << std::endl;
		#endif
		}	

		size_t RemoveBulges(BifurcationStorage & bifStorage, DNASequence & sequence, 
			size_t k, size_t minBranchSize, size_t bifId)
		{	
			size_t ret = 0;
			VertexVisitMap visit;
			VertexVector nowVertex;
			std::vector<StrandIterator> startVertex;
			if(bifStorage.ListPositions(bifId, std::back_inserter(startVertex), sequence) < 2)
			{
				return ret;
			}

			nowVertex.resize(startVertex.size());
			std::vector<char> endChar(startVertex.size(), ' ');
			for(size_t i = 0; i < startVertex.size(); i++)
			{
				if(startVertex[i].ProperKMer(k + 1))
				{                    
					nowVertex[i].first = nowVertex[i].second = startVertex[i];
					nowVertex[i].first.Jump(1);
					nowVertex[i].second.Jump(k);
					endChar[i] = *nowVertex[i].second;
				}
			}

			std::vector<size_t> travelRange(startVertex.size(), 1);
			for(size_t step = 0; step < minBranchSize; step++)
			{
				for(size_t kmerId = 0; kmerId < nowVertex.size(); kmerId++)
				{
					StrandIterator & kmerStart = nowVertex[kmerId].first;
					StrandIterator & kmerEnd = nowVertex[kmerId].second;
					if(kmerEnd.Valid())
					{
						bool collapsed = false;
						size_t bifurcation = bifStorage.GetBifurcation(kmerStart);
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
								if(endChar[kmerId] != endChar[it->second.kmerId])
								{
									ret++;
									collapsed = true;
									CollapseBulge(sequence, bifStorage, visit, k, startVertex, it->second, nowData);
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
							travelRange[kmerId]++;
							++kmerStart;
							++kmerEnd;
						}
					}
				}
			}

			return ret;
		}

		size_t EnumerateBifurcations(DNASequence & sequence, size_t k, BifurcationStorage & bifStorage)
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
				if(i == save[near].first)
				{
					bifStorage.AddPoint(++it, save[near].second);
				}
			}
		}

		size_t RemoveWhirls(BifurcationStorage & bifStorage, DNASequence & sequence, 
			size_t k, size_t minBranchSize)
		{
			size_t ret = 0;
			for(StrandIterator it = sequence.PositiveBegin(); it != sequence.PositiveRightEnd(); ++it)
			{
				size_t bifId = bifStorage.GetBifurcation(it);
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
								ret++;
								WhirlRemovalBifurcationUpdate(bifStorage, sequence, k, it, step);
								sequence.EraseN(AdvanceForward(it, k), step);
								remove = true;
								break;
							}
						}
					}
				}
			}
			
			return ret;
		}
	}
	
	
	void GraphAlgorithm::SerializeGraph(const DNASequence & sequence, size_t k, std::ostream & out)
	{
		out << "digraph G" << std::endl << "{" << std::endl;
		out << "rankdir=LR" << std::endl;
		KMerIndex index(&sequence);
		index.SetupIndex(k);
		KMerSet visit(sequence.Size(), KMerIndex::KMerHashFunction(index.GetK()),
			KMerIndex::KMerEqualTo(index.GetK()));
		for(StrandIterator it = sequence.PositiveBegin(); it.Valid(); ++it)
		{
			ProcessIterator(visit, index, it, out);
		}
		
		for(StrandIterator it = sequence.NegativeBegin(); it.Valid(); ++it)
		{
			ProcessIterator(visit, index, it, out);
		}

		out << "}";
	}

	void GraphAlgorithm::ListNonBranchingPaths(const DNASequence & sequence, size_t k, std::ostream & out, std::ostream & indexOut)
	{
		size_t count = 0;
		KMerIndex index(&sequence);
		index.SetupIndex(k);
		std::vector<std::pair<size_t, StrandIterator> > multiKmer;
		for(StrandIterator it = sequence.PositiveBegin(); it.Valid(); it++)
		{			
			if(it.ProperKMer(k) && (count = index.CountEquivalentKMers(it)) > 1)
			{
				multiKmer.push_back(std::make_pair(count, it));
			}
		}
		
		std::string buf;
		std::vector<StrandIterator> kmer;
		std::vector<Bool> visit(sequence.Size(), false);
		std::sort(multiKmer.begin(), multiKmer.end(), Less);	

		boost::function<StrandIterator& (StrandIterator&)> moveForward = boost::bind(&StrandIterator::operator++, _1);
		boost::function<StrandIterator& (StrandIterator&)> moveBackward = boost::bind(&StrandIterator::operator--, _1);		
		boost::function<bool (const StrandIterator&)> invalid = boost::bind(Invalid, boost::ref(visit), _1);
		boost::function<void (const StrandIterator&)> invalidator = boost::bind(Invalidate, boost::ref(visit), _1);

		for(size_t i = 0; i < multiKmer.size(); i++)
		{
			if(!visit[multiKmer[i].second.GetPosition()])
			{
				index.ListEquivalentKmers(multiKmer[i].second, kmer);
				kmer.erase(std::remove_if(kmer.begin(), kmer.end(), invalid), kmer.end());
				if(kmer.size() > 1)
				{
					size_t forward = Extend(visit, kmer, moveForward, invalidator) + 1;
					size_t backward = Extend(visit, kmer, moveBackward, invalidator);
					if(forward + backward < index.GetK())
					{
						continue;
					}

					std::for_each(kmer.begin(), kmer.end(), invalidator);
					out << "Consensus: " << std::endl;
					StrandIterator end = AdvanceForward(kmer[0], forward);
					StrandIterator start = AdvanceBackward(kmer[0], backward);
					std::copy(start, end, std::ostream_iterator<char>(out));
					out << std::endl;	

					for(size_t j = 0; j < kmer.size(); j++)
					{
						buf.clear();
						end = AdvanceForward(kmer[j], forward);
						start = AdvanceBackward(kmer[j], backward);
						std::pair<size_t, size_t> coord = sequence.SpellOriginal(start, end, std::back_inserter(buf));
						out << (kmer[j].GetDirection() == DNASequence::positive ? '+' : '-') << "s, ";
						out << coord.first << ':' << coord.second << " " << buf << std::endl;
						indexOut << coord.second - coord.first << ' ' << coord.first << ' ' << coord.second << std::endl;
					}

					indexOut << DELIMITER << std::endl;
				}
			}
		}
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
			deletedBulge = 0;
			size_t counter = 0;
			size_t bifurcationCount = EnumerateBifurcations(sequence, k, bifStorage);
			std::cerr << "Total bifurcations: " << bifurcationCount << std::endl;
			std::cerr << "Removing whirls..." << std::endl;
			totalWhirls = RemoveWhirls(bifStorage, sequence, k, minBranchSize);

			std::cerr << "Removing bulges..." << std::endl;
			for(size_t id = 0; id < bifurcationCount; id++)
			{
				totalBulges += RemoveBulges(bifStorage, sequence, k, minBranchSize, id);
			}

			std::cerr << "Total whirls: " << totalWhirls;
			std::cerr << "Total bulges: " << totalBulges << std::endl;		
			std::cerr << "Deleted bpairs by bulge removal: " << deletedBulge << std::endl;
			sequence.Optimize();
		}
		while(totalBulges > 0 || totalWhirls > 0);
	}
	
}
