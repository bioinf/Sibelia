#include "graphalgorithm.h"
//#define _DEBUG
//#undef _DEBUG
namespace SyntenyBuilder
{
	namespace
	{
		typedef char Bool;		
		typedef DNASequence::StrandIterator StrandIterator;
		typedef google::sparse_hash_set<StrandIterator, KMerIndex::KMerHashFunction,
			KMerIndex::KMerEqualTo> KMerSet;	
		typedef std::pair<StrandIterator, StrandIterator> Vertex;
		typedef std::vector<Vertex> VertexVector;

		struct VisitData
		{
			size_t kmerId;
			size_t distance;
			VisitData() {}
			VisitData(size_t kmerId, size_t distance): kmerId(kmerId), distance(distance) {}

		};

		typedef boost::unordered_multimap<StrandIterator, VisitData,
			KMerIndex::KMerHashFunction, KMerIndex::KMerEqualTo> VertexVisitMap;

		size_t totalBulges;
		size_t deletedBulge;
		size_t totalBifurcation;

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

		void ClearVisit(VertexVisitMap & visit, const std::vector<StrandIterator> & startVertex,
			VisitData targetData)
		{			
			StrandIterator target = startVertex[targetData.kmerId];
			for(size_t i = 0; i < targetData.distance; i++)
			{
				std::pair<VertexVisitMap::iterator, VertexVisitMap::iterator>
					range = visit.equal_range(++target);
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

		void CollapseBulge(KMerIndex & index,
			DNASequence & sequence,
			VertexVisitMap & visit,
			const std::vector<StrandIterator> & startVertex,
			VisitData sourceData,
			VisitData targetData)
		{
			size_t k = index.GetK();

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

			ClearVisit(visit, startVertex, targetData);
			StrandIterator sourceIt = startVertex[sourceData.kmerId];
			StrandIterator targetIt = startVertex[targetData.kmerId];
			size_t diff = targetData.distance - sourceData.distance;
			deletedBulge += diff;
			sequence.CopyN(sourceIt, sourceData.distance, targetIt);
			targetIt.Jump(sourceData.distance);
			sequence.EraseN(targetIt, diff);
			
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

		void RemoveBulges(KMerIndex & index, DNASequence & sequence, StrandIterator vertex, size_t minBranchSize)
		{			
			size_t k = index.GetK();
			std::set<size_t> passed;
			VertexVisitMap visit(0, KMerIndex::KMerHashFunction(k), KMerIndex::KMerEqualTo(k));			
			std::vector<StrandIterator> startVertex;
			VertexVector nowVertex;
			index.ListEquivalentKmers(vertex, startVertex);			
			if(startVertex.size() < 2)
			{
				return;
			}

			totalBifurcation++;
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
			for(std::vector<StrandIterator>::iterator it = startVertex.begin(); it != startVertex.end(); it++)
			{
				for(size_t j = 0; j < k; j++)
				{
					passed.insert(it->GetPosition());
					++(*it);
				}

				*it = AdvanceBackward(*it, k);
			}

			KMerIndex::KMerEqualTo equal(k);
			for(size_t step = 0; step < minBranchSize; step++)
			{
				for(size_t kmerId = 0; kmerId < nowVertex.size(); kmerId++)
				{
					StrandIterator & kmerStart = nowVertex[kmerId].first;
					StrandIterator & kmerEnd = nowVertex[kmerId].second;
					if(kmerEnd.Valid() && passed.count(kmerEnd.GetPosition()) == 0)
					{
						if(equal(kmerStart, vertex))
						{
							continue;
						}

						bool collapsed = false;
						passed.insert(kmerEnd.GetPosition());
						VisitData nowData = VisitData(kmerId, travelRange[kmerId]);

						std::pair<VertexVisitMap::iterator, VertexVisitMap::iterator> range = visit.equal_range(kmerStart);						
						for(VertexVisitMap::iterator it = range.first; it != range.second; ++it)
						{								
							if(endChar[kmerId] != endChar[it->second.kmerId])
							{
								totalBulges++;
								collapsed = true;
								CollapseBulge(index, sequence, visit, startVertex, it->second, nowData);
								travelRange[kmerId] = it->second.distance;
								endChar[kmerId] = endChar[it->second.kmerId];
								break;
							}
						}
							
						if(!collapsed)
						{
							travelRange[kmerId]++;
							visit.insert(std::make_pair(kmerStart, nowData));
							++kmerStart;
							++kmerEnd;
						}
					}
				}
			}
		}
	}
	
	
	void GraphAlgorithm::SerializeGraph(const KMerIndex & index, const DNASequence & sequence, std::ostream & out)
	{
		out << "digraph G" << std::endl << "{" << std::endl;
		out << "rankdir=LR" << std::endl;
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

	void GraphAlgorithm::ListNonBranchingPaths(const KMerIndex & index, const DNASequence & sequence, std::ostream & out, std::ostream & indexOut)
	{
		size_t count = 0;
		size_t k = index.GetK();
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
	
	void GraphAlgorithm::SimplifyGraph(KMerIndex & index, DNASequence & sequence, size_t minBranchSize)
	{
		totalBulges = 0;
		deletedBulge = 0;
		size_t counter = 0;
		const size_t MOD = 100000;

		std::cerr << DELIMITER << std::endl;
		std::cerr << "Finding all bifurcations in the graph..." << std::endl;		

		size_t k = index.GetK();
		KMerSet bifurcation(0, KMerIndex::KMerHashFunction(k), KMerIndex::KMerEqualTo(k));
		for(StrandIterator it = sequence.PositiveBegin(); it != sequence.PositiveRightEnd(); ++it)
		{
			if(it.ProperKMer(k + 1) && index.CountEquivalentKMers(it) > 1)
			{
				bifurcation.insert(it);
			}
		}

		for(KMerSet::iterator it = bifurcation.begin(); it != bifurcation.end(); ++it)
		{
			if(it->Valid() && it->ProperKMer(k + 1))
			{
				RemoveBulges(index, sequence, *it, minBranchSize);
			}
		}

		std::cerr << "Total bifurcations: " << totalBifurcation << std::endl;
		std::cerr << "Deleted bpairs by bulge removal: " << deletedBulge << std::endl;
		std::cerr << "Total bulges: " << totalBulges << std::endl;		
		sequence.DropHash();
	}
	
}
