#include "graphalgorithm.h"
//#define _DEBUG
//#undef _DEBUG
namespace SyntenyBuilder
{
	namespace
	{
		typedef char Bool;		
		typedef DNASequence::StrandIterator StrandIterator;

	/*	int pos;
		int deletedWhirl;
		int deletedBulge;
		const int ANY_CLASS = -1;					
		

		struct VertexHashFunction
		{
			VertexHashFunction() {}
			size_t operator () (const DeBruijnGraph::Vertex & v) const
			{
				return v.GetHashCode();
			}
		};		

		struct VisitData
		{
			char classId;
			DeBruijnGraph::Edge * edge;
			int distance;
			VisitData() {}
			VisitData(char classId, DeBruijnGraph::Edge * edge, int distance):
				classId(classId), edge(edge), distance(distance) {}

		};

		struct VisitDataComparer
		{
		public:
			bool operator () (const VisitData & data)
			{
				return data.edge == e;
			}

			VisitDataComparer(DeBruijnGraph::Edge * e): e(e) {}

		private:
			DeBruijnGraph::Edge * e;
		};
		
		struct VertexEqual
		{
			VertexEqual() {}
			bool operator () (const DeBruijnGraph::Vertex & a, const DeBruijnGraph::Vertex & b) const
			{
				return a.Equal(b);
			}
		};*/	

		typedef boost::unordered_set<KMerIndex::StrandIterator, KMerIndex::KMerHashFunction,
			KMerIndex::KMerEqualTo> KMerVisit;

		/*

		typedef boost::unordered_multimap<DeBruijnGraph::Vertex, VisitData,
			VertexHashFunction, VertexEqual> VertexVisitMap; 
		*/

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
		
		void ProcessIterator(KMerVisit & visit, const KMerIndex & index, DNASequence::StrandIterator it, std::ostream & out)
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
		/*
		std::string ToString(const DeBruijnGraph::Vertex & v)
		{
			std::string buf;
			v.Spell(std::back_inserter(buf));
			return buf;
		}

		void PrintRaw(DNASequence & s, std::ostream & out)
		{
			std::string rcomp;
			s.SpellRaw(std::ostream_iterator<char>(out));
			out << std::endl;
			for(int i = 0; i < s.Size(); i++)
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

		void PrintPath(const DeBruijnGraph::Edge & e, int size, std::ostream & out)
		{
			out << e.StartIterator().GetPosition() << " ";
			out << (e.Direction() == DeBruijnGraph::positive ? "s+ " : "s- ");
			CopyN(e.StartIterator(), size, std::ostream_iterator<char>(out));
			std::cerr << std::endl;
		}			

		void CollapseWhirl(DeBruijnGraph & g, DeBruijnGraph::Edge & e, int distance)
		{
		#ifdef _DEBUG
			static int whirl = 0;
			std::cerr << "Whirl #" << whirl++ << std::endl;
			std::cerr << "Before: " << std::endl;
			PrintRaw(g.sequence, std::cerr);
			std::cerr << "Whirl branch: " << std::endl;			
			PrintPath(e, distance + g.EdgeSize(), std::cerr);	
		#endif

			DNASequence::StrandIterator it = --e.EndIterator();
			for(int i = 0; i <= distance; i++)
			{
				deletedWhirl++;
				it.Invalidate();
			}

		#ifdef _DEBUG
			std::cerr << std::endl << "After: " << std::endl;
			PrintRaw(g.sequence, std::cerr);
			std::cerr << std::endl << std::string(80, '-') << std::endl;
		#endif
		}

		int RemoveWhirls(DeBruijnGraph & g, DeBruijnGraph::Vertex & v, int minCycleSize)
		{
			int ret = 0;
			bool keepOn = true;
			std::vector<DeBruijnGraph::Edge> now;
			std::vector<DeBruijnGraph::Edge> start;			
			while(keepOn)
			{
				keepOn = false;
				g.ListEdges(v, start);
				now = start;
				for(int step = 0; step < minCycleSize && !keepOn; step++)
				{
					for(size_t i = 0; i < static_cast<int>(now.size()) && !keepOn; i++)
					{
						if(!now[i].IsNull())
						{
							DeBruijnGraph::Vertex u = now[i].EndVertex();
							if(u.Equal(v))
							{
								ret++;
								keepOn = true;
								CollapseWhirl(g, start[i], step);
							}						
							else
							{
								now[i] = now[i].NextEdge();
							}
						}
					}
				}
			}

			return ret;
		}		

		void ClearVisit(VertexVisitMap & visit, VisitData target)
		{
			DeBruijnGraph::Edge e(*target.edge);
			for(int i = 0; i <= target.distance; i++)
			{
				std::pair<VertexVisitMap::iterator, VertexVisitMap::iterator>
					range = visit.equal_range(e.EndVertex());
				for(VertexVisitMap::iterator it = range.first; it != range.second; ++it)
				{
					if(it->second.edge == target.edge)
					{
						visit.erase(it);
						break;
					}
				}

				e = e.NextEdge();
			}
		}

		DeBruijnGraph::Edge CollapseBulge(DeBruijnGraph & g, VertexVisitMap & visit,
			VisitData source, VisitData target)
		{
			
		#ifdef _DEBUG
			static int bulge = 0;
			std::cerr << "Bulge #" << bulge++ << std::endl;
			std::cerr << "Before: " << std::endl;
			PrintRaw(g.sequence, std::cerr);
			std::cerr << "Source branch: " << std::endl;			
			PrintPath(*source.edge, source.distance + g.EdgeSize(), std::cerr);
			std::cerr << "Target branch: " << std::endl;			
			PrintPath(*target.edge, target.distance + g.EdgeSize(), std::cerr);
		#endif
			
			ClearVisit(visit, target);
			DNASequence::StrandIterator it = --source.edge->EndIterator();
			DNASequence::StrandIterator jt = --target.edge->EndIterator();
			DeBruijnGraph::Edge ret = *source.edge;
			for(int i = 0; i <= source.distance; i++, ++it, ++jt)
			{
				jt.AssignBase(*it);
				ret = ret.NextEdge();
			}

			int diff = target.distance - source.distance;
			for(int i = 0; i < diff; i++)
			{
				deletedBulge++;
				jt.Invalidate();				
			}

		#ifdef _DEBUG
			std::cerr << "After: " << std::endl;
			PrintRaw(g.sequence, std::cerr);
			std::cerr << "Source branch: " << std::endl;			
			PrintPath(*source.edge, source.distance + g.EdgeSize(), std::cerr);
			std::cerr << "Target branch: " << std::endl;			
			PrintPath(*target.edge, source.distance + g.EdgeSize(), std::cerr);
			std::cerr << std::string(80, '-') << std::endl;
		#endif

			return ret;
		}	

		int RemoveBulges(DeBruijnGraph & g, DeBruijnGraph::Vertex & v, int minCycleSize)
		{			
			int ret = 0;
			VertexVisitMap visit;
			std::set<int> passed;
			std::vector<std::vector<DeBruijnGraph::Edge> > now;
			std::vector<std::vector<DeBruijnGraph::Edge> > start;
			g.ListEdgesSeparate(v, start);
			std::vector<std::vector<int> > travelRange(start.size());

			//This is a workaround for some conceptual flaw
			for(size_t i = 0; i < start.size(); i++)
			{
				travelRange[i].assign(start[i].size(), 0);
				for(size_t j = 0; j < start[i].size(); j++)
				{
					if(start[i][j].Direction() == DeBruijnGraph::negative)
					{
						DeBruijnGraph::StrandConstIterator it = start[i][j].StartIterator();
						for(int pos = 0; pos < g.EdgeSize() - 1; pos++, ++it)
						{
							passed.insert(it.GetPosition());
						}
					}
				}
			}
			//Workaround ends

			now = start;
			int nowPos = pos;
			for(int step = 0; step < minCycleSize; step++)
			{
				for(int classId = 0; classId < static_cast<int>(now.size()); classId++)
				{
					for(int instance = 0; instance < static_cast<int>(now[classId].size()); instance++)
					{
						if(!now[classId][instance].IsNull() && passed.count(now[classId][instance].Position()) == 0)
						{
							passed.insert(now[classId][instance].Position());
							DeBruijnGraph::Vertex u = now[classId][instance].EndVertex();					
							if(u.Equal(v) || u.IsNull())
							{
								continue;
							}

							bool collapsed = false;
							std::pair<VertexVisitMap::iterator, VertexVisitMap::iterator> range = visit.equal_range(u);						
							VisitData nowData = VisitData(classId, &start[classId][instance], travelRange[classId][instance]);
							for(VertexVisitMap::iterator it = range.first; it != range.second; ++it)
							{								
								if(classId != it->second.classId)								
								{
									ret++;			
									collapsed = true;
									now[classId][instance] = CollapseBulge(g, visit, it->second, nowData);
									travelRange[classId][instance] = it->second.distance;
									break;
								}
							}
							
							if(!collapsed)
							{
								travelRange[classId][instance]++;
								visit.insert(std::make_pair(u, nowData));
								now[classId][instance] = now[classId][instance].NextEdge();
							}
						}
					}
				}
			}
	
			return ret;
		}*/
	}
	

	void GraphAlgorithm::SerializeGraph(KMerIndex & index, const DNASequence & sequence, std::ostream & out)
	{
		out << "digraph G" << std::endl << "{" << std::endl;
		out << "rankdir=LR" << std::endl;
		KMerVisit visit(sequence.Size(), KMerIndex::KMerHashFunction(index.GetK()),
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
		std::vector<char> visit(sequence.Size(), false);
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

					indexOut << std::string(80, '-') << std::endl;
				}
			}
		}
	}	

	/*
	int GraphAlgorithm::Simplify(DeBruijnGraph & g, int minCycleSize)
	{
	#ifdef _DEBUG
		std::cerr << std::endl;
	#endif
		int bulgeCount = 0;
		int whirlCount = 0;
		int bifurcationCount = 0;
		deletedWhirl = deletedBulge = 0;
		g.sequence.KeepHash(g.VertexSize());
		google::sparse_hash_set<std::string> globalVisit;
		std::vector<std::vector<DeBruijnGraph::Edge> > edge;
		const int MOD = 100000;

		std::cerr << std::string(50, '-') << std::endl;

		int counter = 0;
		std::string buf(g.VertexSize(), ' ');
		for(DeBruijnGraph::StrandConstIterator it = g.sequence.PositiveBegin(); it.Valid(); it++)
		{
			if(++counter % MOD == 0)
			{
				std::cerr << it.GetPosition() << ' ' << bifurcationCount << ' ';
				std::cerr << bulgeCount << ' ' <<  deletedBulge << ' ' << std::endl;
			}
			
 			DeBruijnGraph::Vertex v = g.ConstructVertex(it);
			if(!v.IsNull() && g.ListEdgesSeparate(v, edge) > 1)
			{
				v.Spell(buf.begin());
				if(globalVisit.count(buf) == 0)
				{
					bifurcationCount++;
					globalVisit.insert(buf);
					bulgeCount += RemoveBulges(g, v, minCycleSize);
				}
			}
		}

		std::cerr << "Total bifurcations: " << bifurcationCount << std::endl;
		std::cerr << "Deleted bpairs by bulge removal: " << deletedBulge << std::endl;
		std::cerr << "Total bulges: " << bulgeCount << std::endl;		
		g.sequence.Optimize();
		g.sequence.DropHash();
		return bulgeCount;
	}
	*/
}
