#include "graphalgorithm.h"
//#define _DEBUG
//#undef _DEBUG
namespace SyntenyBuilder
{
	namespace
	{
		typedef char Bool;

		class HashFunction
		{
		public:
			size_t operator () (const DeBruijnGraph::Vertex & v) const
			{
				return v.GetHashCode();
			}
		};		

		struct VisitData
		{
			int classId;
			int instance;
			int distance;

			VisitData() {}
			VisitData(int classId, int instance, int distance):
				classId(classId), instance(instance), distance(distance) {}
		};

		typedef google::dense_hash_set<DeBruijnGraph::Vertex, HashFunction> VertexVisit;
		typedef google::dense_hash_map<DeBruijnGraph::Vertex, std::vector<VisitData>, HashFunction> VertexVisitMap;

		std::ostream& operator << (std::ostream & out, const DeBruijnGraph::Vertex & v)
		{
			v.Spell(std::ostream_iterator<char>(out));
			return out;
		}

		std::ostream& operator << (std::ostream & out, DeBruijnGraph::Edge & e)
		{
			char buf[1 << 8];
			out << e.StartVertex() << " -> " << e.EndVertex() << " ";
			if(e.Direction() == DeBruijnGraph::positive)
			{
				sprintf(buf, "[color=\"%s\", label=\"%i\"];", "blue", e.Position());
			}
			else
			{
				sprintf(buf, "[color=\"%s\", label=\"%i\"];", "red", e.Position());
			}

			return out << buf;
		}

		void OutputProcessIterator(std::ostream & out, VertexVisit & visit, DeBruijnGraph & g, DeBruijnGraph::StrandConstIterator it)
		{
			std::vector<DeBruijnGraph::Edge> edge;
			DeBruijnGraph::Vertex v = g.ConstructVertex(it);
			if(!v.IsNull() && visit.find(v) == visit.end())
			{
				visit.insert(v);				
				g.ListEdges(v, edge);
				for(size_t i = 0; i < edge.size(); i++)
				{
					out << edge[i] << std::endl;
				}
			}
		}

		bool Less(const std::pair<int, DeBruijnGraph::Edge> & p1, const std::pair<int, DeBruijnGraph::Edge> & p2)
		{
			return p1.first < p2.first;
		}

		bool Invalid(const std::vector<Bool> & visit, const DeBruijnGraph::Edge & g)
		{
			return visit[g.Position()] == 1;
		}

		void Invalidate(std::vector<Bool> & visit, const DeBruijnGraph::Edge & g)
		{
			visit[g.Position()] = true;
		}

		void MoveForward(DeBruijnGraph::Edge & edge)
		{
			if(!edge.IsNull())
			{
				edge = edge.NextEdge();
			}
		}

		void MoveBackward(DeBruijnGraph::Edge & edge)
		{
			if(!edge.IsNull())
			{
				edge = edge.PreviousEdge();
			}
		}

		int ExtendForward(std::vector<Bool> & visit, std::vector<DeBruijnGraph::Edge> edge,
			boost::function<void (const DeBruijnGraph::Edge&)> invalidator)
		{
			int ret = 0;
			std::for_each(edge.begin(), edge.end(), MoveForward);			
			for(bool fail = false; !fail; ret++)
			{
				fail = edge[0].IsNull();
				if(!fail)
				{
					char consensus = edge[0].SpellEnd();
					for(size_t i = 1; i < edge.size() && !fail; i++)
					{
						fail = edge[i].IsNull() || visit[edge[i].Position()] || edge[i].SpellEnd() != consensus;
					}
				}

				if(!fail)
				{
					std::for_each(edge.begin(), edge.end(), invalidator);
					std::for_each(edge.begin(), edge.end(), MoveForward);
				}
			}

			return ret - 1;
		}

		int ExtendBackward(std::vector<Bool> & visit, std::vector<DeBruijnGraph::Edge> edge,
			boost::function<void (const DeBruijnGraph::Edge&)> invalidator)
		{
			int ret = 0;
			std::for_each(edge.begin(), edge.end(), MoveBackward);
			for(bool fail = false; !fail; ret++)
			{
				fail = edge[0].IsNull();
				if(!fail)
				{
					char consensus = edge[0].SpellStart();
					for(size_t i = 1; i < edge.size() && !fail; i++)
					{
						fail = edge[i].IsNull() || visit[edge[i].Position()] || edge[i].SpellStart() != consensus;
					}
				}

				if(!fail)
				{
					std::for_each(edge.begin(), edge.end(), invalidator);
					std::for_each(edge.begin(), edge.end(), MoveBackward);
				}
			}

			return ret - 1;
		}				

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

		void FindParallelEdges(DeBruijnGraph & g, DeBruijnGraph::Edge edge, int size, 
			std::vector<DeBruijnGraph::Edge> & ret)
		{
			ret.clear();
			std::vector<DeBruijnGraph::Edge> init;
			g.FindEquivalentEdges(edge, init);
			std::vector<DeBruijnGraph::Edge> now(init);

		#ifdef _DEBUG
			std::cerr << "Path: ";
			PrintPath(edge, size + g.EdgeSize(), std::cerr);
			std::cerr << "Finding parallel paths, candidates: " << std::endl;
			for(size_t i = 0; i < init.size(); i++)
			{
				bool suits = true;
				DeBruijnGraph::Edge test = init[i];
				for(int j = 0; j < size; j++)
				{
					if(test.IsNull())
					{
						suits = false;
						break;
					}

					MoveForward(test);
				}

				if(suits && !test.IsNull())
				{
					PrintPath(init[i], size + g.EdgeSize(), std::cerr);
				}
			}

			std::cerr << "---------------" << std::endl;
		#endif

			for(int step = 0; step < size; step++)
			{
				MoveForward(edge);
				for(size_t i = 0; i < now.size(); i++)
				{	
					if(!now[i].IsNull())
					{
						MoveForward(now[i]);
					}

					if(now[i].IsNull() || now[i].EndVertex() != edge.EndVertex())
					{
						init[i] = DeBruijnGraph::Edge();
					}					
				}								
			}
						
			std::copy(init.begin(), std::remove(init.begin(), init.end(), DeBruijnGraph::Edge()),
				std::back_inserter(ret));
		}

		bool CheckBulge(std::vector<DeBruijnGraph::Edge> edge1, int size1, std::vector<DeBruijnGraph::Edge> edge2, int size2)
		{
			std::vector<int> visit;
			for(int i = 0; i < size1; i++)
			{
				for(size_t j = 0; j < edge1.size(); j++)
				{					
					visit.push_back(edge1[j].Position());
					MoveForward(edge1[j]);
				}
			}

			std::sort(visit.begin(), visit.end());
			for(int i = 0; i < size2; i++)
			{
				for(size_t j = 0; j < edge2.size(); j++)
				{
					int pos = edge2[j].Position();
					if(std::binary_search(visit.begin(), visit.end(), edge2[j].Position()))
					{
						return false;
					}

					MoveForward(edge2[j]);
				}
			}

			return true;
		}

		bool CheckBranch(std::vector<DeBruijnGraph::Edge> edge1, int size1)
		{
			std::vector<int> visit;
			for(int i = 0; i < size1; i++)
			{
				for(size_t j = 0; j < edge1.size(); j++)
				{					
					visit.push_back(edge1[j].Position());
					MoveForward(edge1[j]);
				}
			}

			std::sort(visit.begin(), visit.end());
			for(size_t i = 0; visit.size() > 0 && i < visit.size() - 1; i++)
			{
				if(visit[i] == visit[i + 1])
				{
					return false;
				}
			}

			return true;
		}

		bool CollapseBulge(DeBruijnGraph & g,
			const DeBruijnGraph::Vertex & start,
			const DeBruijnGraph::Vertex & dest, 
			DeBruijnGraph::Edge smallBranch, int smallSize,
			DeBruijnGraph::Edge bigBranch, int bigSize)
		{			
			std::vector<DeBruijnGraph::Edge> parallel;
			std::vector<DeBruijnGraph::Edge> smallParallel;
			FindParallelEdges(g, bigBranch, bigSize, parallel);
			FindParallelEdges(g, smallBranch, smallSize, smallParallel);

			bool bigCheck = CheckBranch(parallel, bigSize);
			bool smallCheck = CheckBranch(smallParallel, smallSize);
			if(!bigCheck && !smallCheck)
			{
				return false;
			}

			if(!bigCheck)
			{
				if(smallSize == bigSize)
				{
					return CollapseBulge(g, start, dest, bigBranch, bigSize, smallBranch, smallSize);
				}

				return false;
			}

			if(!CheckBulge(smallParallel, smallSize, parallel, bigSize))
			{
				return false;
			}

		#ifdef _DEBUG
			static int bulge = 0;
			std::cerr << "Bulge #" << bulge++ << std::endl;
			std::cerr << "Before: " << std::endl;
			PrintRaw(g.sequence, std::cerr);
			std::cerr << "Small branch: " << std::endl;
			for(size_t i = 0; i < smallParallel.size(); i++)
			{
				PrintPath(smallParallel[i], smallSize + g.EdgeSize(), std::cerr);
			}
			
			std::cerr << "Big branch: " << std::endl;
			for(size_t i = 0; i < parallel.size(); i++)
			{
				PrintPath(parallel[i], bigSize + g.EdgeSize(), std::cerr);
			}

		#endif
		
			std::vector<DeBruijnGraph::StrandIterator> it;
			for(size_t i = 0; i < parallel.size(); i++)
			{
				it.push_back(--parallel[i].EndIterator());
			}

			for(int step = 0; step < smallSize; step++)
			{
				for(size_t i = 0; i < it.size(); i++)
				{
					it[i].AssignBase(smallBranch.SpellEnd());
					++it[i];
				}

				MoveForward(smallBranch);
			}

			for(int step = 0; step < bigSize - smallSize; step++)
			{
				for(size_t i = 0; i < it.size(); i++)
				{
					it[i].Invalidate();
				}
			}

		#ifdef _DEBUG
			std::cerr << "After: " << std::endl;
			PrintRaw(g.sequence, std::cerr);
			std::cerr << "Small branch: " << std::endl;
			for(size_t i = 0; i < smallParallel.size(); i++)
			{
				PrintPath(smallParallel[i], smallSize + g.EdgeSize(), std::cerr);
			}
			
			std::cerr << "Big branch: " << std::endl;
			for(size_t i = 0; i < parallel.size(); i++)
			{
				PrintPath(parallel[i], smallSize + g.EdgeSize(), std::cerr);
			}

			std::cerr << "---------------" << std::endl;
		#endif

			return true;
		}		

		bool ProcessBifurcation(DeBruijnGraph & g, DeBruijnGraph::Vertex & v, int minCycleSize)
		{					
			VertexVisitMap visit;
			visit.set_empty_key(DeBruijnGraph::Vertex());			
			std::vector<std::vector<DeBruijnGraph::Edge> > init;
			g.ListEdgesSeparate(v, init);
			std::vector<std::vector<DeBruijnGraph::Edge> > now = init;
			bool keepOn = true;
			for(int step = 1; step <= minCycleSize && keepOn; step++)
			{

			#ifdef _DEBUG
				std::cerr << "step = " << step << std::endl;
			#endif

				keepOn = false;
				for(int classId = 0; classId < static_cast<int>(now.size()); classId++)
				{
					for(int instance = 0; instance < static_cast<int>(now[classId].size()); instance++)
					{
						if(now[classId][instance].IsNull())
						{
							continue;
						}

						keepOn = true;
						DeBruijnGraph::Vertex u = now[classId][instance].EndVertex();

					#ifdef _DEBUG
						u.Spell(std::ostream_iterator<char>(std::cerr));
						std::cerr << std::endl;
					#endif

						VertexVisitMap::iterator it = visit.find(u);
						if(it == visit.end())
						{
							visit.insert(std::make_pair(u, std::vector<VisitData>(1, VisitData(classId, instance, step))));
						}
						else 
						{
							DeBruijnGraph::Edge & nowEdge = init[classId][instance];
							for(size_t i = 0; i < it->second.size(); i++)
							{
								if(it->second[i].classId != classId)
								{
									int prevStep = it->second[i].distance;
									int size = step + prevStep;
									DeBruijnGraph::Edge & prevEdge = init[it->second[i].classId][it->second[i].instance];									
									if(size > 2 && size <= minCycleSize)
									{			
										bool success = false;
										if(step < it->second[i].distance)
										{
											success = CollapseBulge(g, v, u, nowEdge, step - 1, prevEdge, prevStep - 1);
										}
										else
										{
											success = CollapseBulge(g, v, u, prevEdge, prevStep - 1, nowEdge, step - 1);
										}

										if(success)
										{
											return true;
										}
									}
									else
									{
										it->second.push_back(VisitData(classId, instance, step));
									}
								}
							}
						}

						MoveForward(now[classId][instance]);
					}
				}

			#ifdef _DEBUG
				std::cerr << "---------------" << std::endl;
			#endif

			}

			return false;
		}
	}
	
	std::ostream& operator << (std::ostream & out, DeBruijnGraph & g)
	{
		out << "digraph G" << std::endl << "{" << std::endl;
		out << "rankdir=LR" << std::endl;
		VertexVisit visit(g.sequence.Size());
		visit.set_empty_key(DeBruijnGraph::Vertex());
		for(DeBruijnGraph::StrandConstIterator it = g.sequence.PositiveBegin(); it.Valid(); it++)
		{
			OutputProcessIterator(out, visit, g, it);
		}
		
		for(DeBruijnGraph::StrandConstIterator it = g.sequence.NegativeBegin(); it.Valid(); it++)
		{
			OutputProcessIterator(out, visit, g, it);
		}

		return out << "}";
	}
	
	void GraphAlgorithm::ListNonBranchingPaths(DeBruijnGraph & g, std::ostream & out, std::ostream & indexOut)
	{
		int count = 0;
		std::vector<std::pair<int, DeBruijnGraph::Edge> > multiedge;
		for(DeBruijnGraph::StrandIterator it = g.sequence.PositiveBegin(); it.Valid(); it++)
		{
			DeBruijnGraph::Edge edge = g.ConstructPositiveEdge(it);
			if(!edge.IsNull() && (count = g.CountEquivalentEdges(edge)) > 1)
			{
				multiedge.push_back(std::make_pair(count, edge));
			}
		}
		
		std::string buf;
		std::vector<DeBruijnGraph::Edge> edge;
		std::vector<Bool> visit(g.sequence.Size(), false);
		std::sort(multiedge.begin(), multiedge.end(), Less);	
		boost::function<bool (const DeBruijnGraph::Edge&)> isInvalid = boost::bind(Invalid, boost::ref(visit), _1);
		boost::function<void (const DeBruijnGraph::Edge&)> invalidator = boost::bind(Invalidate, boost::ref(visit), _1);
		for(size_t i = 0; i < multiedge.size(); i++)
		{
			if(!visit[multiedge[i].second.Position()])
			{
				g.FindEquivalentEdges(multiedge[i].second, edge);
				edge.erase(std::remove_if(edge.begin(), edge.end(), isInvalid), edge.end());
				if(edge.size() > 1)
				{
					int forward = ExtendForward(visit, edge, invalidator);
					int backward = ExtendBackward(visit, edge, invalidator);
					std::for_each(edge.begin(), edge.end(), invalidator);
					out << "Consensus: " << std::endl;
					DeBruijnGraph::StrandConstIterator end = Advance(edge[0].EndIterator(), forward);
					DeBruijnGraph::StrandConstIterator start = AdvanceBackward(edge[0].StartIterator(), backward);
					std::copy(start, end, std::ostream_iterator<char>(out));
					out << std::endl;					
					for(size_t j = 0; j < edge.size(); j++)
					{
						buf.clear();
						end = Advance(edge[j].EndIterator(), forward);
						start = AdvanceBackward(edge[j].StartIterator(), backward);
						std::pair<int, int> coord = g.sequence.SpellOriginal(start, end, std::back_inserter(buf));
						out << (edge[j].Direction() == DeBruijnGraph::positive ? '+' : '-') << "s, ";
						out << coord.first << ':' << coord.second << " " << buf << std::endl;
						indexOut << coord.second - coord.first << ' ' << coord.first << ' ' << coord.second << std::endl;
					}

					indexOut << std::string(80, '-') << std::endl;
				}
			}
		}
	}	

	int GraphAlgorithm::Simplify(DeBruijnGraph & g, int minCycleSize)
	{
	#ifdef _DEBUG
		std::cerr << std::endl;
	#endif
		int bulgeCount = 0;
		int bifurcationCount = 0;
 		g.sequence.KeepHash(g.VertexSize());
		std::vector<std::vector<DeBruijnGraph::Edge> > edge;
		for(DeBruijnGraph::StrandConstIterator it = g.sequence.PositiveBegin(); it.Valid(); it++)
		{
			if(it.GetPosition() % 1000000 == 0)
			{
				std::cerr << it.GetPosition() << " " << bifurcationCount << " " << bulgeCount << std::endl;
			}

			int pos = it.GetPosition();

			DeBruijnGraph::Vertex v = g.ConstructVertex(it);
			if(!v.IsNull() && g.ListEdgesSeparate(v, edge) > 1)
			{
				bifurcationCount++;
				while(ProcessBifurcation(g, v, minCycleSize))
				{
					bulgeCount++;
				}
			}
		}

		g.sequence.Optimize();
		return bulgeCount;
	}
}
