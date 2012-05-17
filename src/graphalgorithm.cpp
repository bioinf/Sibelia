#include "graphalgorithm.h"
//#define _DEBUG
//#undef _DEBUG
namespace SyntenyBuilder
{
	namespace
	{
		int deletedWhirl;
		int deletedBulge;
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
			char classId;
			DeBruijnGraph::Edge * edge;
			int distance;
			VisitData() {}
			VisitData(char classId, DeBruijnGraph::Edge * edge, int distance):
				classId(classId), edge(edge), distance(distance) {}

			bool operator < (const VisitData & toCompare)
			{
				return std::make_pair(classId, distance) < std::make_pair(toCompare.classId, toCompare.distance);
			}
		};

		const int ANY_CLASS = -1;		
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
		
		void CollapseBulge(DeBruijnGraph & g, const VisitData & source, const VisitData & target)
		{
			
		#ifdef _DEBUG
			static int bulge = 0;
			std::cerr << "Bulge #" << bulge++ << std::endl;
			std::cerr << "Before: " << std::endl;
			PrintRaw(g.sequence, std::cerr);
			std::cerr << "Source branch: " << std::endl;			
			int shift = source.distance == 0 ? 1 : 0;
			PrintPath(*source.edge, source.distance + g.EdgeSize() - shift, std::cerr);
			std::cerr << "Target branch: " << std::endl;			
			PrintPath(*target.edge, target.distance + g.EdgeSize(), std::cerr);
			std::cerr << std::endl;
		#endif
			
			DNASequence::StrandIterator it = --source.edge->EndIterator();
			DNASequence::StrandIterator jt = --target.edge->EndIterator();
			for(int i = 0; i < source.distance; i++, ++it, ++jt)
			{
				it.AssignBase(*jt);
			}

			int diff = target.distance - source.distance;
			for(int i = 0; i < diff; i++)
			{
				deletedBulge++;
				jt.Invalidate();				
			}

		#ifdef _DEBUG
			std::cerr << std::endl << "After: " << std::endl;
			PrintRaw(g.sequence, std::cerr);
			std::cerr << "Source branch: " << std::endl;			
			PrintPath(*source.edge, source.distance + g.EdgeSize() - shift, std::cerr);
			std::cerr << "Target branch: " << std::endl;			
			PrintPath(*target.edge, source.distance + g.EdgeSize() - shift, std::cerr);
			std::cerr << std::endl << std::string(80, '-') << std::endl;
		#endif
		}

		int ProcessVertex(DeBruijnGraph & g, VertexVisitMap::iterator it, 
			std::map<DeBruijnGraph::Edge*, int> & deprecate)
		{
			int ret = 0;
			int source = -1;				
			for(int i = 0; i < static_cast<int>(it->second.size()); i++)
			{
				if(deprecate[it->second[i].edge] <= it->second[i].distance)
				{
					source = i;
					break;
				}
			}

			if(source >= 0)
			{
				char classId = it->second[source].classId;
				for(int i = source + 1; i < static_cast<int>(it->second.size()); i++)
				{
					if(deprecate[it->second[i].edge] <= it->second[i].distance &&
						it->second[i].classId != classId)
					{
						ret++;
						CollapseBulge(g, it->second[source], it->second[i]);
						deprecate[it->second[i].edge] = it->second[i].distance;
					}
				}
			}

			return ret;
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
			std::cerr << std::endl;
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
			std::set<int> deprecate;
			std::vector<DeBruijnGraph::Edge> now;
			std::vector<DeBruijnGraph::Edge> start;			
			int iter = 0;
			while(keepOn)
			{
				keepOn = false;
				deprecate.clear();
				g.ListEdges(v, start);
				now = start;
				for(int step = 0; step < minCycleSize; step++)
				{
					for(size_t i = 0; i < static_cast<int>(now.size()); i++)
					{
						if(!now[i].IsNull() && deprecate.find(now[i].Position()) == deprecate.end())
						{
							DeBruijnGraph::Vertex u = now[i].EndVertex();
							if(u == v)
							{
								ret++;
								keepOn = true;
								MoveForward(now[i]);
								if(!now[i].IsNull())
								{
									deprecate.insert(now[i].Position());
								}

								CollapseWhirl(g, start[i], step);
							}						
							else
							{
								MoveForward(now[i]);
							}
						}
					}
				}
			}

			return ret;
		}

		int RemoveBulges(DeBruijnGraph & g, DeBruijnGraph::Vertex & v, int minCycleSize)
		{			
			int ret = 0;
			VertexVisitMap visit;
			visit.set_empty_key(DeBruijnGraph::Vertex());
			std::map<DeBruijnGraph::Edge*, int> deprecate;

			std::vector<std::vector<DeBruijnGraph::Edge> > now;
			std::vector<std::vector<DeBruijnGraph::Edge> > start;
			g.ListEdgesSeparate(v, start);
			now = start;
			for(int classId = 0; classId < static_cast<int>(now.size()); classId++)
			{
				for(int instance = 0; instance < static_cast<int>(now[classId].size()); instance++)
				{
					deprecate[&now[classId][instance]] = -1;
				}
			}

			for(int step = 0; step < minCycleSize; step++)
			{
				for(int classId = 0; classId < static_cast<int>(now.size()); classId++)
				{
					for(int instance = 0; instance < static_cast<int>(now[classId].size()); instance++)
					{
						if(!now[classId][instance].IsNull())
						{
							DeBruijnGraph::Vertex u = now[classId][instance].EndVertex();
							visit[u].push_back(VisitData(classId, &start[classId][instance], step));
							MoveForward(now[classId][instance]);
						}
					}
				}
			}

			for(VertexVisitMap::iterator it = visit.begin(); it != visit.end(); it++)
			{
				ret += ProcessVertex(g, it, deprecate);
			}

			return ret;
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
		int whirlCount = 0;
		int bifurcationCount = 0;
		deletedWhirl = deletedBulge = 0;
 		g.sequence.KeepHash(g.VertexSize());
		std::vector<std::vector<DeBruijnGraph::Edge> > edge;
		for(DeBruijnGraph::StrandConstIterator it = g.sequence.PositiveBegin(); it.Valid(); it++)
		{
			if(it.GetPosition() % 10000 == 0)
			{
				std::cerr << it.GetPosition() << " " << bifurcationCount << " " << whirlCount << ' ' << bulgeCount << std::endl;
			}

			DeBruijnGraph::Vertex v = g.ConstructVertex(it);
			if(!v.IsNull() && g.ListEdgesSeparate(v, edge) > 1)
			{
				bifurcationCount++;
				whirlCount += RemoveWhirls(g, v, minCycleSize);
				bulgeCount += RemoveBulges(g, v, minCycleSize);
			}
		}

		std::cerr << "Total bifurcations: " << bifurcationCount << std::endl;
		std::cerr << "Deleted bpairs by whirl removal: " << deletedWhirl << std::endl;
		std::cerr << "Deleted bpairs by bulge removal: " << deletedBulge << std::endl;
		std::cerr << "Total whirls: " << whirlCount << std::endl;
		std::cerr << "Total bulges: " << bulgeCount << std::endl;		
		g.sequence.Optimize();
		return bulgeCount;
	}
}
