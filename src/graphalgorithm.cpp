#include "graphalgorithm.h"

namespace SyntenyBuilder
{
	namespace
	{
		template<class E>
			class HashFunctionPtr
			{
			public:
				size_t operator() (E e) const
				{
					return e->GetHashCode();
				}
			};

		template<class E>
			class EqualPtr
			{
			public:
				bool operator() (E e1, E e2) const
				{
					return *e1 == *e2;
				}
			};

		class EdgeInvalid
		{
		public:
			EdgeInvalid(const std::vector<Bool> & invalid): invalid_(&invalid) {}
			bool operator () (const DeBruijnGraph::Edge & e)
			{
				return (*invalid_)[e.GetHashCode()] == 0;
			}
		private:
			const std::vector<Bool> * invalid_;
		};
		
		class InvalidateEdge
		{
		public:
			InvalidateEdge(std::vector<char> & invalid): invalid_(&invalid) {}
			void operator () (const DeBruijnGraph::Edge & e)
			{
				(*invalid_)[e.GetHashCode()] = true;
			}
		private:
			std::vector<Bool> * invalid_;
		};

		bool Less(const std::pair<int, DeBruijnGraph::Edge> & p1, const std::pair<int, DeBruijnGraph::Edge> & p2)
		{
			return p1.first < p2.first;
		}

		template<class F, class G>
			void ExtendForward(DeBruijnGraph::EdgeListPair edge, F checker, G invalidator)
			{
				int ret = 0;
				std::vector<DeBruijnGraph::Edge*> edgePtr;
				edge.TransformIntoPtr(edgePtr);
				return ret;

			}
		template<class F, class G>
			int ExtendBackward(DeBruijnGraph::EdgeListPair edge, F checker, G invalidator)
			{
				int ret = 0;
				
				return ret;
			}
	}


	void GraphAlgorithm::ListNonBranchingPaths(const DeBruijnGraph & g)
	{/*
		//Find all direct multiedges in the graph (edges with multiplicity > 1)
		DeBruijnGraph::EdgeListPair ret;
		std::vector<std::pair<int, DeBruijnGraph::Edge> > multiedge;
		DeBruijnGraph::DirectStrandIterator end = g.LastEdge();
		for(DeBruijnGraph::DirectStrandIterator it = g.SequenceBegin(); it != end; it++)
		{
			DeBruijnGraph::Edge edge = g.ConstructEdge(it);
			if(g.FindEquivalentEdges(edge, ret) > 1 && ret.direct.size() > 0)
			{
				multiedge.push_back(std::make_pair(ret.Size(), edge));
			}
		}

		std::sort(multiedge.begin(), multiedge.end(), Less);	
		std::vector<Bool> visit(g.GetSize(), false);
		EdgeInvalid checker(visit);
		InvalidateEdge invalidator(visit);
		std::vector<DeBruijnGraph::DirectStrandIterator> pathStart;
		std::vector<DeBruijnGraph::DirectStrandIterator> pathEnd;
		for(size_t i = 0; i < multiedge.size(); i++)
		{
			g.FindEquivalentEdges(multiedge[i].second, ret);
			if(ret.Filter(checker) > 1)
			{
		//		int forward = ExtendForward(ret, checker, invalidator);
		//		int backward = ExtendBackward(ret, checker, invalidator);
			}
		}*/
	}

	void GraphAlgorithm::DebugOutput(const DeBruijnGraph & g, std::ostream & output)
	{
		output << "Direct strand: ";
		std::copy(g.SequenceBegin(), g.SequenceEnd(), std::ostream_iterator<char>(output));

		output << std::endl;
		std::string buf(g.GetVertexSize(), ' ');
		std::vector<DeBruijnGraph::EdgePtr> edge;
		google::sparse_hash_set<DeBruijnGraph::VertexPtr,
								HashFunctionPtr<DeBruijnGraph::VertexPtr>,
								EqualPtr<DeBruijnGraph::VertexPtr> > visit;

		DeBruijnGraph::DirectStrandIterator end = g.LastVertex();
		for(DeBruijnGraph::DirectStrandIterator it = g.SequenceBegin(); it != end; it++)
		{
			DeBruijnGraph::VertexPtr v = g.ConstructVertex(it);
			if(visit.find(v) == visit.end())
			{
				visit.insert(v);
				v->Spell(buf.begin());
				output << buf << std::endl;
				v->ListOutEdges(edge);
				for(size_t i = 0; i < edge.size(); i++)
				{
					edge[i]->EndVertex()->Spell(buf.begin());
					std::cout << '\t' << buf << std::endl;
				}
			}
		}
	}
}

