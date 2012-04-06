#ifndef _DE_BRUIJN_GRAPH_H_
#define _DE_BRUIJN_GRAPH_H_

#include "kmermultiset.h"
#include "syntenyutility.h"


namespace SyntenyBuilder
{
	//This class represents de Bruijn graph. It stores edges (k-mer) in a hash
	//table. Further, more efficient implementation can be provided. 
	class DeBruijnGraph
	{
	public:				
		class Edge;
		class Vertex;
		class ComplementaryEdge;		
		typedef boost::shared_ptr<Edge> EdgePtr;
		typedef boost::shared_ptr<Vertex> VertexPtr;
		typedef StringConstIterator DirectStrandIterator;

		template<class E>
			class GraphElementCreator: public std::unary_function<int, E>
			{
			public:
				GraphElementCreator(const DeBruijnGraph * ptr): ptr_(ptr) {}
				E operator() (int shift_)
				{
					return E(shift_, ptr_);
				}
			private:
				const SyntenyBuilder::DeBruijnGraph * ptr_;
			};

		struct EdgeListPair
		{			
			std::vector<Edge> direct;
			std::vector<ComplementaryEdge> comp;
			EdgeListPair() { }
			EdgeListPair(const std::vector<Edge> & direct, const std::vector<ComplementaryEdge> & comp):
				direct(direct), comp(comp) {}

			size_t Assign(const std::vector<int> & directId, const std::vector<int> & compId, const DeBruijnGraph * graph)
			{
				direct.resize(directId.size());
				comp.resize(compId.size());
				std::transform(directId.begin(), directId.end(), direct.begin(), GraphElementCreator<Edge>(graph));
				std::transform(compId.begin(), compId.end(), comp.begin(), GraphElementCreator<ComplementaryEdge>(graph));
				return Size();
			}

			void TransformIntoPtr(std::vector<Edge*> & edgePtr)
			{
				edgePtr.assign(Size(), 0);
				std::transform(direct.begin(), direct.end(), edgePtr.end(), std::addressof<DeBruijnGraph::Edge>);
				std::transform(comp.begin(), comp.end(), edgePtr.end(), std::addressof<DeBruijnGraph::Edge>);
			}

			template<class F>
				size_t Filter(F f)
				{
					comp.erase(std::remove_if(comp.begin(), comp.end(), f), comp.end());
					direct.erase(std::remove_if(direct.begin(), direct.end(), f), direct.end());
					return Size();
				}			

			void Clear()
			{
				comp.clear();
				direct.clear();				
			}

			size_t Size() const
			{
				return comp.size() + direct.size();
			}			
		};								

		class Vertex
		{
		public:
			bool IsNull() const;			
			bool operator == (const Vertex & vertex) const;
			bool operator != (const Vertex & vertex) const;			
			virtual size_t GetHashCode() const;
			virtual void Spell(std::string::iterator it) const;			
			//ListOutEdges returns list of the edges (k-mers) that leave vertex located at 
			//position "shift". Edges are grouped by equivalency -- if two edges are in the
			//same structure EdgeList, they are equivalent.		
			size_t ListOutEdges(std::vector<EdgePtr> & edgeList) const;
			size_t ListOutEdges(std::vector<EdgeListPair> & edgeList) const;
			Vertex(): shift_(0), graph_(0) {}
		protected:
			friend class DeBruijnGraph;			
			Vertex(int shift, const DeBruijnGraph * graph);
			Vertex(int shift, int lastChar, const DeBruijnGraph * graph): shift_(shift), lastChar_(lastChar), graph_(graph) {}
			int shift_;
			int lastChar_;
			const DeBruijnGraph * graph_;
		};

		class ComplementaryVertex: public Vertex
		{
		public:
			size_t GetHashCode() const;
			ComplementaryVertex() {}
			void Spell(std::string::iterator it) const;
		private:
			friend class DeBruijnGraph;			
			ComplementaryVertex(int shift, const DeBruijnGraph * graph): Vertex(shift, graph) {}
			ComplementaryVertex(int shift, int lastChar, const DeBruijnGraph * graph): Vertex(shift, lastChar, graph) {}
		};

		class Edge
		{
		public:			
			Edge(): graph_(0) {} 
			bool IsNull() const;
			size_t GetHashCode() const;
			virtual VertexPtr EndVertex() const;
			virtual VertexPtr StartVertex() const;		
			virtual EdgePtr NextEdge() const;
			virtual DirectStrandIterator GetPosition() const;
			virtual char Spell() const;
		protected:					
			friend class DeBruijnGraph;
			Edge(int shift, const DeBruijnGraph * graph);
			int startPos_;
			int endPos_;
			const DeBruijnGraph * graph_;
		};

		class ComplementaryEdge: public Edge
		{
		public:
			ComplementaryEdge() {}
			VertexPtr EndVertex() const;
			VertexPtr StartVertex() const;						
			EdgePtr NextEdge() const;
			char Spell() const;
		protected:			
			friend class DeBruijnGraph;
			friend class GraphElementCreator<ComplementaryEdge>;
			void Assign(int shift);
			ComplementaryEdge(int shift, const DeBruijnGraph * graph);
		};

		DeBruijnGraph(const std::string & sequence, int edgeSize);
		size_t GetSize() const;
		size_t GetEdgeSize() const;
		size_t GetVertexSize() const;
		//FindEquivalentEdges returns list of the equivalent edges (k-mers). Reverse
		//complementary variants of the same edge are considered equivalent since we
		//glue reverse complementary vertices.
		size_t FindEquivalentEdges(const Edge & edge, EdgeListPair & edgeList) const;
		//Method for clearing sequence_ from deleted characters.
		void ClearSequence();

		DirectStrandIterator SequenceBegin() const;
		DirectStrandIterator SequenceEnd() const;
		DirectStrandIterator LastVertex() const;
		DirectStrandIterator LastEdge() const;

		VertexPtr ConstructVertex(DirectStrandIterator it) const;
		Edge ConstructEdge(DirectStrandIterator it) const;
	private:
		DISALLOW_COPY_AND_ASSIGN(DeBruijnGraph);				
		static const int END_OF_SEQUENCE = -1;
		//Current version of the sequence (after possible simplification)		
		std::string sequence_;
		//Original version of the sequence (before doing any modifications)
		std::string original_;
		//Length of the k-mer that corresponds to the edge in graph
		const int edgeSize_;
		//Length of the (k - 1)-mer that corresponds to the vertex, equals edgeSize_ - 1
		const int vertexSize_;
		//Number of the vertices and edges in the graph
		int edgeCount_;
		int vertexCount_;		
		//Map from each position in the current sequence to the original sequence
		std::vector<int> position_;
		//Container in which we store all edges
		KMerMultiSet vertexSet_;
		KMerMultiSet edgeMultiSet_;		

		void Init();
		StringIterator SequenceIterator(std::string::iterator it);
		StringConstIterator SequenceIterator(std::string::const_iterator it) const;
		StringReverseIterator SequenceReverseIterator(std::string::iterator it);
		StringConstReverseIterator SequenceReverseIterator(std::string::const_iterator it) const;
		int GetNextPosForward(int shift, int steps = 1) const;
		int GetNextPosBackward(int shift, int steps = 1) const;
	};
}

#endif
