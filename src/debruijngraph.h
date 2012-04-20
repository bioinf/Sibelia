#ifndef _DE_BRUIJN_GRAPH_H_
#define _DE_BRUIJN_GRAPH_H_

#include "dnasequence.h"
#include "indexmultiset.h"

namespace SyntenyBuilder
{
	//This class represents de Bruijn graph. It stores edges (k-mer) in a hash table
	class DeBruijnGraph
	{
	public:				
		class Edge;
		class Vertex;
		typedef DNASequence::StrandIterator StrandIterator;
		typedef DNASequence::StrandConstIterator StrandConstIterator;
		
		enum DirectionTag
		{
			positive,
			negative
		};

		class Vertex
		{
		public:
			bool IsNull() const
			{
				return graph_ == 0;
			}

			size_t GetHashCode() const
			{
				KMerHashFunction f(graph_->vertexSize_);
				return f(it_);
			}

			bool operator == (const Vertex & toCompare) const
			{
				if(graph_ != toCompare.graph_)
				{
					return false;
				}

				if(IsNull() || toCompare.IsNull())
				{
					return IsNull() && toCompare.IsNull();
				}

				if(toCompare.GetHashCode() != GetHashCode())
				{
					return false;
				}

				StrandConstIterator end = Advance(it_, graph_->vertexSize_);
				return std::mismatch(it_, end, toCompare.it_).first == end;
			}

			bool operator != (const Vertex & toCompare) const
			{
				return !(*this == toCompare);
			}

			template<class OutputIterator>
				void Spell(OutputIterator out) const
				{
					std::copy_n(it_, graph_->vertexSize_, out);
				}

			Vertex(): graph_(0) {}
		private:			
			Vertex(const DeBruijnGraph * graph, StrandConstIterator it): graph_(graph), it_(it) {}			
			const DeBruijnGraph * graph_;
			StrandConstIterator it_;
			friend class DeBruijnGraph;
		};
		
		class Edge
		{
		public:
			bool IsNull() const
			{
				return graph_ == 0;
			}

			char SpellEnd() const
			{
				return *end_;
			}

			char SpellStart() const
			{
				return *start_;
			}

			DirectionTag Direction() const
			{
				return tag_;
			}

			Vertex StartVertex() const
			{
				return graph_->ConstructVertex(start_);
			}

			Vertex EndVertex() const
			{
				return graph_->ConstructVertex(++StrandConstIterator(start_));
			}
			
			Edge NextEdge() const
			{
				if(!(++StrandIterator(end_)).Valid())
				{
					return Edge();
				}

				return Edge(graph_, ++StrandIterator(start_), ++StrandIterator(end_), tag_);
			}

			Edge PreviousEdge() const
			{
				StrandIterator prev = --StrandIterator(start_);
				if(!prev.Valid())
				{
					return Edge();
				}

				return Edge(graph_, prev, --StrandIterator(end_), tag_);
			}

			int Position() const
			{
				if(tag_ == positive)
				{
					return start_.GetPosition();
				}

				return end_.GetPosition();
			}
			
			StrandIterator StartIterator() const
			{
				return start_;
			}

			StrandIterator EndIterator() const
			{
				return ++StrandIterator(end_);
			}

			Edge(): graph_(0) {}
		private:
			Edge(DeBruijnGraph * graph, StrandIterator start, DirectionTag tag):
				graph_(graph),
				start_(start),
				end_(Advance(start_, graph_->edgeSize_ - 1)),
				tag_(tag) {}
			Edge(DeBruijnGraph * graph, StrandIterator start, StrandIterator end, DirectionTag tag):
				graph_(graph),
				start_(start),
				end_(end),
				tag_(tag) {}
			DeBruijnGraph * graph_;
			StrandIterator start_;
			StrandIterator end_;
			DirectionTag tag_;
			friend class DeBruijnGraph;
		};
				
		int VertexSize() const
		{
			return vertexSize_;
		}

		int EdgeSize() const
		{
			return edgeSize_;
		}

		Vertex ConstructVertex(StrandConstIterator it) const
		{
			if((it.GetDirection() == DNASequence::positive && it.GetPosition() > positiveVertexBound_) ||
				(it.GetDirection() == DNASequence::negative && it.GetPosition() < negativeVertexBound_))
			{
				return Vertex();
			}

			return Vertex(this, it);
		}

		Edge ConstructPositiveEdge(StrandConstIterator it)
		{
			if(Advance(it, edgeSize_ - 1).Valid())
			{
				return MakePositiveEdge(it.GetPosition());
			}

			return Edge();
		}
		
		int CountEquivalentEdges(const Edge & edge) const;
		int ListEdges(const Vertex & v, std::vector<Edge> & edge);
		int ListEdgesSeparate(const Vertex & v, std::vector<std::vector<Edge> > & edge);
		int FindEquivalentEdges(const Edge & edge, std::vector<Edge> & ret);
		DeBruijnGraph(const std::string & sequence, int edgeSize);

		DNASequence sequence;
	private:
		DISALLOW_COPY_AND_ASSIGN(DeBruijnGraph);				
		
		//Length of the k-mer that corresponds to the edge in graph
		const int edgeSize_;
		//Length of the (k - 1)-mer that corresponds to the vertex, equals edgeSize_ - 1
		const int vertexSize_;

		class KMerHashFunction
		{
		public:
			size_t operator ()(StrandConstIterator it) const
			{
				return it.GetHashCode(k_);
			}

			KMerHashFunction(int k): k_(k) {}
		private:
			int k_;
		};
		
		class KMerEqualTo
		{
		public:
			bool operator()(StrandConstIterator it1, StrandConstIterator it2) const
			{
				StrandConstIterator end1(Advance(it1, k_));
				return std::mismatch(it1, end1, it2).first == end1;
			}

			KMerEqualTo(int k): k_(k) {}
		private:
			int k_;
		};

		class IndexTransformer
		{
		public:
			StrandConstIterator operator()(int index) const
			{
				return sequence_->PositiveByIndex(index);
			}

			IndexTransformer(const DNASequence * sequence): sequence_(sequence){}
		private:
			const DNASequence * sequence_;
		};
		
		typedef IndexMultiSet<StrandConstIterator, IndexTransformer, KMerHashFunction, KMerEqualTo> KMerMultiSet;		
		//Container in which we store all edges
		KMerMultiSet edge_;
		int positiveVertexBound_;
		int negativeVertexBound_;

		void Init();
		void CalcBound();
		void PrintSet() const;
		void InvalidateAfter(int pos, bool erase);
		void InvalidateBefore(int pos, bool eraase);
		Edge MakePositiveEdge(int shift);
		Edge MakeNegativeEdge(int shift);
	};
}

#endif
