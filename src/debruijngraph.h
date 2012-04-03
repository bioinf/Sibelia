#ifndef _DE_BRUIJN_GRAPH_H_
#define _DE_BRUIJN_GRAPH_H_

#include "kmermultiset.h"
#include "syntenyutility.h"

namespace SyntenyBuilder
{
	//This class represents de Bruijn graph. It stores edges (k-mer) in a hash
	//table. Further, more efficient implementation can be provided. The code
	//is quite messy, maybe later it will be remastered with more incapsulation
	//and OOP to remove some duplications from the code.
	class DeBruijnGraph
	{
	public:
		DeBruijnGraph(const std::string & sequence, int edgeSize);
		void Simplify(int minCycleSize);
		void DebugOutput(std::ostream & out) const;
		void ListNonBranchingPaths(std::ostream & out) const;
		~DeBruijnGraph();
	private:
		DISALLOW_COPY_AND_ASSIGN(DeBruijnGraph);				
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
		KMerMultiSet edgeMultiSet_;

		//FindEquivalentEdges returns list of the equivalent edges (k-mers). Reverse
		//complementary variants of the same edge are considered equivalent since we
		//glue reverse complementary vertices.
		size_t FindEquivalentEdges(int shift, EdgeList & edgeList) const;
		//ListOutEdges returns list of the edges (k-mers) that leave vertex located at 
		//position "shift". Edges are grouped by equivalency -- if two edges are in the
		//same structure EdgeList, they are equivalent
		size_t ListOutEdges(int shift, std::vector<EdgeList> & edgeList) const;				
		//This method returns the length of the nonbranching path and marks edges that
		//are visited.
		int PassNonBranchingPathForward(const EdgeList & edgeList, std::vector<char> & visit) const;
		int PassNonBranchingPathBackward(const EdgeList & edgeList, std::vector<char> & visit) const;
		//Procedures that remove small cycles.
		int RemoveBulges(int minCycleSize);
		int ProcessBifurcation(int shift, const std::vector<EdgeList> & edge, int minCycleSize);
		void CollapseBulge(int startVertex, int endVertex, const VisitData & smallBranch, const VisitData & bigBranch);
		//Method for clearing sequence_ from deleted characters.
		void ClearSequence();
	};
}

#endif
