#ifndef _DE_BRUIJN_GRAPH_H_
#define _DE_BRUIJN_GRAPH_H_

#include "common.h"

namespace SyntenyBuilder
{
	//This class represents de Bruijn graph. It stores vertices in a hash
	//table with direct access. Further, more efficient implementation
	//can be provided.
	class DeBruijnGraph
	{
	public:
		DeBruijnGraph(const std::string & sequence, int k);
		void DebugOutput(std::ostream & out) const;
		void ListNonBranchingPaths(std::ostream & out) const;
		~DeBruijnGraph();
	private:
		DISALLOW_COPY_AND_ASSIGN(DeBruijnGraph);

		static const uint64_t HASH_BASE = 57;
		static const int VERTEX_NOT_FOUND = -1;
		static const size_t HASH_TABLE_MAX_SIZE;		

		//Each vertex is represented as a k-mer
		const int k_;
		//Whole sequence
		std::string sequence_;
		//List of the vertices. vertex_[i] = shift in original sequence,
		//that represents contents (k-mer string) of the i-th vertex
		std::vector<size_t> vertex_;
		//Hash table for storing vertices
		std::vector<std::vector<int> > vertexTable_;
		uint64_t basePowK;
		//Adjacency list of the graph
		std::vector<std::vector<std::pair<int, int> > > adjacencyList_;

		//Methods for dealing with vertices/vertex numbers
		int InsertVertex(uint64_t hashValue, size_t shift);		
		int FindVertex(uint64_t hashValue, std::string::const_iterator it) const;
		uint64_t CalculateKMerHash(std::string::const_iterator it) const;
		uint64_t ShiftKMerHash(uint64_t hashValue, char firstChar, char nextChar) const;
		//Other methods
	};
}

#endif
