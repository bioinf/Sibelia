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
		static const std::string alphabet_;

		class KMerEqualTo
		{
		public:
			KMerEqualTo(int k): k_(k) { }
			bool operator()(std::string::iterator it1, std::string::iterator it2) const
			{
				return std::mismatch(it1, it1 + k_, it2).first == it1 + k_;
			}
		private:
			int k_;
		};

		class KMerHashFunction
		{
		public:
			KMerHashFunction(int k): k_(k){ }
			size_t operator()(std::string::iterator it) const;
		private:
			int k_;
			static const size_t HASH_BASE = 57;
		};
		
		void ListOutEdges(std::string::iterator shift, std::vector<std::string::iterator> & edge) const;

		const int k_;
		std::string sequence_;
		typedef google::sparse_hash_map<std::string::iterator, int, KMerHashFunction, KMerEqualTo> ShiftTable; 
		mutable ShiftTable kMerTable_;
	};
}

#endif
