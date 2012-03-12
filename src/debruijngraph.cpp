#include "debruijngraph.h"

namespace SyntenyBuilder
{
	const std::string DeBruijnGraph::alphabet_ = "acgt";

	size_t DeBruijnGraph::KMerHashFunction::operator()(std::string::iterator it) const
	{
		size_t base = 1;
		size_t hash = 0;
		std::string::const_reverse_iterator end(it);
		for(std::string::const_reverse_iterator iter(it + k_); iter != end; iter++)
		{
			hash += *iter * base;
			base *= HASH_BASE;
		}

		return hash;
	}

	DeBruijnGraph::DeBruijnGraph(const std::string & sequence, int k):
		sequence_(sequence), k_(k), kMerTable_(sequence.size() - k_, KMerHashFunction(k_), KMerEqualTo(k_))
	{
		//Hash all k-mers
		for(std::string::iterator it = sequence_.begin(); it < sequence_.end() - k_ + 1; it++)
		{
			kMerTable_[it]++;
		}
	}

	DeBruijnGraph::~DeBruijnGraph()
	{
		
	}

	void DeBruijnGraph::ListOutEdges(std::string::iterator shift, std::vector<std::string::iterator> & edge) const
	{
		edge.clear();
		std::string temp(shift, shift + k_ - 1);
		temp.push_back(' ');
		//Check all possible edges
		for(size_t i = 0; i < alphabet_.size(); i++)
		{
			temp.back() = alphabet_[i];
			ShiftTable::iterator found = kMerTable_.find(temp.begin());
			if(found != kMerTable_.end())
			{
				edge.push_back(found->first);
			}
		}
	}

	//travel through the graph
	void DeBruijnGraph::DebugOutput(std::ostream & out) const
	{
		size_t vertexLength = k_ - 1;
		std::vector<std::string::iterator> edge;
		for(ShiftTable::iterator it = kMerTable_.begin(); it != kMerTable_.end(); it++)
		{
			//if this edge has not been visited
			if(it->second > 0)
			{
				std::string::iterator shift = it->first;
				std::copy(shift, shift + vertexLength, std::ostream_iterator<char>(out));
				out << std::endl;
				ListOutEdges(shift, edge);
				for(size_t i = 0; i < edge.size(); i++)
				{
					ShiftTable::iterator jt = kMerTable_.find(edge[i]);
					out << '\t';
					std::copy(jt->first + 1, jt->first + 1 + vertexLength, std::ostream_iterator<char>(out));
					out << ' ' << jt->second << std::endl;
					//Mark this edge as visited
					jt->second = -jt->second;
				}
			}
		}

		//Return everything to the initial state
		for(ShiftTable::iterator it = kMerTable_.begin(); it != kMerTable_.end(); it++)
		{
			it->second = -it->second;
		}
	}
}