#include "debruijngraph.h"

namespace SyntenyBuilder
{
	const size_t DeBruijnGraph::HASH_TABLE_MAX_SIZE = 1 << 20;

	DeBruijnGraph::DeBruijnGraph(const std::string & sequence, int z):
		sequence_(sequence), k_(z - 1), vertexTable_(std::min(sequence.size(), HASH_TABLE_MAX_SIZE)), basePowK(1)
	{
		for(uint64_t i = 0; i < k_ - 1; i++)
		{
			basePowK *= HASH_BASE;
		}

		uint64_t nowHash = CalculateKMerHash(sequence_.begin());
		int nowVertex = InsertVertex(nowHash, 0);
		for(size_t i = 0; i < sequence_.size() - k_ + 1; i++)
		{			
			if(i != sequence_.size() - k_)
			{
				uint64_t nextHash = ShiftKMerHash(nowHash, sequence_[i], sequence_[i + k_]);
				int nextVertex = InsertVertex(nextHash, i + 1);
				bool found = false;
				for(size_t l = 0; l < adjacencyList_[nowVertex].size(); l++)
				{
					if(adjacencyList_[nowVertex][l].first == nextVertex)
					{
						found = true;
						adjacencyList_[nowVertex][l].second++;
						break;
					}
				}

				if(!found)
				{
					adjacencyList_[nowVertex].push_back(std::make_pair(nextVertex, 1));
				}

				nowHash = nextHash;
				nowVertex = nextVertex;
			}
		}
	}

	DeBruijnGraph::~DeBruijnGraph()
	{
		
	}

	void DeBruijnGraph::DebugOutput(std::ostream & out) const
	{
		for(size_t i = 0; i < vertex_.size(); i++)
		{
			std::string::const_iterator it = sequence_.begin() + vertex_[i];
			std::copy(it, it + k_, std::ostream_iterator<char>(out));
			out << ":" << std::endl;
			for(size_t j = 0; j < adjacencyList_[i].size(); j++)
			{
				out << '\t';
				it = sequence_.begin() + vertex_[adjacencyList_[i][j].first];
				std::copy(it, it + k_, std::ostream_iterator<char>(out));
				out << ' ' << adjacencyList_[i][j].second << std::endl;			
			}

			out << "---" << std::endl;
		}
	}

	int DeBruijnGraph::FindVertex(uint64_t hash, std::string::const_iterator it) const
	{
		int ret = VERTEX_NOT_FOUND;
		size_t bucket = hash % vertexTable_.size();
		for(size_t i = 0; i < vertexTable_[bucket].size(); i++)
		{
			size_t nextShift = vertex_[vertexTable_[bucket][i]];
			if(std::mismatch(it, it + k_, sequence_.begin() + nextShift).first == it + k_)
			{
				ret = vertexTable_[bucket][i];
				break;
			}
		}

		return ret;
	}

	int DeBruijnGraph::InsertVertex(uint64_t hash, size_t shift)
	{
		int found = FindVertex(hash, sequence_.begin() + shift);
		if(found == VERTEX_NOT_FOUND)
		{
			found = static_cast<int>(vertex_.size());
			vertex_.push_back(static_cast<int>(shift));
			adjacencyList_.push_back(std::vector<std::pair<int, int> >());
			vertexTable_[hash % vertexTable_.size()].push_back(found);
		}

		return found;
	}

	uint64_t DeBruijnGraph::CalculateKMerHash(std::string::const_iterator it) const
	{
		uint64_t base = 1;
		uint64_t hash = 0;

		std::string::const_reverse_iterator end(it);
		for(std::string::const_reverse_iterator iter(it + k_); iter != end; iter++ )
		{
			hash += *iter * base;
			base *= HASH_BASE;
		}

		return hash;
	}

	uint64_t DeBruijnGraph::ShiftKMerHash(uint64_t hash, char firstChar, char nextChar) const
	{
		return (hash - firstChar * basePowK) * HASH_BASE + nextChar;
	}

}