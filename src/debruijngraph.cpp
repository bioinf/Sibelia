#include "debruijngraph.h"


namespace SyntenyBuilder
{
	DeBruijnGraph::DeBruijnGraph(const std::string & sequence, int edgeSize):
		sequence_(sequence), 
		original_(sequence),
		edgeSize_(edgeSize), 
		vertexSize_(edgeSize - 1),
		edgeCount_(static_cast<int>(sequence_.size() - edgeSize_ + 1)),
		vertexCount_(static_cast<int>(sequence_.size() - edgeSize_ + 1)),
		edgeMultiSet_(edgeSize_, sequence.begin(), DELETED_CHARACTER) 
	{
		//Hash all k-mers (edges)
		for(int i = 0; i < edgeCount_; i++)
		{
			edgeMultiSet_.InsertShift(i);
		}

		//Store all positions in the original string
		position_.resize(sequence_.size());
		for(size_t i = 0; i < sequence_.size(); i++)
		{
			position_[i] = i;
		}
	}

	DeBruijnGraph::~DeBruijnGraph()
	{
		
	}

	size_t DeBruijnGraph::ListOutEdges(int shift, std::vector<EdgeList> & edgeList) const
	{		
		EdgeList ret;
		edgeList.clear();
		std::string directEdge(edgeSize_, ' ');
		std::string revCompEdge(edgeSize_, ' ');
		CopySequence(sequence_.begin() + shift, vertexSize_, directEdge.begin());		
		//Check all possible edges, including reverse complementary
		for(size_t i = 0; i < alphabet_.size(); i++)
		{
			*directEdge.rbegin() = alphabet_[i];			
			BuildReverseComplementary(directEdge.begin(), edgeSize_, revCompEdge.begin());			
			edgeMultiSet_.FindEquivalentShift(directEdge.begin(), ret.direct);
			edgeMultiSet_.FindEquivalentShift(revCompEdge.begin(), ret.revComp);
			if(ret.Size() > 0)
			{
				edgeList.push_back(ret);
			}
		}

		return edgeList.size();
	}	

	size_t DeBruijnGraph::FindEquivalentEdges(int shift, EdgeList & edgeList) const
	{		
		edgeList.Clear();
		std::string revCompEdge(edgeSize_, ' ');
		BuildReverseComplementary(sequence_.begin() + shift, edgeSize_, revCompEdge.begin());
		edgeMultiSet_.FindEquivalentShift(shift, edgeList.direct);
		edgeMultiSet_.FindEquivalentShift(revCompEdge.begin(), edgeList.revComp);
		return edgeList.Size();
	}	
	
	void DeBruijnGraph::DebugOutput(std::ostream & out) const
	{
		std::vector<int> temp;
		std::vector<EdgeList> edge;
		KMerMultiSet passed(vertexSize_, sequence_.begin(), DELETED_CHARACTER);
		for(int i = 0; i < vertexCount_; i++)
		{
			if(passed.FindEquivalentShift(i, temp) == 0)
			{
				passed.InsertShift(i);
				CopySequence(sequence_.begin() + i, vertexSize_, std::ostream_iterator<char>(out));
				out << std::endl;
				ListOutEdges(i, edge);
				for(size_t j = 0; j < edge.size(); j++)
				{
					for(size_t k = 0; k < edge[j].direct.size(); k++)
					{
						out << '\t';
						CopySequence(sequence_.begin() + edge[j].direct[k] + 1, vertexSize_, std::ostream_iterator<char>(out));
						out << std::endl;
					}
				}
			}
		}

		out << sequence_ << std::endl;
	}

	void DeBruijnGraph::CollapseBulge(int startVertex, int endVertex, const VisitData & smallBranch, const VisitData & bigBranch)
	{
		
	}

	int DeBruijnGraph::ProcessBifurcation(int vertex, const std::vector<EdgeList> & edge, int minCycleSize)
	{
		typedef google::sparse_hash_map<int,
										VisitData,
										KMerMultiSet::KMerHashFunction,
										KMerMultiSet::KMerEqualTo>
										VertexVisitMap;
		VertexVisitMap reach(minCycleSize, 
			KMerMultiSet::KMerHashFunction(vertexSize_, sequence_.begin(), DELETED_CHARACTER),
			KMerMultiSet::KMerEqualTo(vertexSize_, sequence_.begin(), DELETED_CHARACTER));

		int ret = 0;
	/*	for(int step = 1; step < minCycleSize && ret == 0; step++)
		{
			for(size_t i = 0; i < edge.size() && ret == 0; i++)
			{
				for(size_t j = 0; j < edge[i].direct.size() && ret == 0; j++)
				{
					int nextVertex = edge[i].direct[j] + step;
					VisitData nowVisit(i, edge[i].direct[j], true, step);
					VertexVisitMap::iterator it = reach.find(nextVertex);
					if(it != reach.end())
					{
						if(it->second != nowVisit && it->second.step + step < minCycleSize)
						{
							CollapseBulge(vertex, nextVertex, it->second, nowVisit);
						}
					}
					else
					{
						reach.insert(std::make_pair(nextVertex, nowVisit));
					}
				}

				for(size_t j = 0; j < edge[i].revComp.size() && ret == 0; j++)
				{
					int nextVertex = edge[i].direct[j] - step;
					VisitData nowVisit(i, edge[i].direct[j], false, step);
					VertexVisitMap::iterator it = reach.find(nextVertex);
					if(it != reach.end())
					{
						if(it->second != nowVisit && it->second.step + step < minCycleSize)
						{
							CollapseBulge(vertex, nextVertex, it->second, nowVisit);
						}
					}
					else
					{
						reach.insert(std::make_pair(nextVertex, nowVisit));
					}
				}
			}
		}*/

		return ret;
	}

	int DeBruijnGraph::RemoveBulges(int minCycleSize)
	{		
		int cyclesRemoved = 0;
		std::vector<EdgeList> edge;		
		for(int i = 0; i < vertexCount_; i++)
		{
			while(ListOutEdges(i, edge) > 1) //A bifurcation is found
			{
				cyclesRemoved += ProcessBifurcation(i, edge, minCycleSize);
			}
		}

		return cyclesRemoved;
	}

	void DeBruijnGraph::ClearSequence()
	{
		sequence_.erase(std::remove(sequence_.begin(), sequence_.end(), DELETED_CHARACTER), sequence_.end());
		position_.erase(std::remove(position_.begin(), position_.end(), DELETED_CHARACTER), position_.end());

		edgeMultiSet_.Clear();
		edgeCount_ = static_cast<int>(sequence_.size() - edgeSize_ + 1);
		vertexCount_ = static_cast<int>(sequence_.size() - vertexSize_ + 1);
		for(int i = 0; i < edgeCount_; i++)
		{
			edgeMultiSet_.InsertShift(i);
		}
	}

	void DeBruijnGraph::Simplify(int minCycleSize)
	{
		RemoveBulges(minCycleSize);
		ClearSequence();
	}

	int DeBruijnGraph::PassNonBranchingPathForward(const EdgeList & edge, std::vector<char> & visit) const
	{
		int travel = 1;
		bool fail = false;
		int size = static_cast<int>(sequence_.size());
		do
		{
			int consensus = edge.direct[0] + edgeSize_ + travel - 1;
			fail = !(consensus < size);
			for(size_t i = 0; i < edge.direct.size(); i++)
			{				
				int shift = edge.direct[i] + edgeSize_ + travel - 1;				
				fail = fail || !(shift < size) || visit[edge.direct[i] + travel] || (sequence_[shift] != sequence_[consensus]);
			}

			for(size_t i = 0; i < edge.revComp.size(); i++)
			{
				int shift = edge.revComp[i] - travel;				
				fail = fail || shift < 0 || visit[edge.revComp[i] - travel] || (complement_[sequence_[shift]] != sequence_[consensus]);
			}
			
			if(!fail)
			{
				edge.Invalidate(visit, travel++);
			}
		}
		while(!fail);

		return travel - 1;
	}

	int DeBruijnGraph::PassNonBranchingPathBackward(const EdgeList & edge, std::vector<char> & visit) const
	{
		int travel = 1;
		bool fail = false;
		int size = static_cast<int>(sequence_.size());
		do
		{
			int consensus = edge.direct[0] - travel;
			fail = !(consensus >= 0);
			for(size_t i = 0; i < edge.direct.size(); i++)
			{				
				int shift = edge.direct[i] - travel;				
				fail = fail || shift < 0 || visit[edge.direct[i] - travel] || (sequence_[shift] != sequence_[consensus]);
			}

			for(size_t i = 0; i < edge.revComp.size(); i++)
			{
				int shift = edge.revComp[i] + edgeSize_ + travel - 1;			
				fail = fail || shift < 0 || !(shift < size) || visit[edge.revComp[i] + travel] ||
					(complement_[sequence_[shift]] != sequence_[consensus]);
			}

			if(!fail)
			{
				edge.Invalidate(visit, -travel++);
			}
		}
		while(!fail);

		return travel - 1;
	}

	void DeBruijnGraph::ListNonBranchingPaths(std::ostream & out) const
	{
		EdgeList edgeList;
		std::vector<char> visit(sequence_.size(), false);
		std::vector<std::pair<size_t, int> > interestEdge;		
		for(int i = 0; i < edgeCount_; i++)
		{	
			if(FindEquivalentEdges(i, edgeList) > 1 && edgeList.direct.size() > 0) //We have found and edge that occur many times
			{
				interestEdge.push_back(std::make_pair(edgeList.Size(), i));
			}
		}

		//Sort any possible point of interest to process paths with less multiplicity first
		std::sort(interestEdge.begin(), interestEdge.end());									
		for(size_t j = 0; j < interestEdge.size(); j++)
		{
			int i = interestEdge[j].second;
			if(FindEquivalentEdges(i, edgeList) > 0 && edgeList.Filter(visit) > 1)			
			{	
				int pathLengthForward = PassNonBranchingPathForward(edgeList, visit);
				int pathLengthBackward = PassNonBranchingPathBackward(edgeList, visit);
				int pathLength = pathLengthForward + pathLengthBackward;
				int consensusLength = pathLength  + edgeSize_;
				edgeList.Invalidate(visit, 0);
				
				out << "Consensus: ";
				std::copy(sequence_.begin() + (i - pathLengthBackward), sequence_.begin() + (i - pathLengthBackward + consensusLength),
					std::ostream_iterator<char>(out));
				out << std::endl << "Original sequences: " << std::endl;	
				//List found sequences that occur on the direct strand
				for(size_t l = 0; l < edgeList.direct.size(); l++)
				{
					int startPos = position_[edgeList.direct[l] - pathLengthBackward];
					int endPos = position_[edgeList.direct[l] + pathLengthForward] + edgeSize_;
					out << startPos + 1 << ":" << endPos + 1 << " ";
					std::copy(original_.begin() + startPos, original_.begin() + endPos, std::ostream_iterator<char>(out));
					out << std::endl;
				}

				//List found sequences that occur on the reverse-complementary strand
				for(size_t l = 0; l < edgeList.revComp.size(); l++)
				{
					int startPos = position_[edgeList.revComp[l] - pathLengthForward];
					int endPos = position_[edgeList.revComp[l] + pathLengthBackward] + edgeSize_;
					//Compute coordinates of the string on the second strand
					int revStartPos = original_.size() - endPos;
					int revEndPos = revStartPos + endPos - startPos;

					out << - (revStartPos + 1) << ":" << -(revEndPos + 1) << " ";
					std::string buf(endPos - startPos, ' ');
					BuildReverseComplementary(original_.begin() + startPos, endPos - startPos, buf.begin());
					std::copy(buf.begin(), buf.end(), std::ostream_iterator<char>(out));
					out << std::endl;
				}
			}
		}
	}
}