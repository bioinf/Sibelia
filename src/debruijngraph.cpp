#include "debruijngraph.h"

namespace SyntenyBuilder
{
	DeBruijnGraph::DeBruijnGraph(const std::string & sequence, int edgeSize):
		sequence_(sequence), 
		original_(sequence),
		edgeSize_(edgeSize), 
		vertexSize_(edgeSize - 1),
		vertexSet_(vertexSize_, sequence.begin(), DELETED_CHARACTER),
		edgeMultiSet_(edgeSize_, sequence.begin(), DELETED_CHARACTER)
	{
		Init();
		//Store all positions in the original string
		position_.resize(sequence_.size());
		for(size_t i = 0; i < sequence_.size(); i++)
		{
			position_[i] = i;
		}
	}

	StringIterator DeBruijnGraph::SequenceIterator(std::string::iterator it)
	{
		return StringIterator(it, sequence_.end(), DELETED_CHARACTER);
	}

	StringConstIterator DeBruijnGraph::SequenceIterator(std::string::const_iterator it) const
	{
		return StringConstIterator(it, sequence_.end(), DELETED_CHARACTER);
	}

	StringReverseIterator DeBruijnGraph::SequenceReverseIterator(std::string::iterator it)
	{
		return StringReverseIterator(std::string::reverse_iterator(it), sequence_.rend(), DELETED_CHARACTER);
	}

	StringConstReverseIterator DeBruijnGraph::SequenceReverseIterator(std::string::const_iterator it) const
	{
		return StringConstReverseIterator(std::string::const_reverse_iterator(it), sequence_.rend(), DELETED_CHARACTER);
	}

	size_t DeBruijnGraph::GetSize() const
	{
		return sequence_.size();
	}

	size_t DeBruijnGraph::GetEdgeSize() const
	{
		return edgeSize_;
	}

	size_t DeBruijnGraph::GetVertexSize() const
	{
		return vertexSize_;
	}
	
	DeBruijnGraph::DirectStrandIterator DeBruijnGraph::SequenceBegin() const
	{
		return SequenceIterator(sequence_.begin());
	}

	DeBruijnGraph::DirectStrandIterator DeBruijnGraph::SequenceEnd() const
	{
		return SequenceIterator(sequence_.end());
	}

	DeBruijnGraph::DirectStrandIterator DeBruijnGraph::LastVertex() const
	{
		int shift = GetNextPosBackward(sequence_.size(), vertexSize_ - 1);
		return SequenceIterator(sequence_.begin() + shift);
	}

	DeBruijnGraph::DirectStrandIterator DeBruijnGraph::LastEdge() const
	{
		int shift = GetNextPosBackward(sequence_.size(), vertexSize_ - 1);
		return SequenceIterator(sequence_.begin() + shift);
	}

	DeBruijnGraph::VertexPtr DeBruijnGraph::ConstructVertex(DirectStrandIterator it) const
	{
		std::vector<int> ret;
		if(vertexSet_.FindEquivalentShift(it.CurrentPosition() - sequence_.begin(), ret) > 0)
		{
			return VertexPtr(new Vertex(it.CurrentPosition() - sequence_.begin(), this));
		}

		return VertexPtr(new ComplementaryVertex(it.CurrentPosition() - sequence_.begin(), this));
	}

	DeBruijnGraph::Edge DeBruijnGraph::ConstructEdge(DirectStrandIterator it) const
	{
		return Edge(it.CurrentPosition() - sequence_.begin(), this);
	}

	size_t DeBruijnGraph::FindEquivalentEdges(const Edge & edge, EdgeListPair & edgeList) const
	{	
		edgeList.Clear();
		std::vector<int> direct;
		std::vector<int> comp;
		std::string compEdge(edgeSize_, ' ');
		BuildReverseComplementary(SequenceIterator(sequence_.begin() + edge.startPos_), edgeSize_, compEdge.begin());
		
		edgeMultiSet_.FindEquivalentShift(edge.startPos_, direct);
		edgeMultiSet_.FindEquivalentShift(compEdge.begin(), comp);
		return edgeList.Assign(direct, comp, this);
	}	

	void DeBruijnGraph::Init()
	{
		edgeMultiSet_.Clear();
		edgeCount_ = static_cast<int>(sequence_.size() - edgeSize_ + 1);
		vertexCount_ = static_cast<int>(sequence_.size() - vertexSize_ + 1);
		for(int i = 0; i < edgeCount_; i++)
		{
			edgeMultiSet_.InsertShift(i);
		}

		std::vector<int> ret;
		std::string buf(vertexSize_, ' ');
		for(int i = 0; i < vertexCount_; i++)
		{
			BuildReverseComplementary(sequence_.begin() + i, vertexSize_, buf.begin());
			if(vertexSet_.FindEquivalentShift(i, ret) == 0 && vertexSet_.FindEquivalentShift(buf.begin(), ret) == 0)
			{
				vertexSet_.InsertShift(i);
			}			
		}
	}

	void DeBruijnGraph::ClearSequence()
	{
		sequence_.erase(std::remove(sequence_.begin(), sequence_.end(), DELETED_CHARACTER), sequence_.end());
		position_.erase(std::remove(position_.begin(), position_.end(), DELETED_CHARACTER), position_.end());
		Init();
	}

	int DeBruijnGraph::GetNextPosForward(int shift, int steps) const
	{
 		StringConstIterator it = SequenceIterator(sequence_.begin() + shift);
		std::advance(it, steps);
		if(it.CurrentPosition() != sequence_.end())
		{
			return it.CurrentPosition() - sequence_.begin();
		}

		return END_OF_SEQUENCE;
	}

	int DeBruijnGraph::GetNextPosBackward(int shift, int steps) const
	{
		if(steps > 0)
		{
			StringConstReverseIterator it = SequenceReverseIterator(sequence_.begin() + shift);
			std::advance(it, steps - 1);
			if(it.CurrentPosition() != sequence_.rend())
			{
				return sequence_.rend() - it.CurrentPosition() - 1;
			}

			return END_OF_SEQUENCE;
		}		

		return shift;
	}
	
	//int DeBruijnGraph::ProcessBifurcation(int vertex, const std::vector<EdgeList> & edge, int minCycleSize)
	//{
	//	typedef google::sparse_hash_map<int,
	//									VisitData,
	//									KMerMultiSet::KMerHashFunction,
	//									KMerMultiSet::KMerEqualTo>
	//									VertexVisitMap;
	//	VertexVisitMap reach(minCycleSize, 
	//		KMerMultiSet::KMerHashFunction(vertexSize_, sequence_.begin(), DELETED_CHARACTER),
	//		KMerMultiSet::KMerEqualTo(vertexSize_, sequence_.begin(), DELETED_CHARACTER));

	//	int ret = 0;
	///*	for(int step = 1; step < minCycleSize && ret == 0; step++)
	//	{
	//		for(size_t i = 0; i < edge.size() && ret == 0; i++)
	//		{
	//			for(size_t j = 0; j < edge[i].direct.size() && ret == 0; j++)
	//			{
	//				int nextVertex = edge[i].direct[j] + step;
	//				VisitData nowVisit(i, edge[i].direct[j], true, step);
	//				VertexVisitMap::iterator it = reach.find(nextVertex);
	//				if(it != reach.end())
	//				{
	//					if(it->second != nowVisit && it->second.step + step < minCycleSize)
	//					{
	//						CollapseBulge(vertex, nextVertex, it->second, nowVisit);
	//					}
	//				}
	//				else
	//				{
	//					reach.insert(std::make_pair(nextVertex, nowVisit));
	//				}
	//			}

	//			for(size_t j = 0; j < edge[i].revComp.size() && ret == 0; j++)
	//			{
	//				int nextVertex = edge[i].direct[j] - step;
	//				VisitData nowVisit(i, edge[i].direct[j], false, step);
	//				VertexVisitMap::iterator it = reach.find(nextVertex);
	//				if(it != reach.end())
	//				{
	//					if(it->second != nowVisit && it->second.step + step < minCycleSize)
	//					{
	//						CollapseBulge(vertex, nextVertex, it->second, nowVisit);
	//					}
	//				}
	//				else
	//				{
	//					reach.insert(std::make_pair(nextVertex, nowVisit));
	//				}
	//			}
	//		}
	//	}*/

	//	return ret;
	//}

	//int DeBruijnGraph::RemoveBulges(int minCycleSize)
	//{		
	//	int cyclesRemoved = 0;
	//	std::vector<EdgeList> edge;		
	//	for(int i = 0; i < vertexCount_; i++)
	//	{			
	//		while(ListOutEdges(i, edge) > 1) //A bifurcation has been found
	//		{
	//			cyclesRemoved += ProcessBifurcation(i, edge, minCycleSize);
	//		}
	//	}

	//	return cyclesRemoved;
	//}

}