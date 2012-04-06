#include "debruijngraph.h"

namespace SyntenyBuilder
{
	DeBruijnGraph::Vertex::Vertex(int shift, const DeBruijnGraph * graph): shift_(shift), graph_(graph)
	{		
		lastChar_ = graph_->GetNextPosForward(shift_, graph_->vertexSize_ - 1);
		if(lastChar_ == END_OF_SEQUENCE)
		{
			graph_ = 0;
		}
	}

	size_t DeBruijnGraph::Vertex::GetHashCode() const
	{
		KMerMultiSet::KMerHashFunction f(graph_->vertexSize_, graph_->sequence_.begin(), DELETED_CHARACTER);
		return f(shift_);
	}

	void DeBruijnGraph::Vertex::Spell(std::string::iterator it) const
	{
		std::copy_n(graph_->SequenceIterator(graph_->sequence_.begin() + shift_), graph_->vertexSize_, it);
	}

	bool DeBruijnGraph::Vertex::IsNull() const
	{
		return graph_ == 0;
	}

	bool DeBruijnGraph::Vertex::operator != (const Vertex & vertex) const
	{
		return !(*this == vertex);
	}

	bool DeBruijnGraph::Vertex::operator == (const Vertex & vertex) const
	{
		if(IsNull() || vertex.IsNull())
		{
			return false;
		}

		std::string buf1(graph_->vertexSize_, ' ');
		std::string buf2(graph_->vertexSize_, ' ');
		Spell(buf1.begin());
		Spell(buf2.begin());
		return buf1 == buf2;
	}

	size_t DeBruijnGraph::Vertex::ListOutEdges(std::vector<DeBruijnGraph::EdgePtr> & edgeList) const
	{	
		edgeList.clear();
		std::vector<int> direct;
		std::vector<int> revComp;
		std::string directEdge(graph_->edgeSize_, ' ');
		std::string revCompEdge(graph_->edgeSize_, ' ');
		Spell(directEdge.begin());
		//Check all possible edges, including reverse complementary
		for(size_t i = 0; i < alphabet_.size(); i++)
		{
			*directEdge.rbegin() = alphabet_[i];			
			BuildReverseComplementary(directEdge.begin(), graph_->edgeSize_, revCompEdge.begin());			
			graph_->edgeMultiSet_.FindEquivalentShift(directEdge.begin(), direct);
			graph_->edgeMultiSet_.FindEquivalentShift(revCompEdge.begin(), revComp);
			if(direct.size() + revComp.size() > 0)
			{
				for(size_t i = 0; i < direct.size(); i++)
				{
					edgeList.push_back(EdgePtr(new Edge(direct[i], graph_)));
				}

				for(size_t i = 0; i < revComp.size(); i++)
				{
					edgeList.push_back(EdgePtr(new ComplementaryEdge(revComp[i], graph_)));
				}
			}
		}

		return edgeList.size();
	}
	/*
	size_t DeBruijnGraph::Vertex::ListOutEdges(std::vector<DeBruijnGraph::EdgeListPair> & edgeList) const
	{	
		edgeList.clear();
		std::vector<int> direct;
		std::vector<int> revComp;
		std::string directEdge(graph_->edgeSize_, ' ');
		std::string revCompEdge(graph_->edgeSize_, ' ');
		Spell(directEdge.begin());
		//Check all possible edges, including reverse complementary
		for(size_t i = 0; i < alphabet_.size(); i++)
		{
			*directEdge.rbegin() = alphabet_[i];			
			BuildReverseComplementary(directEdge.begin(), graph_->edgeSize_, revCompEdge.begin());			
			graph_->edgeMultiSet_.FindEquivalentShift(directEdge.begin(), direct);
			graph_->edgeMultiSet_.FindEquivalentShift(revCompEdge.begin(), revComp);
			if(direct.size() + revComp.size() > 0)
			{
				edgeList.push_back(EdgeListPair());
				edgeList.back().Assign(direct, revComp, graph_);
			}
		}

		return edgeList.size();
	}*/

	size_t DeBruijnGraph::ComplementaryVertex::GetHashCode() const
	{
		std::string buf(graph_->vertexSize_, ' ');
		BuildReverseComplementary(graph_->SequenceIterator(graph_->sequence_.begin() + shift_), graph_->vertexSize_, buf.begin());
		KMerMultiSet::KMerHashFunction f(graph_->vertexSize_, buf.begin(), DELETED_CHARACTER);
		return f(0);
	}

	void DeBruijnGraph::ComplementaryVertex::Spell(std::string::iterator it) const
	{
		StringConstIterator next = ++graph_->SequenceIterator(graph_->sequence_.begin() + lastChar_);
		StringConstReverseIterator rev = graph_->SequenceReverseIterator(next.CurrentPosition());
		std::copy_n(RevCompIterator(rev), graph_->vertexSize_, it);
	}
}