#include "debruijngraph.h"

namespace SyntenyBuilder
{
	bool DeBruijnGraph::Edge::IsNull() const
	{
		return graph_ == 0;
	}

	size_t DeBruijnGraph::Edge::GetHashCode() const
	{
		return startPos_;
	}

	DeBruijnGraph::DirectStrandIterator DeBruijnGraph::Edge::GetPosition() const
	{
		return graph_->SequenceIterator(graph_->sequence_.begin() + startPos_);
	}

	DeBruijnGraph::Edge::Edge(int shift, const DeBruijnGraph * graph): graph_(graph)
	{
		startPos_ = shift;
		endPos_ = graph_->GetNextPosForward(startPos_, graph_->edgeSize_ - 1);
	}

	DeBruijnGraph::VertexPtr DeBruijnGraph::Edge::StartVertex() const
	{
		return VertexPtr(new Vertex(startPos_, graph_));
	}

	DeBruijnGraph::VertexPtr DeBruijnGraph::Edge::EndVertex() const
	{		
		return VertexPtr(new Vertex(graph_->GetNextPosForward(startPos_, 1), graph_));
	}

	char DeBruijnGraph::Edge::Spell() const
	{		
		return graph_->sequence_[endPos_];
	}

	/*
	void DeBruijnGraph::Edge::Spell(std::string::iterator it) const
	{
		this->StartVertex().SpellDirect(it);
		*it++ = this->Spell();
	}*/

	DeBruijnGraph::EdgePtr DeBruijnGraph::Edge::NextEdge() const
	{
		if(graph_->GetNextPosForward(endPos_) != END_OF_SEQUENCE)
		{
			return EdgePtr(new Edge(graph_->GetNextPosForward(startPos_), graph_));
		}

		return EdgePtr(new Edge());
	}

	DeBruijnGraph::ComplementaryEdge::ComplementaryEdge(int shift, const DeBruijnGraph * graph): Edge(shift, graph)
	{
		
	}

	DeBruijnGraph::VertexPtr DeBruijnGraph::ComplementaryEdge::StartVertex() const
	{
		return VertexPtr(new ComplementaryVertex(graph_->GetNextPosForward(startPos_, 1), graph_));
	}

	DeBruijnGraph::VertexPtr DeBruijnGraph::ComplementaryEdge::EndVertex() const
	{
		return VertexPtr(new ComplementaryVertex(startPos_, graph_));
	}

	char DeBruijnGraph::ComplementaryEdge::Spell() const
	{
		return complement_[graph_->sequence_[startPos_]];
	}

	DeBruijnGraph::EdgePtr DeBruijnGraph::ComplementaryEdge::NextEdge() const
	{
		if(graph_->GetNextPosBackward(startPos_) != END_OF_SEQUENCE)
		{
			return EdgePtr(new Edge(graph_->GetNextPosBackward(endPos_), graph_));
		}

		return EdgePtr(new ComplementaryEdge());
	}
}