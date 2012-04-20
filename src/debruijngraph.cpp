#include "debruijngraph.h"

namespace SyntenyBuilder
{
	DeBruijnGraph::DeBruijnGraph(const std::string & sequence, int edgeSize):
		sequence(sequence,
			boost::bind(&DeBruijnGraph::Init, this),
			boost::bind(&DeBruijnGraph::InvalidateBefore, this, _1, _2),
			boost::bind(&DeBruijnGraph::InvalidateAfter, this, _1, _2)),
		edgeSize_(edgeSize), vertexSize_(edgeSize_ - 1),
		edge_(sequence.size(),
			IndexTransformer(&this->sequence), 
			KMerHashFunction(edgeSize),
			KMerEqualTo(edgeSize))
	{
		Init();
	}

	void DeBruijnGraph::CalcBound()
	{
		negativeVertexBound_ = Advance(sequence.PositiveBegin(), vertexSize_ - 1).GetPosition();
		positiveVertexBound_ = Advance(sequence.NegativeBegin(), vertexSize_ - 1).GetPosition();		
	}

	void DeBruijnGraph::Init()
	{
		CalcBound();
		edge_.Clear();
		for(int i = 0; i < static_cast<int>(sequence.Size() - edgeSize_ + 1); i++)
		{
			edge_.Insert(i);
		}
	}

	size_t DeBruijnGraph::ListEdges(const Vertex & v, std::vector<Edge> & edge)
	{
		edge.clear();
		std::vector<std::vector<Edge> > ret;
		ListEdgesSeparate(v, ret);
		for(size_t i = 0; i < ret.size(); i++)
		{
			std::copy(ret[i].begin(), ret[i].end(), std::back_inserter(edge));
		}

		return edge.size();
	}

	size_t DeBruijnGraph::CountEquivalentEdges(const Edge & edge) const
	{
		StrandConstIterator negative = edge.Direction() == positive ? 
			sequence.NegativeByIndex(edge.end_.GetPosition()) :
			sequence.PositiveByIndex(edge.end_.GetPosition());
		return edge_.Count(edge.start_) + edge_.Count(negative);
	}

	size_t DeBruijnGraph::FindEquivalentEdges(const Edge & edge, std::vector<Edge> & ret)
	{
		ret.clear();
		StrandConstIterator negIt = edge.Direction() == positive ? 
			sequence.NegativeByIndex(edge.end_.GetPosition()) :
			sequence.PositiveByIndex(edge.end_.GetPosition());
		std::vector<int> positive;
		std::vector<int> negative;
		edge_.Find(edge.start_, std::back_inserter(positive));
		edge_.Find(negIt, std::back_inserter(negative));
		std::transform(positive.begin(), positive.end(), std::back_inserter(ret), 
			boost::bind(&DeBruijnGraph::MakePositiveEdge, this, _1));
		std::transform(negative.begin(), negative.end(), std::back_inserter(ret), 
			boost::bind(&DeBruijnGraph::MakeNegativeEdge, this, _1));
		return ret.size();
	}
	
	size_t DeBruijnGraph::ListEdgesSeparate(const Vertex & v, std::vector<std::vector<Edge> > & edge)
	{
		edge.clear();
		std::string buf(edgeSize_, 't');
		v.Spell(buf.begin());
		DNASequence temp(buf);
		std::vector<int> positive;
		std::vector<int> negative;
		for(size_t i = 0; i < DNASequence::alphabet.size(); i++)
		{
			positive.clear();
			negative.clear();
			temp.NegativeBegin().AssignBase(DNASequence::alphabet[i]);
			edge_.Find(temp.PositiveBegin(), std::back_inserter(positive));
			edge_.Find(temp.NegativeBegin(), std::back_inserter(negative));
			if(positive.size() + negative.size() > 0)
			{
				edge.push_back(std::vector<Edge>());
				std::transform(positive.begin(), positive.end(), std::back_inserter(edge.back()),
					boost::bind(&DeBruijnGraph::MakePositiveEdge, this, _1));
				std::transform(negative.begin(), negative.end(), std::back_inserter(edge.back()),
					boost::bind(&DeBruijnGraph::MakeNegativeEdge, this, _1));
			}
		}

		return edge.size();
	}

	DeBruijnGraph::Edge DeBruijnGraph::MakePositiveEdge(int shift)
	{
		return Edge(this, sequence.PositiveByIndex(shift), positive);
	}

	DeBruijnGraph::Edge DeBruijnGraph::MakeNegativeEdge(int shift)
	{
		StrandConstIterator it = Advance(sequence.PositiveByIndex(shift), edgeSize_ - 1);
		return Edge(this, sequence.NegativeByIndex(it.GetPosition()), negative);
	}

	void DeBruijnGraph::InvalidateBefore(int pos, bool erase)
	{
		StrandIterator it = sequence.NegativeByIndex(pos);
		for(int i = 0; i < edgeSize_ && it.Valid(); i++, it++)
		{
			StrandIterator end = sequence.PositiveByIndex(it.GetPosition());
			if(Advance(end, edgeSize_ - 1).Valid())
			{
				edge_.Erase(it.GetPosition());
			}
		}
	}

	void DeBruijnGraph::InvalidateAfter(int pos, bool erase)
	{
		StrandIterator it = sequence.NegativeByIndex(pos);
		for(int i = 0; i < edgeSize_ && it.Valid(); i++, it++)
		{
			StrandIterator end = sequence.PositiveByIndex(it.GetPosition());
			if(Advance(end, edgeSize_ - 1).Valid())
			{
				edge_.Insert(it.GetPosition());
			}
		}

		if(erase)
		{
			CalcBound();
		}
	}
}