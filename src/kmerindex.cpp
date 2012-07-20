#include "kmerindex.h"
#include "debruijngraph.h"

namespace SyntenyBuilder
{
	void KMerIndex::SetupIndex(size_t k)
	{
		k_ = k;
		delete kmer_;
		kmer_ = new KMerMultiSet(sequence_->Size(),
			IndexTransformer(sequence_), 
			KMerHashFunction(k),
			KMerEqualTo(k));

		IndexKMers(sequence_->PositiveBegin(), sequence_->PositiveRightEnd());
	}

	void KMerIndex::IndexKMers(StrandIterator start, StrandIterator end)
	{
		for(; start != end; ++start)
		{
			kmer_->Insert(start.GetPosition());
		}
	}

	size_t KMerIndex::GetK() const
	{
		return k_;
	}

	/*

	int DeBruijnGraph::CountEquivalentEdges(const Edge & edge) const
	{
		StrandConstIterator negative = edge.Direction() == positive ? 
			sequence.NegativeByIndex(edge.end_.GetPosition()) :
			sequence.PositiveByIndex(edge.end_.GetPosition());
		return static_cast<int>(edge_->Count(edge.start_) + edge_->Count(negative));
	}
	*/

	size_t KMerIndex::ListEquivalentKmers(StrandIterator it, std::vector<StrandIterator> & ret) const
	{
		ret.clear();
		std::vector<size_t> positive;
		std::vector<size_t> negative;
		kmer_->Find(it, std::back_inserter(positive));
		it.Jump(k_);
		kmer_->Find(it.Invert(), std::back_inserter(negative));
		
		std::transform(positive.begin(), positive.end(), std::back_inserter(ret),
			boost::bind(

		return static_cast<int>(ret.size());
	}
	
	/*
	int DeBruijnGraph::ListEdgesSeparate(const Vertex & v, std::vector<std::vector<Edge> > & edge)
	{
		edge.clear();
		std::string buf(edgeSize_, 't');
		v.Spell(buf.begin());
		v.Spell(buf.begin());
		DNASequence temp(buf);
		std::vector<int> positive;
		std::vector<int> negative;
		for(size_t i = 0; i < DNASequence::alphabet.size(); i++)
		{
			positive.clear();
			negative.clear();
			temp.NegativeBegin().AssignBase(DNASequence::alphabet[i]);
			edge_->Find(temp.PositiveBegin(), std::back_inserter(positive));
			edge_->Find(temp.NegativeBegin(), std::back_inserter(negative));
			if(positive.size() + negative.size() > 0)
			{
				edge.push_back(std::vector<Edge>());
				std::transform(positive.begin(), positive.end(), std::back_inserter(edge.back()),
					boost::bind(&DeBruijnGraph::MakePositiveEdge, this, _1));
				std::transform(negative.begin(), negative.end(), std::back_inserter(edge.back()),
					boost::bind(&DeBruijnGraph::MakeNegativeEdge, this, _1));
			}
		}

		return static_cast<int>(edge.size());
	}

	DeBruijnGraph::Edge DeBruijnGraph::MakePositiveEdge(int shift)
	{
		return Edge(this, sequence.PositiveByIndex(shift), positive);
	}

	DeBruijnGraph::Edge DeBruijnGraph::MakeNegativeEdge(int shift)
	{
		StrandConstIterator it = AdvanceForward(sequence.PositiveByIndex(shift), edgeSize_ - 1);
		return Edge(this, sequence.NegativeByIndex(it.GetPosition()), negative);
	}

	void DeBruijnGraph::InvalidateBefore(int pos, bool erase)
	{
		StrandIterator it = sequence.NegativeByIndex(pos);
		for(int i = 0; i < edgeSize_ && it.Valid(); i++, it++)
		{
			StrandIterator end = sequence.PositiveByIndex(it.GetPosition());
			if(AdvanceForward(end, edgeSize_ - 1).Valid())
			{
				edge_->Erase(it.GetPosition());
			}
		}
	}

	void DeBruijnGraph::InvalidateAfter(int pos, bool erase)
	{
		StrandIterator it = sequence.NegativeByIndex(pos);
		for(int i = 0; i < edgeSize_ && it.Valid(); i++, it++)
		{
			StrandIterator end = sequence.PositiveByIndex(it.GetPosition());
			if(AdvanceForward(end, edgeSize_ - 1).Valid())
			{
				edge_->Insert(it.GetPosition());
			}
		}

		if(erase)
		{
			CalcBound();
		}
	}*/
}