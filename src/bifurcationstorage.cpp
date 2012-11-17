//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "bifurcationstorage.h"
namespace SyntenyFinder
{
	const BifurcationStorage::BifurcationId BifurcationStorage::NO_BIFURCATION = -1;
	
	bool BifurcationStorage::IteratorProxy::Valid() const
	{
		return ptr_ != 0;
	}

	StrandIterator BifurcationStorage::IteratorProxy::GetIterator() const
	{
		return StrandIterator(*ptr_, direction_);
	}

	void BifurcationStorage::Clear()
	{
		maxId_ = 0;
		for(size_t strand = 0; strand < 2; strand++)
		{
			posBifurcation_[strand].clear();
			bifurcationPos_[strand].assign(maxId_ + 1, IteratorVector());
		}
	}

	BifurcationStorage::BifurcationStorage(size_t maxId): maxId_(static_cast<BifurcationId>(maxId))
	{
		for(size_t strand = 0; strand < 2; strand++)
		{
			bifurcationPos_[strand].assign(maxId_ + 1, IteratorVector());
		}
	}

	size_t BifurcationStorage::GetMaxId() const
	{
		return maxId_;
	}

	size_t BifurcationStorage::TotalElements() const
	{
		return bifurcationPos_[0].size() + bifurcationPos_[1].size();
	}

	size_t BifurcationStorage::CountBifurcations(size_t inBifId) const
	{
		BifurcationId bifId = static_cast<BifurcationId>(inBifId);
		return bifurcationPos_[0][bifId].size() + bifurcationPos_[1][bifId].size();
	}

	void BifurcationStorage::Dump(const DNASequence & sequence, size_t k, std::ostream & out) const
	{
		std::string strandName[] = {"Positive", "Negative"};
		StrandIterator start[] = {sequence.PositiveBegin(0), sequence.NegativeBegin(0)};
		for(size_t strand = 0; strand < 2; strand++)
		{
			out << strandName[strand] << ", bif:" ;
			for(size_t bifId = 0; bifId < bifurcationPos_[strand].size(); bifId++)
			{
				for(size_t instance = 0; instance < bifurcationPos_[strand][bifId].size(); instance++)
				{
					StrandIterator jt(*bifurcationPos_[strand][bifId][instance], static_cast<DNASequence::Direction>(strand));
					size_t pos = sequence.GlobalIndex(jt);
					out << " {" << bifId << ", " << pos << ", ";
					CopyN(jt, k, std::ostream_iterator<char>(out));
					out << "}";
				}
			}
			
			out << std::endl << strandName[strand] << ", pos:" ;
			for(IteratorMap::const_iterator it = posBifurcation_[strand].begin(); it != posBifurcation_[strand].end(); ++it)
			{
				StrandIterator jt(*it->first, static_cast<DNASequence::Direction>(strand));
				size_t pos = sequence.GlobalIndex(jt);
				out << " {" << pos << ", ";
				CopyN(jt, k, std::ostream_iterator<char>(out));
				out << ", " << it->second << "}";
			}			

			out << std::endl;
		}
	}
	
	void BifurcationStorage::AddPoint(DNASequence::StrandIterator it, size_t inBifId)
	{
		BifurcationId bifId = static_cast<BifurcationId>(inBifId);
		if(GetBifurcation(it) == NO_BIFURCATION)
		{
			size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
			BaseIterator kt = it.Base();
			IteratorPtr newPtr = IteratorPtr(new BaseIterator(kt));
			bifurcationPos_[strand][bifId].push_back(newPtr);
			posBifurcation_[strand].insert(std::make_pair(newPtr, bifId));
			assert(GetBifurcation(it) == bifId);
		}
	}

	
	void BifurcationStorage::ErasePoint(DNASequence::StrandIterator it)
	{
		size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
		IteratorPtr lookUp(new BaseIterator(it.Base()));
		IteratorMap::iterator kt = posBifurcation_[strand].find(lookUp);
		if(kt != posBifurcation_[strand].end())
		{
			IteratorPtr ptr = kt->first;
			posBifurcation_[strand].erase(kt);
			ptr.reset();
		}
	}
	
	size_t BifurcationStorage::GetBifurcation(DNASequence::StrandIterator it) const
	{
		size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
		IteratorPtr lookUp(new BaseIterator(it.Base()));
		IteratorMap::const_iterator kt = posBifurcation_[strand].find(lookUp);
		return kt == posBifurcation_[strand].end() ? NO_BIFURCATION : kt->second;
	}

	void BifurcationStorage::NotifyBefore(StrandIterator begin, StrandIterator end)
	{
		size_t pos = 0;
		nowInvalid_ = 0;
		invalid_.push_back(std::vector<BifurcationRecord>());
		for(StrandIterator it = begin; it != end; ++it, ++pos)
		{			
			BifurcationId bifId = static_cast<BifurcationId>(GetBifurcation(it));
			if(bifId != NO_BIFURCATION)
			{
				invalid_.back().push_back(BifurcationRecord(pos, bifId));
				ErasePoint(it);
			}
		}		

	}

	void BifurcationStorage::NotifyAfter(StrandIterator begin, StrandIterator end)
	{
		size_t pos = 0;
		size_t record = 0;
		for(StrandIterator it = begin; it != end; ++it, ++pos)
		{
			if(record < invalid_[nowInvalid_].size() && invalid_[nowInvalid_][record].first == pos)
			{
				AddPoint(it, invalid_[nowInvalid_][record++].second);
			}
		}

		if(++nowInvalid_ == invalid_.size())
		{
			invalid_.clear();
		}
	}

	void BifurcationStorage::FormDictionary(boost::unordered_map<std::string, size_t> & dict, size_t k) const
	{
		dict.clear();
		for(size_t dir = 0; dir < 2; dir++)
		{
			DNASequence::Direction type = static_cast<DNASequence::Direction>(dir);
			for(BifurcationStore::const_iterator it = bifurcationPos_[dir].begin(); it != bifurcationPos_[dir].end(); ++it)
			{
				for(IteratorVector::const_iterator jt = it->begin(); jt != it->end(); ++jt)
				{
					StrandIterator begin(**jt, type);
					std::string body(begin, AdvanceForward(begin, k));
					dict[body] = it - bifurcationPos_[dir].begin();
				}				
			}
		}		
	}
}