//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "bifurcationstorage.h"
namespace SyntenyFinder
{
	const size_t BifurcationStorage::POSITIVE_BIT = 0;
	const size_t BifurcationStorage::NEGATIVE_BIT = 1;
	const BifurcationStorage::BifurcationId BifurcationStorage::NO_BIFURCATION = -1;	

	namespace
	{
		template<class T>
			bool PtrEqual(T * ptr1, T * ptr2)
			{
				return ptr1 != 0 && ptr2 != 0 && *ptr1 == *ptr2;
			}
	}
	
	bool BifurcationStorage::IteratorProxy::Valid() const
	{
		return (*ptr_).get_padding_int() != NO_BIFURCATION;
	}
	
	StrandIterator BifurcationStorage::IteratorProxy::operator*() const
	{
		return StrandIterator(*ptr_, direction_);
	}

	void BifurcationStorage::Cleanup()
	{
		for(size_t i = 0; i < toClear_.size(); i++)
		{
			toClear_[i].first->erase(toClear_[i].second);
		}

		toClear_.clear();
	}
	
	void BifurcationStorage::Clear()
	{
		maxId_ = 0;
		for(size_t strand = 0; strand < 2; strand++)
		{
			posBifurcation_[strand].clear();
			bifurcationPos_[strand].assign(maxId_ + 1, IteratorList());
		}
	}
	
	BifurcationStorage::BifurcationStorage(size_t maxId): maxId_(static_cast<BifurcationId>(maxId))
	{
		for(size_t strand = 0; strand < 2; strand++)
		{
			bifurcationPos_[strand].assign(maxId_ + 1, IteratorList());
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
				for(IteratorList::const_iterator it = bifurcationPos_[strand][bifId].begin(); it != bifurcationPos_[strand][bifId].end(); ++it)
				{
					if((*it).get_padding_int() != NO_BIFURCATION)
					{
						StrandIterator jt(*it, static_cast<DNASequence::Direction>(strand));
						size_t pos = sequence.GlobalIndex(jt);
						out << " {" << bifId << ", " << pos << ", ";
						CopyN(jt, k, std::ostream_iterator<char>(out));
						out << "}";
					}
				}
			}
			
			out << std::endl << strandName[strand] << ", pos:" ;
			for(IteratorMap::const_iterator it = posBifurcation_[strand].begin(); it != posBifurcation_[strand].end(); ++it)
			{
				StrandIterator jt(**it, static_cast<DNASequence::Direction>(strand));
				size_t pos = sequence.GlobalIndex(jt);
				out << " {" << pos << ", ";
				CopyN(jt, k, std::ostream_iterator<char>(out));
				out << ", " << (**it).get_padding_int() << "}";
			}			

			out << std::endl;
		}
	}
	
	void BifurcationStorage::AddPoint(DNASequence::StrandIterator it, size_t inBifId)
	{
		BifurcationId bifId = static_cast<BifurcationId>(inBifId);
		if(GetBifurcation(it) == NO_BIFURCATION && inBifId != NO_BIFURCATION)
		{
			size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
			it.SetInfoBit(strand, true);
			BaseIterator newIt = it.Base();
			newIt.get_padding_int() = bifId;
			IteratorPtr newPtr = bifurcationPos_[strand][bifId].insert(bifurcationPos_[strand][bifId].begin(), newIt);
			posBifurcation_[strand].insert(newPtr);
			assert(GetBifurcation(it) == bifId);
		}
	}

	
	BifurcationStorage::BifurcationId BifurcationStorage::ErasePointInternal(DNASequence::StrandIterator it, IteratorPtr & ret)
	{				
		BifurcationId bifId = NO_BIFURCATION;
		size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
		IteratorMap::iterator kt = LookUp(it);
		if(kt != posBifurcation_[strand].end())
		{
			ret = *kt;
			bifId = (**kt).get_padding_int();
			posBifurcation_[strand].erase(kt);
		}

		return bifId;
	}

	void BifurcationStorage::ErasePoint(DNASequence::StrandIterator it)
	{
		IteratorPtr buf;
		size_t bifId = ErasePointInternal(it, buf);
		if(bifId != NO_BIFURCATION)
		{				
			size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
			it.SetInfoBit(strand, false);
			(*buf).get_padding_int() = NO_BIFURCATION;
			toClear_.push_back(std::make_pair(&bifurcationPos_[strand][bifId], buf));
		}
	}
	
	size_t BifurcationStorage::GetBifurcation(DNASequence::StrandIterator it) const
	{
		size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;		
		IteratorMap::const_iterator kt = LookUp(it);
		return kt == posBifurcation_[strand].end() ? NO_BIFURCATION : (**kt).get_padding_int();
	}

	void BifurcationStorage::NotifyBefore(StrandIterator begin, StrandIterator end)
	{
		size_t pos = 0;
		nowInvalid_ = 0;
		invalid_.push_back(std::vector<BifurcationRecord>());
		for(StrandIterator it = begin; it != end; ++it, ++pos)
		{			
			IteratorPtr ret;
			BifurcationId bifId = ErasePointInternal(it, ret);
			if(bifId != NO_BIFURCATION)
			{
				invalid_.back().push_back(BifurcationRecord(pos, ret, bifId));
			}
		}		
	}

	void BifurcationStorage::NotifyAfter(StrandIterator begin, StrandIterator end)
	{
		size_t pos = 0;
		size_t record = 0;
		for(StrandIterator it = begin; it != end; ++it, ++pos)
		{
			if(record < invalid_[nowInvalid_].size() && invalid_[nowInvalid_][record].pos == pos)
			{
				size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
				BifurcationId bifId = invalid_[nowInvalid_][record].bifId;
				IteratorPtr jt = invalid_[nowInvalid_][record].ptrIt;				
				BaseIterator newIt = it.Base();
				newIt.get_padding_int() = bifId;
				*jt = newIt;
				posBifurcation_[strand].insert(jt);
				record++;
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
			for(ListVector::const_iterator it = bifurcationPos_[dir].begin(); it != bifurcationPos_[dir].end(); ++it)
			{
				for(IteratorList::const_iterator jt = it->begin(); jt != it->end(); ++jt)
				{
					StrandIterator begin(*jt, type);
					std::string body(begin, AdvanceForward(begin, k));
					dict[body] = it - bifurcationPos_[dir].begin();
				}				
			}
		}		
	}

	BifurcationStorage::IteratorMap::const_iterator BifurcationStorage::LookUp(StrandIterator it) const
	{
		size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
		if(it.GetInfoBit(strand))
		{
			IteratorList temp(1, it.Base());
			BifurcationId bifId = NO_BIFURCATION;
			IteratorPtr lookUp = temp.begin();
			return posBifurcation_[strand].find(lookUp);		
		}

		return posBifurcation_[strand].end();
	}

	BifurcationStorage::IteratorMap::iterator BifurcationStorage::LookUp(StrandIterator it)
	{
		size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
		if(it.GetInfoBit(strand))
		{
			IteratorList temp(1, it.Base());
			BifurcationId bifId = NO_BIFURCATION;
			IteratorPtr lookUp = temp.begin();
			return posBifurcation_[strand].find(lookUp);		
		}

		return posBifurcation_[strand].end();
	}
}