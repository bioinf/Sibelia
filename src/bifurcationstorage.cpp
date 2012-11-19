//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "bifurcationstorage.h"
namespace SyntenyFinder
{
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
		return (*adjList_)[index_] != 0;
	}

	StrandIterator BifurcationStorage::IteratorProxy::operator*() const
	{
		return StrandIterator(*(*adjList_)[index_], direction_);
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

	BifurcationStorage::BifurcationStorage(size_t maxId): maxId_(static_cast<BifurcationId>(maxId)), empty_(0)
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
					if(bifurcationPos_[strand][bifId][instance] != 0)
					{
						StrandIterator jt(*bifurcationPos_[strand][bifId][instance], static_cast<DNASequence::Direction>(strand));
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
				out << ", " << (*it)->get_padding_int() << "}";
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
			BaseIterator kt = it.Base();
			IteratorPtr newPtr = IteratorPtr(new BaseIterator(kt));		
			newPtr->get_padding_int() = bifId;
			IteratorVector & line = bifurcationPos_[strand][bifId];
			typedef std::equal_to<IteratorPtr> EqualOp;
			IteratorVector::iterator jt = std::find_if(line.begin(), line.end(), boost::bind(&EqualOp::operator(), boost::cref(EqualOp()), IteratorPtr(0), _1));
			if(jt == line.end())
			{
				bifurcationPos_[strand][bifId].push_back(newPtr);
			}
			else
			{
				empty_--;
				*jt = newPtr;
			}

			posBifurcation_[strand].insert(newPtr);
			assert(GetBifurcation(it) == bifId);
		}
	}

	
	BifurcationStorage::BifurcationId BifurcationStorage::ErasePointInternal(DNASequence::StrandIterator it, IteratorVector::iterator & ret)
	{		
		BifurcationId bifId = NO_BIFURCATION;
		BaseIterator base = it.Base();
		IteratorWeakPtr lookUp = &base;
		size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
		IteratorMap::iterator kt = posBifurcation_[strand].find(lookUp);		
		if(kt != posBifurcation_[strand].end())
		{
			bifId = (*kt)->get_padding_int();
			IteratorVector & line = bifurcationPos_[strand][bifId];
			ret = std::find_if(line.begin(), line.end(), boost::bind(PtrEqual<BaseIterator>, *kt, _1));
			posBifurcation_[strand].erase(kt);
			empty_++;
			delete *ret;
			*ret = 0;
		}

		return bifId;
	}

	void BifurcationStorage::ErasePoint(DNASequence::StrandIterator it)
	{
		IteratorVector::iterator buf;
		ErasePointInternal(it, buf);
	}
	
	size_t BifurcationStorage::GetBifurcation(DNASequence::StrandIterator it) const
	{
		size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
		BaseIterator base = it.Base();
		IteratorWeakPtr lookUp = &base;
		IteratorMap::const_iterator kt = posBifurcation_[strand].find(lookUp);
		return kt == posBifurcation_[strand].end() ? NO_BIFURCATION : (*kt)->get_padding_int();
	}

	void BifurcationStorage::NotifyBefore(StrandIterator begin, StrandIterator end)
	{
		size_t pos = 0;
		nowInvalid_ = 0;
		invalid_.push_back(std::vector<BifurcationRecord>());
		for(StrandIterator it = begin; it != end; ++it, ++pos)
		{			
			IteratorVector::iterator ret;
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
				BifurcationId bifId = invalid_[nowInvalid_][record].bifId;
				IteratorVector::iterator jt = invalid_[nowInvalid_][record].ptrIt;
				size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
				IteratorPtr ptr(new BaseIterator(it.Base()));
				ptr->get_padding_int() = bifId;
				*jt = ptr;
				posBifurcation_[strand].insert(ptr);
				record++;
			}
		}

		if(++nowInvalid_ == invalid_.size())
		{
			invalid_.clear();
		}
	}

	BifurcationStorage::~BifurcationStorage()
	{
		for(size_t strand = 0; strand < 2; strand++)
		{
			for(BifurcationStore::iterator it = bifurcationPos_[strand].begin(); it != bifurcationPos_[strand].end(); ++it)
			{
				for(IteratorVector::iterator jt = it->begin(); jt != it->end(); ++jt)
				{
					delete *jt;
				}
			}
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