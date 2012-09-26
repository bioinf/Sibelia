#include "bifurcationstorage.h"
namespace SyntenyBuilder
{
	const BifurcationStorage::BifurcationId BifurcationStorage::NO_BIFURCATION = -1;
	
	void BifurcationStorage::Clear()
	{
		maxId_ = 0;
		for(size_t strand = 0; strand < 2; strand++)
		{
			posBifurcation_[strand].clear();
			bifurcationPos_[strand].clear();
		}
	}

	BifurcationStorage::BifurcationStorage(): maxId_(0),
		bifurcationPos_(2), posBifurcation_(2)
	{
	}

	size_t BifurcationStorage::GetMaxId() const
	{
		return maxId_ + 1;
	}

	size_t BifurcationStorage::CountBifurcations(size_t inBifId) const
	{
		BifurcationId bifId = static_cast<BifurcationId>(inBifId);
		return bifurcationPos_[0].count(bifId) + bifurcationPos_[1].count(bifId);
	}

	void BifurcationStorage::Dump(const DNASequence & sequence, size_t k, std::ostream & out) const
	{		
		std::string strandName[] = {"Positive", "Negative"};
		StrandIterator start[] = {sequence.PositiveBegin(), sequence.NegativeBegin()};
		for(size_t strand = 0; strand < 2; strand++)
		{
			out << strandName[strand] << ", bif:" ;
			for(CBifMapIterator it = bifurcationPos_[strand].begin();
				it != bifurcationPos_[strand].end(); ++it)
			{
				StrandIterator jt(it->second, static_cast<DNASequence::Direction>(strand));
				size_t pos = 0; //std::distance(start[strand], jt);
				out << " {" << it->first << ", " << pos << ", ";
				CopyN(jt, k, std::ostream_iterator<char>(out));
				out << "}";
			}
			
			out << std::endl << strandName[strand] << ", pos:" ;
			for(PosBifurcation::const_iterator it = posBifurcation_[strand].begin();
				it != posBifurcation_[strand].end(); ++it)
			{
				StrandIterator jt((*it)->second, static_cast<DNASequence::Direction>(strand));
				size_t pos = 0; //std::distance(start[strand], jt);
				out << " {" << pos << ", ";
				CopyN(jt, k, std::ostream_iterator<char>(out));
				out << ", " << (*it)->first << "}";
			}			

			out << std::endl;
		}
	}

	void BifurcationStorage::AddPoint(DNASequence::StrandIterator it, size_t inBifId)
	{
		BifurcationId bifId = static_cast<BifurcationId>(inBifId);
		if(GetBifurcation(it) == NO_BIFURCATION)
		{
			maxId_ = std::max(maxId_, bifId);
			size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
			BaseIterator kt = it.Base();
			BifMapIterator place = bifurcationPos_[strand].insert(std::make_pair(bifId, it.Base()));			
			posBifurcation_[strand].insert(place);
			assert(GetBifurcation(it) == bifId);
		}
	}

	void BifurcationStorage::ErasePoint(DNASequence::StrandIterator it)
	{
		size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
		BifMapIterator jt = temp_.insert(std::make_pair(-1, it.Base()));
		PosBifurcation::iterator kt = posBifurcation_[strand].find(jt);
		if(kt != posBifurcation_[strand].end())
		{
			bifurcationPos_[strand].erase(*kt);
			posBifurcation_[strand].erase(kt);
		}

		temp_.clear();
	}

	size_t BifurcationStorage::GetBifurcation(DNASequence::StrandIterator it) const
	{
		size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
		BifMapIterator jt = temp_.insert(std::make_pair(-1, it.Base()));
		PosBifurcation::const_iterator kt = posBifurcation_[strand].find(jt);
		temp_.clear();
		return kt == posBifurcation_[strand].end() ? NO_BIFURCATION : (*kt)->first;
	}
}