#include "bifurcationstorage.h"
namespace SyntenyBuilder
{
	const size_t BifurcationStorage::NO_BIFURCATION = -1;
	
	void BifurcationStorage::Clear()
	{
		for(size_t strand = 0; strand < 2; strand++)
		{
			posBifurcation_[strand].clear();
			bifurcationPos_[strand].clear();
		}
	}

	size_t BifurcationStorage::GetMaxId() const
	{
		return maxId_ + 1;
	}

	size_t BifurcationStorage::CountBifurcations(size_t bifId) const
	{
		return bifurcationPos_[0].count(bifId) + bifurcationPos_[1].count(bifId);
	}

	void BifurcationStorage::Dump(std::ostream & out) const
	{/*
		std::string strandName[] = {"Positive", "Negative"};
		for(size_t strand = 0; strand < 2; strand++)
		{
			out << strandName[strand] << ", pos:" ;
			for(posBifurcation_::const_iterator it = posBifurcation_[strand].begin();
				it != posBifurcation_[strand].end(); ++it)
			{
				out << " {" << (*it)->second << ", " << (*it)->first << "}";
			}

			out << std::endl;

			out << strandName[strand] << ", bif:" ;
			for(CBifMapIterator it = bifurcationPos_[strand].begin();
				it != bifurcationPos_[strand].end(); ++it)
			{
				out << " {" << it->first << ", " << it->second << "}";
			}

			out << std::endl;
		}*/
	}

	void BifurcationStorage::AddPoint(DNASequence::StrandIterator it, size_t bifId)
	{
		if(GetBifurcation(it) == NO_BIFURCATION)
		{
			maxId_ = std::max(maxId_, bifId);
			size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
			BifMapIterator place = bifurcationPos_[strand].insert(std::make_pair(bifId, it));			
			posBifurcation_[strand].insert(place);
			assert(GetBifurcation(it) == bifId);
		}
	}

	void BifurcationStorage::ErasePoint(DNASequence::StrandIterator it)
	{
		size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
		BifMapIterator jt = temp_.insert(std::make_pair(-1, it));
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
		BifMapIterator jt = temp_.insert(std::make_pair(-1, it));
		PosBifurcation::const_iterator kt = posBifurcation_[strand].find(jt);
		temp_.clear();
		return kt == posBifurcation_[strand].end() ? NO_BIFURCATION : (*kt)->first;
	}
}