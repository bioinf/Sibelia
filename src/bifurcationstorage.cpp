#include "bifurcationstorage.h"
namespace SyntenyBuilder
{
	const size_t BifurcationStorage::NO_BIFURCATION = -1;

	void BifurcationStorage::Clear()
	{
		for(size_t strand = 0; strand < 2; strand++)
		{
			posBifurcation[strand].clear();
			bifurcationPos[strand].clear();
		}
	}

	void BifurcationStorage::Dump(std::ostream & out) const
	{
		std::string strandName[] = {"Positive:", "Negative:"};
		for(size_t strand = 0; strand < 2; strand++)
		{
			out << strandName[strand];
			for(PosBifurcation::const_iterator it = posBifurcation[strand].begin();
				it != posBifurcation[strand].end(); ++it)
			{
				out << " {" << (*it)->second << ", " << (*it)->first << "}";
			}

			out << std::endl;
		}
	}

	void BifurcationStorage::AddPoint(DNASequence::StrandIterator it, size_t bifId)
	{
		size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
		BifMapIterator place = bifurcationPos[strand].insert(std::make_pair(bifId, it.GetPosition()));
		posBifurcation[strand].insert(place);
	}

	void BifurcationStorage::ErasePoint(DNASequence::StrandIterator it)
	{
		size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
		BifMapIterator jt = temp.insert(std::make_pair(-1, it.GetPosition()));
		PosBifurcation::iterator kt = posBifurcation[strand].find(jt);
		if(kt != posBifurcation[strand].end())
		{
			bifurcationPos[strand].erase(*kt);
			posBifurcation[strand].erase(kt);
		}

		temp.clear();
	}

	size_t BifurcationStorage::GetBifurcation(DNASequence::StrandIterator it) const
	{
		size_t strand = it.GetDirection() == DNASequence::positive ? 0 : 1;
		BifMapIterator jt = temp.insert(std::make_pair(-1, it.GetPosition()));
		PosBifurcation::const_iterator kt = posBifurcation[strand].find(jt);
		temp.clear();
		return kt == posBifurcation[strand].end() ? NO_BIFURCATION : (*kt)->first;
	}
}