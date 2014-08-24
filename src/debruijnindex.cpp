//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "debruijnindex.h"

namespace SyntenyFinder
{
	const char DeBruijnIndex::SEPARATION_CHAR = '#';

	size_t DeBruijnIndex::GetStrand(FastaRecord::Direction dir)
	{
		return dir == FastaRecord::positive ? 0 : 1;
	}

	DeBruijnIndex::Location::Location(size_t chrId, size_t index):
		chrId_(static_cast<uint32_t>(chrId)), index_(static_cast<uint32_t>(index))
	{
	}

	size_t DeBruijnIndex::Location::GetIndex() const
	{
		return index_;
	}
	
	size_t DeBruijnIndex::Location::GetChromosomeId() const
	{
		return chrId_;
	}

	size_t DeBruijnIndex::BifurcationData::GetBifurcationId() const
	{
		return bifId_;
	}

	size_t DeBruijnIndex::BifurcationData::GetProjection() const
	{
		return projection_;
	}

	DeBruijnIndex::BifurcationData::BifurcationData()
	{
	}

	DeBruijnIndex::BifurcationData::BifurcationData(size_t pos, size_t bifId, size_t projection, char inMark, char outMark):
		pos_(static_cast<uint32_t>(pos)), bifId_(static_cast<uint32_t>(bifId)), projection_(static_cast<uint32_t>(projection)),
		inMark_(inMark), outMark_(outMark)
	{
	}

	char DeBruijnIndex::BifurcationData::GetInMark() const
	{
		return inMark_;
	}

	char DeBruijnIndex::BifurcationData::GetOutMark() const
	{
		return outMark_;
	}

	bool DeBruijnIndex::BifurcationData::IsValid() const
	{
		return valid_;
	}

	void DeBruijnIndex::BifurcationData::Invalidate()
	{
		valid_ = false;
	}

	size_t DeBruijnIndex::BifurcationData::GetPosition() const
	{
		return pos_;
	}

	DeBruijnIndex::DeBruijnIndex(const std::vector<ChrBifVector> & bifurcation, const std::vector<std::string> & record, size_t k, size_t bifurcationNumber)
	{
		revCompDictionary_.resize(bifurcationNumber);
		bifurcationData_.resize(bifurcation[0].size());
		for(size_t chr = 0; chr < bifurcation[0].size(); chr++)
		{
			bifurcationData_[chr].reserve(bifurcation[0][chr].size());
			for(size_t i = 0; i < bifurcation[0][chr].size(); i++)
			{
				size_t j = bifurcation[0][chr].size() - i - 1;
				size_t bifId = bifurcation[0][chr][i].GetId();
				size_t revBifId = bifurcation[1][chr][j].GetId();
				size_t bifPos = bifurcation[0][chr][i].GetPostion();
				revCompDictionary_[bifId] = revBifId;
				revCompDictionary_[revBifId] = bifId;
				char inMark = bifPos == 0 ? SEPARATION_CHAR : record[chr][bifPos - 1];
				char outMark = bifPos + k >= record[chr].size() ? SEPARATION_CHAR : record[chr][bifPos + 1];
				bifurcationData_[chr].push_back(BifurcationData(bifPos, bifId, bifurcation[0][chr][i].GetProjection(), inMark, outMark));
				bifurcationPlace_[bifId].push_back(Location(chr, i));
			}
		}
	}

	//size_t DeBruijnIndex::GetEdgesOfVertex(size_t bifId, std::vector<Edge> & e) const
	//{
	//	e.clear();
	//	for(size_t i = 0; i < bifurcationPosition_[bifId].size(); i++)
	//	{
	//		Location location = bifurcationPosition_[bifId][i];
	//		e.push_back(GetEdgeAtPosition(location.GetChromosomeId(), location.GetPosition(), location.GetDirection()));
	//	}

	//	return e.size();
	//}
}
