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

	bool DeBruijnIndex::BifurcationData::operator < (const BifurcationData & data) const
	{
		return pos_ < data.pos_;
	}

	DeBruijnIndex::BifurcationData::BifurcationData()
	{
	}

	DeBruijnIndex::BifurcationData::BifurcationData(size_t pos): pos_(static_cast<uint32_t>(pos))
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

	DeBruijnIndex::DeBruijnIndex(const std::vector<ChrBifVector> & bifurcation, const std::vector<std::string> & record,
		size_t k, size_t bifurcationNumber, const std::vector<size_t> & originalChrSize): k_(k), originalChrSize_(originalChrSize)
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

	size_t DeBruijnIndex::GetBifurcationsNumber() const
	{
		return revCompDictionary_.size();
	}
	
	size_t DeBruijnIndex::CountInstances(size_t bifId) const
	{
		size_t revBifId = revCompDictionary_[bifId];
		return bifurcationPlace_[bifId].size() + bifurcationPlace_[revBifId].size();
	}

	size_t DeBruijnIndex::GetBifurcationInstances(size_t bifId, std::vector<BifurcationIterator> & ret) const
	{
		ret.clear();
		for(size_t i = 0; i < bifurcationPlace_[bifId].size(); i++)
		{
			Location place = bifurcationPlace_[bifId][i];
			ret.push_back(BifurcationIterator(this, place.GetChromosomeId(), place.GetIndex(), FastaRecord::positive));
		}

		size_t revBifId = revCompDictionary_[bifId];
		for(size_t i = 0; i < bifurcationPlace_[revBifId].size(); i++)
		{
			Location place = bifurcationPlace_[revBifId][i];
			ret.push_back(BifurcationIterator(this, place.GetChromosomeId(), place.GetIndex(), FastaRecord::negative));
		}

		return ret.size();
	}

	DeBruijnIndex::BifurcationIterator DeBruijnIndex::Begin(size_t chr, FastaRecord::Direction dir) const
	{
		return BifurcationIterator(this, chr, 0, dir);
	}

	DeBruijnIndex::BifurcationIterator DeBruijnIndex::End(size_t chr, FastaRecord::Direction dir) const
	{
		return BifurcationIterator(this, chr, bifurcationData_[chr].size(), dir);
	}
	
	DeBruijnIndex::BifurcationIterator::BifurcationIterator(): parent_(0)
	{
	}

	DeBruijnIndex::BifurcationIterator::BifurcationIterator(const DeBruijnIndex * parent, size_t chrId, size_t index, FastaRecord::Direction dir):
		parent_(parent), chrId_(chrId), index_(index), dir_(dir)	
	{
	}

	DeBruijnIndex::BifurcationData DeBruijnIndex::BifurcationIterator::MyData() const
	{
		if(dir_ == FastaRecord::positive)
		{
			return *(parent_->bifurcationData_[chrId_].begin() + index_);
		}

		return *(parent_->bifurcationData_[chrId_].rbegin() + index_);
	}

	bool DeBruijnIndex::BifurcationIterator::IsValid() const
	{
		return MyData().IsValid();
	}
	
	bool DeBruijnIndex::BifurcationIterator::AtEnd() const
	{
		return index_ == parent_->bifurcationData_[chrId_].size();
	}
	
	char DeBruijnIndex::BifurcationIterator::GetOutMark() const
	{
		if(dir_ == FastaRecord::positive)
		{
			return MyData().GetOutMark();
		}

		return FastaRecord::Translate(MyData().GetInMark());
	}
		
	size_t DeBruijnIndex::BifurcationIterator::GetPosition() const
	{
		size_t pos = MyData().GetPosition();
		if(dir_ == FastaRecord::positive)
		{
			return pos;
		}

		return parent_->bifurcationData_[chrId_].back().GetPosition() - pos;
	}

	size_t DeBruijnIndex::BifurcationIterator::GetChromosomeId() const
	{
		return chrId_;
	}

	size_t DeBruijnIndex::BifurcationIterator::GetBifurcationId() const
	{
		size_t bifId = MyData().GetBifurcationId();
		if(dir_ == FastaRecord::positive)
		{
			return bifId;
		}

		return parent_->revCompDictionary_[bifId];
	}
	
	size_t DeBruijnIndex::BifurcationIterator::GetPositivePosition() const
	{
		size_t pos = MyData().GetPosition();
		if(dir_ == FastaRecord::positive)
		{
			return pos;
		}

		return pos + parent_->k_ - 1;
	}

	size_t DeBruijnIndex::BifurcationIterator::GetProjection() const
	{
		size_t proj = MyData().GetProjection();
		if(dir_ == FastaRecord::positive)
		{
			return proj;
		}

		if(index_ == 0)
		{
			return parent_->originalChrSize_[chrId_];
		}
		
		size_t pos = GetPositivePosition();
		size_t index = parent_->bifurcationData_[chrId_].size() - index_ - 1;		
		std::vector<BifurcationData>::const_iterator it = std::lower_bound(
			parent_->bifurcationData_[chrId_].begin() + index + 1,
			parent_->bifurcationData_[chrId_].end(),
			BifurcationData(pos));

		if(pos - (it - 1)->GetPosition() < it->GetPosition() - pos)
		{
			--it;
		}

		double coeff = double(it->GetPosition() - MyData().GetPosition()) / (it->GetProjection() - MyData().GetProjection());
		return MyData().GetProjection() + static_cast<size_t>(parent_->k_ * coeff);
	}

	size_t DeBruijnIndex::BifurcationIterator::GetPositiveEndingPosition() const
	{
		size_t pos = MyData().GetPosition();
		if(dir_ == FastaRecord::positive)
		{
			return pos + parent_->k_ - 1;
		}

		return pos;
	}

	DeBruijnIndex::BifurcationIterator DeBruijnIndex::BifurcationIterator::operator + (size_t shift) const
	{
		return BifurcationIterator(parent_, chrId_, index_ + shift, dir_);
	}

	DeBruijnIndex::BifurcationIterator& DeBruijnIndex::BifurcationIterator::operator++()
	{
		index_++;
		return *this;
	}

	DeBruijnIndex::BifurcationIterator DeBruijnIndex::BifurcationIterator::operator++(int)
	{
		BifurcationIterator ret(parent_, chrId_, index_++, dir_);
		return ret;
	}
			
	bool DeBruijnIndex::BifurcationIterator::operator == (const BifurcationIterator & it) const
	{
		return parent_ == it.parent_ && chrId_ == it.chrId_ && dir_ == it.dir_ && index_ == it.index_;
	}
	
	bool DeBruijnIndex::BifurcationIterator::operator != (const BifurcationIterator & it) const
	{
		return !(*this == it);
	}
	
}
