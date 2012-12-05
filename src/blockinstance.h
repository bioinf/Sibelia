//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _BLOCK_INSTANCE_H_
#define _BLOCK_INSTANCE_H_

#include "common.h"
#include "dnasequence.h"


namespace SyntenyFinder
{
	inline int Abs(int x)
	{
		return x > 0 ? x : -x;
	}

	class BlockInstance
	{
	public:
		BlockInstance() {}
		BlockInstance(int id, const FASTARecord * chr, size_t start, size_t end): id_(id), chr_(chr), start_(start), end_(end) {}
		int GetSignedBlockId() const
		{
			return id_;
		}

		DNASequence::Direction GetDirection() const
		{
			return id_ > 0 ? DNASequence::positive : DNASequence::negative;
		}

		int GetBlockId() const
		{
			return Abs(id_);
		}

		const FASTARecord& GetChrInstance() const
		{
			return *chr_;
		}

		size_t GetChrId() const
		{
			return GetChrInstance().GetId();
		}

		size_t GetStart() const
		{
			return start_;
		}

		size_t GetEnd() const
		{
			return end_;
		}

		size_t GetConventionalStart() const
		{
			if(GetDirection() == DNASequence::positive)
			{
				return start_ + 1;
			}

			return end_;
		}

		size_t GetConventionalEnd() const
		{
			if(GetDirection() == DNASequence::positive)
			{
				return end_;
			}

			return start_ + 1;
		}

		size_t GetLength() const
		{
			return end_ - start_;
		}

		bool operator < (const BlockInstance & toCompare) const
		{
			return std::make_pair(GetChrInstance().GetId(), GetStart()) < std::make_pair(toCompare.GetChrInstance().GetId(), toCompare.GetStart());
		}

	private:
		int id_;		
		size_t start_;
		size_t end_;
		const FASTARecord * chr_;
	};

	template<class F>
		bool CompareBlocks(const BlockInstance & a, const BlockInstance & b, F f)
		{
			return f(a) < f(b);
		}

	template<class F>
		bool EqualBlocks(const BlockInstance & a, const BlockInstance & b, F f)
		{
			return f(a) == f(b);
		}
	
	inline bool CompareBlocksNaturally(const BlockInstance & a, const BlockInstance & b)
	{
		return a < b;
	}

	typedef boost::function<bool(const BlockInstance&, const BlockInstance&)> BlockComparer;
	extern const BlockComparer compareById;
	extern const BlockComparer compareByChrId;
	extern const BlockComparer compareByStart;
}


#endif