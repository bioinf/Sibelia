//****************************************************************************
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
		void Reverse();
		int GetSignedBlockId() const;
		DNASequence::Direction GetDirection() const;
		int GetBlockId() const;
		const FASTARecord& GetChrInstance() const;
		size_t GetChrId() const;
		size_t GetStart() const;
		size_t GetEnd() const;
		size_t GetLength() const;
		size_t GetConventionalStart() const;
		size_t GetConventionalEnd() const;
		std::pair<size_t, size_t> CalculateOverlap(const BlockInstance & instance) const;		
		bool operator < (const BlockInstance & toCompare) const;
		bool operator == (const BlockInstance & toCompare) const;
		bool operator != (const BlockInstance & toCompare) const;
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