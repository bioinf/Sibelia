//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _VARIANT_H_
#define _VARIANT_H_

#include "../blockfinder.h"

namespace SyntenyFinder
{
	class Variant
	{
	public:	
		void ToString(std::string & buf) const;
		Variant(size_t refPos, size_t blockId, bool collinear, const std::string & refAllele, const std::string & altAllele);
		bool operator < (const Variant & toCompare) const;
	private:
		static const size_t UNKNOWN_BLOCK;
		size_t refPos_;
		size_t blockId_;
		std::string refAllele_;
		std::string altAllele_;
		bool collinear_;
	};
}

#endif