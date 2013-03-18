//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _REARRANGEMENT_H_
#define _REARRANGEMENT_H_

#include "../common.h"

namespace SyntenyFinder
{
	class Reversal
	{
	public:
		Reversal() {}
		Reversal(size_t start, size_t end): start_(start), end_(end) {}
	private:
		size_t start_;
		size_t end_;
		friend std::ostream& operator << (std::ostream & out, const Reversal & reversal);
	};

	class Translocation
	{
	public:
		Translocation() {}
		Translocation(size_t start, size_t end, size_t insert): start_(start), end_(end), insert_(insert) {}
	private:
		size_t start_;
		size_t end_;
		size_t insert_;
		friend std::ostream& operator << (std::ostream & out, const Translocation & translocation);
	};

	std::ostream& operator << (std::ostream & out, const Reversal & reversal);
	std::ostream& operator << (std::ostream & out, const Translocation & translocation);
}

#endif