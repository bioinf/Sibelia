//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************
#include "rearrangement.h"

namespace SyntenyFinder
{
	std::ostream& operator << (std::ostream & out, const Reversal & reversal)
	{
		return out << "Reversal\t" << reversal.start_ << '\t' << reversal.end_;
	}

	std::ostream& operator << (std::ostream & out, const Translocation & translocation)
	{
		return out << "Translocation\t" << translocation.start_ << '\t' << translocation.end_ << '\t' << translocation.insert_;
	}
}