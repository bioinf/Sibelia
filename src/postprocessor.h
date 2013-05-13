//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _POSTPROCESSOR_H_
#define _POSTPROCESSOR_H_

#include "fasta.h"
#include "blockinstance.h"

namespace SyntenyFinder
{

	class Postprocessor
	{
	public:
		static void GlueStripes(std::vector<BlockInstance> & block, const std::vector<FASTARecord> & chr);
	};

}

#endif