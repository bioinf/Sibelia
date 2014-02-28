//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _UTIL_H_
#define _UTIL_H_

#include "blockfinder.h"
#include "outputgenerator.h"

struct ParameterSet
{
	int k;
	int maxBranchSize;
	int minPathLength;
	ParameterSet() {}
	ParameterSet(int k, int maxBranchSize, int minPathLength);
};

void SignalHandler(int sig);
std::vector<ParameterSet> FineStageFile();
std::vector<ParameterSet> LooseStageFile();
std::vector<ParameterSet> ReadStageFile(const std::string & fileName);
void PutProgressChr(size_t progress, SyntenyFinder::BlockFinder::State state);

#endif 