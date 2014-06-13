//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _UTIL_H_
#define _UTIL_H_

#include "blockfinder.h"
#include "outputgenerator.h"

void SignalHandler(int sig);
std::vector<std::pair<int, int> > FarStageFile();
std::vector<std::pair<int, int> > FineStageFile();
std::vector<std::pair<int, int> > LooseStageFile();
std::vector<std::pair<int, int> > ReadStageFile(const std::string & fileName);
void PutProgressChr(size_t progress, SyntenyFinder::BlockFinder::State state);

#endif
