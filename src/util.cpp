//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "util.h"

const std::string DELIMITER(80, '-');

ParameterSet::ParameterSet(int k, int maxBranchSize, int minPathLength):
	k(k), maxBranchSize(maxBranchSize), minPathLength(minPathLength)
{
}

std::vector<ParameterSet> ReadStageFile(const std::string & fileName)
{
	int count = 0;
	std::ifstream stageFile(fileName.c_str());
	if(!stageFile)
	{
		throw std::runtime_error("cannot open stage file");
	}

	if(!(stageFile >> count))
	{
		throw std::runtime_error("cannot read stage file");
	}

	if(count <= 0)
	{
		throw std::runtime_error("number of stages must be positive ");
	}

	std::vector<ParameterSet> ret(count);
	for(int i = 0; i < count; i++)
	{
		if(!(stageFile >> ret[i].k >> ret[i].maxBranchSize >> ret[i].minPathLength))
		{
			throw std::runtime_error("too few records in the stage file");
		}

		if(ret[i].k < 2)
		{
			throw std::runtime_error("vertex size in stage record must be at least 2");
		}

		if(ret[i].maxBranchSize < 0)
		{
			throw std::runtime_error("maximum branch size in stage record must be nonnegative");
		}
	}

	return ret;
}

std::vector<ParameterSet> LooseStageFile()
{
	ParameterSet stage[] = 
	{
		ParameterSet(30, 150, 0),
		ParameterSet(100, 1000, 0),
		ParameterSet(1000, 5000, 0),
		ParameterSet(5000, 15000, 0)		
	};

	return std::vector<ParameterSet>(stage, stage + sizeof(stage) / sizeof(stage[0]));
}

std::vector<ParameterSet> FineStageFile()
{
	ParameterSet stage[] = 
	{
		ParameterSet(30, 150, 0),
		ParameterSet(100, 500, 0),
		ParameterSet(500, 1500, 0)
	};

	return std::vector<ParameterSet>(stage, stage + sizeof(stage) / sizeof(stage[0]));
}

void PutProgressChr(size_t progress, SyntenyFinder::BlockFinder::State state)
{
	static size_t prev = 0;
	while(prev < progress)
	{
		prev++;
		std::cout << '.';
		std::cout.flush();
	}

	switch(state)
	{
	case SyntenyFinder::BlockFinder::start:
		prev = 0;
		std::cout << '[';
		break;
	case SyntenyFinder::BlockFinder::end:
		std::cout << "]" << std::endl;
		break;
	case SyntenyFinder::BlockFinder::run:
		;
	}
}


void SignalHandler(int sig)
{
	int entered = 0;
	if(entered++ == 0)
	{
		SyntenyFinder::TempFile::Cleanup();
	}

	exit(1);
}
