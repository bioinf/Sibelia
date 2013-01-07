//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "util.h"

const std::string DELIMITER(80, '-');

std::vector<std::pair<int, int> > ReadStageFile(const std::string & fileName)
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

	if(count < 0)
	{
		throw std::runtime_error("number of stages must be nonnegative");
	}

	std::vector<std::pair<int, int> > ret(count);
	for(int i = 0; i < count; i++)
	{
		if(!(stageFile >> ret[i].first >> ret[i].second))
		{
			throw std::runtime_error("too few records in the stage file");
		}

		if(ret[i].first < 2)
		{
			throw std::runtime_error("vertex size in stage record must be at least 2");
		}

		if(ret[i].second < 0)
		{
			throw std::runtime_error("minimum branch size in stage record must be nonnegative");
		}
	}

	return ret;
}

std::vector<std::pair<int, int> > LooseStageFile()
{
	std::pair<int, int> stage[] = 
	{
		std::make_pair(30, 150),
		std::make_pair(100, 1000),
		std::make_pair(1000, 5000),		
		std::make_pair(5000, 15000)
	};

	return std::vector<std::pair<int, int> >(stage, stage + sizeof(stage) / sizeof(stage[0]));
}

std::vector<std::pair<int, int> > FineStageFile()
{
	std::pair<int, int> stage[] = 
	{
		std::make_pair(30, 150),
		std::make_pair(100, 1000),
		std::make_pair(1000, 2500),		
	};

	return std::vector<std::pair<int, int> >(stage, stage + sizeof(stage) / sizeof(stage[0]));
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
