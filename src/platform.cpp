//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "platform.h"

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

bool CreateDirectory(const std::string & path)
{
	int result = 0;
#ifdef _WIN32
	result = _mkdir(path.c_str());
#else
	result = mkdir(path.c_str(), 0755);
#endif
	return result == 0;
}
