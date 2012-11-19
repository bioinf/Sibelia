//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "platform.h"

std::vector<std::string> GetResourceDirs()
{
	std::vector<std::string> dirs;
	//relative path to resources
	dirs.push_back("resources");
#ifdef __gnu_linux__
	//TODO: change with real installation dir
	dirs.push_back("/usr/share/sibelia");
#endif
	return dirs;
}
