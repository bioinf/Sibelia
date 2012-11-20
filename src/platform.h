//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _PLATFORM_H_
#define _PLATFORM_H_

#include <string>
#include <vector>

namespace SyntenyFinder
{
	std::vector<std::string> GetResourceDirs();
	bool CreateDirectory(const std::string & path);
}

#endif