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
	dirs.push_back("Resources");
	//TODO: add system resource path
	return dirs;
}
