//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _PLATFORM_H_
#define _PLATFORM_H_

#include "common.h"

namespace SyntenyFinder
{
	std::vector<std::string> GetResourceDirs();
	void CreateOutDirectory(const std::string & path);

	class TempFile
	{
	public:
		TempFile();
		TempFile(const std::string & directory);
		TempFile(const TempFile & toCopy);
		~TempFile();
		void Rewind();
		void Write(const void * ptr, size_t size, size_t count);
		void Read(void * ptr, size_t size, size_t count);
		static void Cleanup();
	private:
		FILE * handle_;
		std::string path_;
		static std::map<std::string, FILE*> register_;
	};
}

#endif