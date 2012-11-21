//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "platform.h"

namespace SyntenyFinder
{
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

	#include <errno.h>
	#include <sys/types.h>
	#include <sys/stat.h>
	
	#ifdef _WIN32
		#include <direct.h>
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

	TempFile::TempFile()
	{
	}

	TempFile::TempFile(const std::string & directory): handle_(0)
	{
		srand(static_cast<unsigned int>(time(0)));
		for(size_t attempt = 0; attempt < 250; attempt++)
		{
			std::string fileName = "Sa";
			while(fileName.size() < L_tmpnam - 3)
			{
				fileName += 'a' + rand() % ('z' - 'a' + 1);
			}

			path_ = directory + "/" + fileName + ".tmp";
		#ifdef _WIN32
			struct __stat64 buf;
			int res = _stat64(path_.c_str(), &buf);			
		#else
			struct stat buf;
			int res = stat(path_.c_str(), &buf);
		#endif

			bool notExists = res == -1 && errno == ENOENT;
			if(notExists && (handle_ = fopen(path_.c_str(), "w+b")) != 0)
			{
				break;
			}
		}

		if(handle_ == 0)
		{
			throw std::runtime_error("Can't create a temporary file, see USAGE how to resolve this");
		}
	}

	TempFile::TempFile(const TempFile & toCopy): handle_(toCopy.handle_)
	{
	}

	TempFile::~TempFile()
	{
		if(handle_ != 0)
		{
			fclose(handle_);
			remove(path_.c_str());
		}
	}

	void TempFile::Rewind()
	{
		rewind(handle_);
	}

	void TempFile::Write(const void * ptr, size_t size, size_t count)
	{
		if(fwrite(ptr, size, count, handle_) != count)
		{
			throw std::runtime_error("Error while writing to a temporary file");
		}
	}

	void TempFile::Read(void * ptr, size_t size, size_t count)
	{
		if(fread(ptr, size, count, handle_) != count)
		{
			throw std::runtime_error("Error while reading from a temporary file");
		}
	}

}