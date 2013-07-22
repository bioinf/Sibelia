//****************************************************************************
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
		dirs.push_back("/usr/share/sibelia/resources");
	#endif
		return dirs;
	}

	#include <errno.h>
	#include <sys/types.h>
	#include <sys/stat.h>
	
	#ifdef _WIN32
		#include <direct.h>
	#endif

	void CreateOutDirectory(const std::string & path)
	{
		int result = 0;
	#ifdef _WIN32
		result = _mkdir(path.c_str());
	#else
		result = mkdir(path.c_str(), 0755);
	#endif
		if (result != 0 && errno != EEXIST)
		{
			throw std::runtime_error(("Cannot create dir " + path).c_str());
		}
	}

	std::map<std::string, FILE*> TempFile::register_;

	TempFile::TempFile()
	{
	}

	TempFile::TempFile(const std::string & directory): handle_(0)
	{
		for(size_t attempt = 0; attempt < 1000; attempt++)
		{
			std::string fileName = "Sib_";
			for(;fileName.size() < 16;)
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
				register_[path_] = handle_;
				break;
			}
		}

		if(handle_ == 0)
		{
			throw std::runtime_error("Can't create a temporary file, see USAGE how to resolve this");
		}
	}

	void TempFile::Cleanup()
	{
		for(std::map<std::string, FILE*>::iterator it = register_.begin(); it != register_.end(); ++it)
		{
			if(it->second != 0)
			{
				fclose(it->second);
				remove(it->first.c_str());
			}
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
			register_.erase(path_);
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