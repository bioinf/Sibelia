#ifndef _FASTA_READER_H_
#define _FASTA_READER_H_

#include "lib/kseq.h"
#include "common.h"

namespace SyntenyBuilder
{
	//This class provides functionality for reading FASTA files.
	//Currently this is an OO wrapper around tiny library kseq.h.
	//It looks little ugly, maybe later it will be replaced with
	//some other library or written by myself.

	class FASTAReader
	{
	public:
		bool IsOk() const
		{
			return fileHandler_ != 0 && !feof(fileHandler_) && !ferror(fileHandler_);
		}

		explicit FASTAReader(const std::string & fileName): fileHandler_(fopen(fileName.c_str(), "r"))
		{

		}

		~FASTAReader()
		{
			if(fileHandler_ != 0)
			{
				fclose(fileHandler_);
			}
		}

		void GetSequence(std::string & buf);
		
	private:
		DISALLOW_COPY_AND_ASSIGN(FASTAReader);
		FILE * fileHandler_;
	};
}

#endif 
