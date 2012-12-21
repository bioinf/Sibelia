//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _FASTA_READER_H_
#define _FASTA_READER_H_

#include <kseq.h>
#include "common.h"

namespace SyntenyFinder
{	
	struct FASTARecord
	{
	public:
		FASTARecord() {}
		FASTARecord(const std::string & sequence, const std::string & description, size_t id):
			description_(description), id_(id)
		{
			sequence_.assign(CFancyIterator(sequence.begin(), tolower, ' '), CFancyIterator(sequence.end(), tolower, ' '));
		}

		size_t GetId() const
		{
			return id_;
		}

		size_t GetConventionalId() const
		{
			return id_ + 1;
		}

		const std::string& GetSequence() const
		{
			return sequence_;
		}

		const std::string& GetDescription() const
		{
			return description_;
		}

	private:
		size_t id_;
		std::string sequence_;
		std::string description_;		
	};

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

		size_t GetSequences(std::vector<FASTARecord> & record);
		
	private:
		DISALLOW_COPY_AND_ASSIGN(FASTAReader);
		FILE * fileHandler_;
	};

	typedef std::vector<FASTARecord> ChrList;

	class FASTAWriter
	{
	public:
		static void WriteSequence(const std::string & fileName, const std::string & header, const std::string & sequence)
		{
			std::ofstream out(fileName.c_str());
			out << ">" << header << std::endl;
			for(size_t i = 0; i < sequence.size(); i += 80)
			{
				size_t j = std::min(i + 80, sequence.size());
				std::copy(sequence.begin() + i, sequence.begin() + j, std::ostream_iterator<char>(out));
				out << std::endl;
			}
		}
	};
}

#endif 
