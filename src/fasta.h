#ifndef _FASTA_READER_H_
#define _FASTA_READER_H_

#include "lib/kseq.h"
#include "common.h"

namespace SyntenyFinder
{
	//This class provides functionality for reading FASTA files.
	//Currently this is an OO wrapper around tiny library kseq.h.
	//It looks little ugly, maybe later it will be replaced with
	//some other library or written by myself.

	class FASTAReader
	{
	public:
		struct FASTARecord
		{
			std::string sequence;
			std::string description;
			FASTARecord() {}
			FASTARecord(const std::string & sequence, const std::string & description):
				sequence(sequence), description(description) {}
		};

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

		void GetSequences(std::vector<FASTARecord> & record);
		
	private:
		DISALLOW_COPY_AND_ASSIGN(FASTAReader);
		FILE * fileHandler_;
	};

	typedef std::vector<FASTAReader::FASTARecord> ChrList;

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
