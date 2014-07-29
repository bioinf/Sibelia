//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _FASTA_READER_H_
#define _FASTA_READER_H_

#include "common.h"

namespace SyntenyFinder
{	
	struct FASTARecord
	{
	public:
		FASTARecord() {}
		FASTARecord(const std::string & sequence, const std::string & description, size_t id):
			description_(description), id_(id), sequence_(sequence)
		{
			//sequence_.assign(CFancyIterator(sequence.begin(), toupper, ' '), CFancyIterator(sequence.end(), toupper, ' '));
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

		std::string GetStripedId() const
		{
			std::string description(description_);
			for(std::string::iterator it = description.begin(); it != description.end(); ++it)
			{
				if(*it == '|' || *it == '.')
				{
					*it = ' ';
				}
			}
			
			std::string buf;
			std::stringstream ss(description);
			std::vector<std::string> token;
			while(ss >> buf)
			{
				token.push_back(buf);
			}

			return token.size() == 5 ? token[3] : description;
		}		

		enum Direction
		{
			positive,
			negative
		};

		class Iterator
		{
		public:			
			Iterator& operator++();
			Iterator& operator--();
			Iterator operator++(int);
			Iterator operator--(int);
			size_t GetPosition() const;			
			Direction GetDirection() const;
			const FASTARecord GetSequence() const;
			Iterator operator + (size_t shift) const;
			Iterator operator - (size_t shift) const;			
		private:
			const FASTARecord * chr;
			size_t pos;
			Direction dir;			
		};


		Iterator Begin(Direction dir) const;
		Iterator End(Direction dir) const;

	private:
		size_t id_;
		std::string sequence_;
		std::string description_;		
	};

	class FASTAReader
	{
	public:
		explicit FASTAReader(const std::string & fileName): 
			inputStream_(fileName.c_str()),
			fileName_(fileName) {}
		size_t 	GetSequences(std::vector<FASTARecord>& record);
		bool 	IsOk() const;
	private:
		DISALLOW_COPY_AND_ASSIGN(FASTAReader);

		void ValidateSequence(std::string & sequence);
		void ValidateHeader(std::string & header);

		std::ifstream inputStream_;
		std::string fileName_;
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
