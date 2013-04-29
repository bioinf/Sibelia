//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "fasta.h"

namespace SyntenyFinder
{	
	namespace
	{
		class ParseException : public std::runtime_error 
		{
		public:
			ParseException(const std::string & what):
				std::runtime_error(what)
			{}
		};
	}

	size_t FASTAReader::GetSequences(std::vector<FASTARecord> & record)
	{
		std::string buffer;
		std::string sequence;
		std::string header;
		int line = 1;
		size_t seqId = 0;

		try
		{
			while(!inputStream_.eof())
			{
				std::getline(inputStream_, buffer);

				if (buffer.empty()) continue;

				if (buffer[0] == '>')
				{
					if (!header.empty())
					{
						if (sequence.empty()) throw ParseException("empty sequence");

						record.push_back(FASTARecord(sequence, header, seqId));
						++seqId;
						sequence.clear();
						header.clear();
					}
					ValidateHeader(buffer);
					header = buffer;
				}
				else
				{
					ValidateSequence(buffer);
					sequence += buffer;
				}

				++line;
			}
			
			if (sequence.empty()) throw ParseException("empty sequence");
			record.push_back(FASTARecord(sequence, header, seqId));
		}
		catch (ParseException & e)
		{
			std::stringstream ss;
			ss << "parse error in " << fileName_ << " on line " << line << ": " << e.what();
			throw std::runtime_error(ss.str());
		}

		return record.size();
	}

	void FASTAReader::ValidateHeader(std::string & header)
	{
		int delim = header.find(' ');
		if (delim == std::string::npos)
		{
			delim = header.length() - 1;
		}
		header = header.substr(1, delim - 1);
		if (header.empty()) throw ParseException("empty header");
	}

	void FASTAReader::ValidateSequence(std::string & sequence)
	{
		const std::string VALID_CHARS = "ACGTURYKMSWBDHWNX-";
		for (size_t i = 0; i < sequence.length(); ++i)
		{
			sequence[i] = toupper(sequence[i]);
			if (VALID_CHARS.find(sequence[i]) == std::string::npos) 
			{
				throw ParseException((std::string("illegal character: ") + sequence[i]).c_str());
			}
		}
	}

	bool FASTAReader::IsOk() const
	{
		return inputStream_.good();
	}
}
