//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "fasta.h"
#include <boost/algorithm/string.hpp>

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

	size_t FastaReader::GetSequences(std::vector<FastaRecord> & record)
	{
		std::string buffer;
		std::string sequence;
		std::string header;
		int line = 1;
		size_t seqId = record.size();

		try
		{
			while(!inputStream_.eof())
			{
				std::getline(inputStream_, buffer);
				boost::algorithm::trim(buffer);

				if (buffer.empty()) continue;

				if (buffer[0] == '>')
				{
					if (!header.empty())
					{
						if (sequence.empty()) throw ParseException("empty sequence");

						record.push_back(FastaRecord(sequence, header, seqId));
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
			record.push_back(FastaRecord(sequence, header, seqId));
		}
		catch (ParseException & e)
		{
			std::stringstream ss;
			ss << "parse error in " << fileName_ << " on line " << line << ": " << e.what();
			throw std::runtime_error(ss.str());
		}

		return record.size();
	}

	void FastaReader::ValidateHeader(std::string & header)
	{
		size_t delim = header.find(' ');
		if (delim == std::string::npos)
		{
			delim = header.length() - 1;
		}
		else
		{
			--delim;
		}

		header = header.substr(1, delim);
		if (header.empty()) throw ParseException("empty header");
	}

	void FastaReader::ValidateSequence(std::string & sequence)
	{
		const std::string VALID_CHARS = "ACGTURYKMSWBDHWNX-";
		for (size_t i = 0; i < sequence.length(); ++i)
		{
			char orig = sequence[i];
			sequence[i] = toupper(sequence[i]);
			if (VALID_CHARS.find(sequence[i]) == std::string::npos) 
			{
				throw ParseException((std::string("illegal character: ") + orig).c_str());
			}
		}
	}

	bool FastaReader::IsOk() const
	{
		return inputStream_.good();
	}
	 
	size_t FastaRecord::Iterator::GetPosition() const
	{
		return pos_;
	}

	FastaRecord::Iterator& FastaRecord::Iterator::operator++()
	{
		pos_++;
		return *this;
	}

	
	FastaRecord::Iterator& FastaRecord::Iterator::operator--()
	{
		pos_--;
		return *this;
	}
	
	FastaRecord::Iterator FastaRecord::Iterator::operator++(int)
	{
		Iterator ret(*this);
		++(*this);
		return ret;
	}

	FastaRecord::Iterator::Iterator(const FastaRecord * chr, size_t pos, Direction dir)
	{
	}

	FastaRecord::Iterator FastaRecord::Iterator::operator--(int)
	{
		Iterator ret(*this);
		--(*this);
		return ret;
	}
	
	FastaRecord::Direction FastaRecord::Iterator::GetDirection() const
	{
		return dir_;
	}
	
	const FastaRecord* FastaRecord::Iterator::GetSequence() const
	{
		return chr_;
	}
	
	FastaRecord::Iterator FastaRecord::Iterator::operator + (size_t shift) const
	{
		return Iterator(chr_, pos_ + shift, dir_);
	}
	
	FastaRecord::Iterator FastaRecord::Iterator::operator - (size_t shift) const
	{
		return Iterator(chr_, pos_ - shift, dir_);
	}

	FastaRecord::Iterator FastaRecord::Begin(Direction dir) const
	{
		return Iterator(this, 0, dir);
	}

	FastaRecord::Iterator FastaRecord::End(Direction dir) const
	{
		return Iterator(this, sequence_.size(), dir);
	}
}
