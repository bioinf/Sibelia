//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "fasta.h"

namespace SyntenyFinder
{	
	static size_t cread(FILE * stream, char * buf, size_t size)
	{
		return fread(buf, sizeof(char), size, stream);
	}

	KSEQ_INIT(FILE*, cread)

	void FASTAReader::GetSequences(std::vector<FASTARecord> & record)
	{
		kseq_t * sequence = kseq_init(fileHandler_);
		while(kseq_read(sequence) >= 0)
		{
			record.push_back(FASTARecord(sequence->seq.s, sequence->name.s));
			for(size_t i = 0; i < record.back().sequence.size(); i++)
			{
				record.back().sequence[i] = tolower(record.back().sequence[i]);
			}	
		}

		kseq_destroy(sequence);
	}
}