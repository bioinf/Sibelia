#include "fasta.h"

namespace SyntenyBuilder
{	
	static size_t cread(FILE * stream, char * buf, size_t size)
	{
		return fread(buf, sizeof(char), size, stream);
	}

	KSEQ_INIT(FILE*, cread)

	void FASTAReader::GetSequences(std::vector<FASTARecord> & record)
	{
		record.clear();
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