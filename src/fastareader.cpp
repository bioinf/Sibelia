#include "fastareader.h"

namespace SyntenyBuilder
{	
	static size_t cread(FILE * stream, char * buf, size_t size)
	{
		return fread(buf, sizeof(char), size, stream);
	}

	KSEQ_INIT(FILE*, cread)

	void FASTAReader::GetSequence(std::string & buffer)
	{
		buffer.clear();
		kseq_t * sequence = kseq_init(fileHandler_);
		if(kseq_read(sequence) >= 0)
		{
			buffer = sequence->seq.s;
		}

		kseq_destroy(sequence);
	}
}