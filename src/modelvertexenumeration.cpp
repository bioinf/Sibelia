//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "indexedsequence.h"

namespace SyntenyFinder
{
	namespace
	{
		size_t BaseMask(char ch)
		{
			switch(ch)
			{
			case 'A':
				return 0;
			case 'C':
				return 1;
			case 'G':
				return 2;
			case 'T':
				return 3;
			}

			assert(false);
		}

		size_t KMerToNumber(std::string::const_iterator it, size_t k)
		{
			size_t ret = 0;
			for(size_t pos = 0; pos < k; ++it, ++pos)
			{
				ret = ret | (BaseMask(*it) << (2 * pos));
			}

			return ret;
		}

		size_t Shift(size_t kmer, char next, size_t k)
		{
			return (kmer >> 2) | (BaseMask(next) << (2 * (k - 1)));
		}		
	}
	
	size_t IndexedSequence::EnumerateBifurcationsModel(const std::vector<std::string> & data, std::vector<BifurcationInstance> & posBifurcation, std::vector<BifurcationInstance> & negBifurcation, size_t model)
	{		
		posBifurcation.clear();
		negBifurcation.clear();
		std::map<size_t, CharSet> prevCharSet;
		std::map<size_t, CharSet> nextCharSet;
		for(std::vector<std::string>::const_iterator it = data.begin(); it < data.end(); ++it)
		{
			if(it->size() >= k_)
			{
				std::string revcomp(CFancyIterator(it->begin(), DNASequence::Translate, ' '), CFancyIterator(it->end(), DNASequence::Translate, ' '));
				const std::string * str[] = {&(*it), &revcomp};
				for(size_t ptr = 0; ptr < 2; ptr++)
				{
					size_t kmer = KMerToNumber(str[ptr]->begin(), k_);
					for(size_t pos = 0; pos + k_ <= str[ptr]->size(); ++pos)
					{
						char prevChar = pos == 0 ? SEPARATION_CHAR : (*str[ptr])[pos - 1];
						char nextChar = pos + k_ < str[ptr]->size() ? (*str[ptr])[pos + k_] : SEPARATION_CHAR;
						size_t eqKmer = kmer & model;
						prevCharSet[eqKmer].Add(prevChar);
						nextCharSet[eqKmer].Add(nextChar);
						if(nextChar != SEPARATION_CHAR)
						{
							kmer = Shift(kmer, nextChar, k_);
							assert(kmer == KMerToNumber(str[ptr]->begin() + pos + 1, k_));
						}
					}
				}
			}
		}

		std::vector<size_t> bifurcation;
		std::map<size_t, CharSet> * charSet[2] = {&prevCharSet, &nextCharSet};
		for(size_t location = 0; location < 2; location++)
		{
			for(std::map<size_t, CharSet>::iterator it = charSet[location]->begin(); it != charSet[location]->end(); ++it)
			{
				if(it->second.Bifurcation())
				{
					bifurcation.push_back(it->first);
				}
			}
		}

		std::sort(bifurcation.begin(), bifurcation.end());
		bifurcation.erase(std::unique(bifurcation.begin(), bifurcation.end()), bifurcation.end());

		for(std::vector<std::string>::const_iterator it = data.begin(); it < data.end(); ++it)
		{
			if(it->size() >= k_)
			{
				std::string revcomp(CFancyIterator(it->begin(), DNASequence::Translate, ' '), CFancyIterator(it->end(), DNASequence::Translate, ' '));
				const std::string * str[] = {&(*it), &revcomp};
				std::vector<BifurcationInstance> * bifurcationPosition[] = {&posBifurcation, &negBifurcation};
				for(size_t ptr = 0; ptr < 2; ptr++)
				{
					size_t kmer = KMerToNumber(str[ptr]->begin(), k_);
					for(size_t pos = 0; pos + k_ <= str[ptr]->size(); ++pos)
					{
						size_t eqKmer = kmer & model;
						size_t bifId = std::lower_bound(bifurcation.begin(), bifurcation.end(), eqKmer) - bifurcation.begin();
						if(bifId != bifurcation.size() && bifurcation[bifId] == eqKmer)
						{
							bifurcationPosition[ptr]->push_back(BifurcationInstance(bifId, it - data.begin(), pos));
						}

						char nextChar = pos + k_ < str[ptr]->size() ? (*str[ptr])[pos + k_] : SEPARATION_CHAR;
						if(nextChar != SEPARATION_CHAR)
						{
							kmer = Shift(kmer, nextChar, k_);
						}
					}
				}
			}
		}

		return bifurcation.size();
	}
}