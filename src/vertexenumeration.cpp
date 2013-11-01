//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "indexedsequence.h"

namespace SyntenyFinder
{
	const char IndexedSequence::SEPARATION_CHAR = '#';
	namespace
	{
		void GetHeight(const std::string & str, const std::vector<saidx_t> & order, const std::vector<saidx_t> & pos, std::vector<Size> & height)
		{
			size_t n = pos.size();
			height.assign(n, 0);
			for (size_t i = 0, h = 0; i < n; ++i)
			{
				if (pos[i] > 0)
				{
					size_t j = order[pos[i] - 1];
					while (i + h < n && j + h < n && str[i + h] == str[j + h])
					{
						h++;
					}

					height[pos[i]] = static_cast<saidx_t>(h);
					if (h > 0)
					{
						h--;
					}
				}
			}
		}
		
		size_t StupidLCP(const std::string & str, size_t suffix1, size_t suffix2)
		{
			size_t ret = 0;
			for(; suffix1 + ret < str.size() && suffix2 + ret < str.size() && str[suffix1 + ret] == str[suffix2 + ret]; ++ret);
			return ret;
		}

		typedef boost::shared_ptr<TempFile> FilePtr;

		FilePtr CreateFileWithSA(const std::string & superGenome, const std::string & tempDir)
		{
			FilePtr posFile(new TempFile(tempDir));
			std::vector<saidx_t> pos(superGenome.size());
			divsufsort(reinterpret_cast<const sauchar_t*>(superGenome.c_str()), &pos[0], static_cast<saidx_t>(pos.size()));
			posFile->Write(&pos[0], sizeof(pos[0]), pos.size());
			posFile->Rewind();
			return posFile;
		}

		void FindPhi(std::vector<saidx_t> & phi, FilePtr posFile)
		{
			saidx_t pos[2];
			posFile->Read(pos, sizeof(pos[0]), 1);
			for(size_t i = 1; i < phi.size(); i++)
			{
				posFile->Read(pos + 1, sizeof(pos[0]), 1);
				phi[pos[1]] = pos[0];
				pos[0] = pos[1];
			}

			posFile->Rewind();
		}

		FilePtr CalculateLCP(const std::string & superGenome, std::vector<saidx_t> & lcp, const std::string & tempDir)
		{
			FilePtr lcpFile(new TempFile(tempDir));
			FilePtr posFile = CreateFileWithSA(superGenome, tempDir);
			{
				std::vector<saidx_t> phi(superGenome.size(), 0);
				FindPhi(phi, posFile);

				saidx_t l = 0;
				for(size_t i = 0; i < superGenome.size(); i++)
				{
					while(superGenome[i + l] == superGenome[phi[i] + l])
					{
						++l;
					}

					phi[i] = l;
					l = std::max(l - 1, 0);
				}

				saidx_t pos;
				for(size_t i = 0; i < superGenome.size(); i++)
				{
					posFile->Read(&pos, sizeof(pos), 1);
					saidx_t lcp = phi[pos];
					lcpFile->Write(&lcp, sizeof(lcp), 1);
				}
			}
			
			lcpFile->Rewind();
			lcp.assign(superGenome.size(), 0);
			lcpFile->Read(&lcp[0], sizeof(lcp[0]), lcp.size());			
			posFile->Rewind();
			return posFile;
		}		
	}

	size_t IndexedSequence::EnumerateBifurcationsSArray(const std::vector<std::string> & data, const std::string & tempDir, std::vector<BifurcationInstance> & positiveBif, std::vector<BifurcationInstance> & negativeBif)
	{
		positiveBif.clear();
		negativeBif.clear();
		Size bifurcationCount = 0;
		std::vector<size_t> cumSize;
		std::string superGenome(1, SEPARATION_CHAR);
		for(size_t chr = 0; chr < data.size(); chr++)
		{
			cumSize.push_back(superGenome.size());
			superGenome += data[chr];
			superGenome += SEPARATION_CHAR;
		}

		for(size_t chr = 0; chr < data.size(); chr++)
		{
			cumSize.push_back(superGenome.size());
			std::string::const_reverse_iterator it1 = data[chr].rbegin();
			std::string::const_reverse_iterator it2 = data[chr].rend();
			superGenome.insert(superGenome.end(), CFancyIterator(it1, DNASequence::Translate, ' '), CFancyIterator(it2, DNASequence::Translate, ' '));
			superGenome += SEPARATION_CHAR;
		}
		
		std::vector<saidx_t> pos;
		std::vector<saidx_t> lcp;
		CreateOutDirectory(tempDir);
		FilePtr posFile = CalculateLCP(superGenome, lcp, tempDir);
		CharSet prev;
		CharSet next;
		std::vector<BifurcationInstance> * ret[] = {&positiveBif, &negativeBif};
		std::vector<std::pair<DNASequence::Direction, BifurcationInstance> > candidate;
		for(size_t start = 0; start < superGenome.size(); )
		{
			pos.assign(1, 0);
			posFile->Read(&pos[0], sizeof(pos[0]), 1);
			if(superGenome[pos[0]] == SEPARATION_CHAR || !IsDefiniteBase(superGenome[pos[0]]))
			{
				start++;
				continue;
			}

			prev.Clear();
			next.Clear();
			size_t match = 0;
			do
			{				
				if(match > 0)
				{
					pos.push_back(saidx_t());
					posFile->Read(&pos.back(), sizeof(pos[0]), 1);
					assert(lcp[start + match] == StupidLCP(superGenome, pos[match], pos[match - 1]));
				}

				if(pos[match] > 0)
				{
					prev.Add(superGenome[pos[match] - 1]);
				}

				if(pos[match] + k_ < superGenome.size())
				{
					next.Add(superGenome[pos[match] + k_]);
				}
			}
			while(++match + start < superGenome.size() && lcp[match + start] >= static_cast<size_t>(k_));

			if(prev.Bifurcation() || next.Bifurcation())
			{
				candidate.clear();
				bool terminal = false;
				for(size_t i = 0; i < pos.size(); i++)
				{
					size_t suffix = pos[i];
					size_t chr = std::upper_bound(cumSize.begin(), cumSize.end(), suffix) - cumSize.begin() - 1;
					DNASequence::Direction strand = chr < data.size() ? DNASequence::positive : DNASequence::negative;
					size_t pos = suffix - cumSize[chr];
					chr = chr < data.size() ? chr : chr - data.size();
					if(pos + k_ <= data[chr].size())
					{
						terminal = terminal || superGenome[suffix - 1] == SEPARATION_CHAR || superGenome[suffix + k_] == SEPARATION_CHAR;
						candidate.push_back(std::make_pair(strand, BifurcationInstance(bifurcationCount, static_cast<Size>(chr), static_cast<Size>(pos))));
					}
				}

				if(candidate.size() > 1 || terminal)
				{
					bifurcationCount++;
					for(size_t i = 0; i < candidate.size(); i++)
					{
						ret[candidate[i].first == DNASequence::positive ? 0 : 1]->push_back(candidate[i].second);
					}
				}
			}
			
			start += pos.size();
		}
				
		std::sort(positiveBif.begin(), positiveBif.end());
		std::sort(negativeBif.begin(), negativeBif.end());
		return bifurcationCount;
	}

	size_t IndexedSequence::EnumerateBifurcationsSArrayInRAM(const std::vector<std::string> & data, std::vector<BifurcationInstance> & positiveBif, std::vector<BifurcationInstance> & negativeBif)
	{
		positiveBif.clear();
		negativeBif.clear();
		Size bifurcationCount = 0;
		std::vector<size_t> cumSize;
		std::string superGenome(1, SEPARATION_CHAR);
		for(size_t chr = 0; chr < data.size(); chr++)
		{
			cumSize.push_back(superGenome.size());
			superGenome += data[chr];
			superGenome += SEPARATION_CHAR;
		}

		for(size_t chr = 0; chr < data.size(); chr++)
		{
			cumSize.push_back(superGenome.size());
			std::string::const_reverse_iterator it1 = data[chr].rbegin();
			std::string::const_reverse_iterator it2 = data[chr].rend();
			superGenome.insert(superGenome.end(), CFancyIterator(it1, DNASequence::Translate, ' '), CFancyIterator(it2, DNASequence::Translate, ' '));
			superGenome += SEPARATION_CHAR;
		}

		std::vector<Size> lcp;
		std::vector<saidx_t> order(superGenome.size());
		{
			std::vector<saidx_t> pos(superGenome.size());
			divsufsort(reinterpret_cast<const sauchar_t*>(superGenome.c_str()), &order[0], static_cast<saidx_t>(order.size()));
			for(size_t i = 0; i < order.size(); i++)
			{
				pos[order[i]] = static_cast<saidx_t>(i);
			}

			GetHeight(superGenome, order, pos, lcp);	
		}

		CharSet prev;
		CharSet next;
		std::vector<BifurcationInstance> * ret[] = {&positiveBif, &negativeBif};
		std::vector<std::pair<DNASequence::Direction, BifurcationInstance> > candidate;
		for(size_t start = 0; start < superGenome.size(); )
		{
			if(superGenome[order[start]] == SEPARATION_CHAR  || !IsDefiniteBase(superGenome[order[start]]))
			{
				start++;
				continue;
			}

			prev.Clear();
			next.Clear();
			size_t end = start;
			do
			{
				if(order[end] > 0)
				{
					prev.Add(superGenome[order[end] - 1]);
				}

				if(order[end] + k_ < superGenome.size())
				{
					next.Add(superGenome[order[end] + k_]);
				}
			}
			while(++end < superGenome.size() && lcp[end] >= k_);

			if(prev.Bifurcation() || next.Bifurcation())
			{
				candidate.clear();
				bool terminal = false;
				for(size_t i = start; i < end; i++)
				{
					size_t suffix = order[i];
					size_t chr = std::upper_bound(cumSize.begin(), cumSize.end(), suffix) - cumSize.begin() - 1;
					DNASequence::Direction strand = chr < data.size() ? DNASequence::positive : DNASequence::negative;
					size_t pos = suffix - cumSize[chr];
					chr = chr < data.size() ? chr : chr - data.size();
					if(pos + k_ <= data[chr].size())
					{
						terminal = terminal || superGenome[suffix - 1] == SEPARATION_CHAR || superGenome[suffix + k_] == SEPARATION_CHAR;
						candidate.push_back(std::make_pair(strand, BifurcationInstance(bifurcationCount, static_cast<Size>(chr), static_cast<Size>(pos))));
					}
				}

				if(candidate.size() > 1 || terminal)
				{
					bifurcationCount++;
					for(size_t i = 0; i < candidate.size(); i++)
					{
						ret[candidate[i].first == DNASequence::positive ? 0 : 1]->push_back(candidate[i].second);
					}
				}
			}
			
			start = end;
		}
		
		std::sort(positiveBif.begin(), positiveBif.end());
		std::sort(negativeBif.begin(), negativeBif.end());
		return bifurcationCount;
	}
}