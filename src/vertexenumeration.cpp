//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "blockbuilder.h"

namespace SyntenyFinder
{
	namespace
	{
		class CharSet
		{
		public:
			void Clear()
			{
				set.clear();
			}

			void Add(char ch)
			{
				if(std::find(set.begin(), set.end(), ch) == set.end())
				{
					set.push_back(ch);
				}
			}

			bool In(char ch) const
			{
				return std::find(set.begin(), set.end(), ch) != set.end();
			}

			size_t Size() const
			{
				return set.size();
			}

		private:
			std::string set;
		};

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

		bool Bifurcation(const CharSet & set)
		{
			return set.Size() > 1 || set.In(DeBruijnIndex::SEPARATION_CHAR);
		}

		void Flank(std::string & str, size_t start, size_t end, size_t k, char sepChar)
		{
			if(end - start > k)
			{
				if(!FastaRecord::IsDefiniteBase(str[start]))
				{
					for(; start < end && !FastaRecord::IsDefiniteBase(str[start]); ++start);
					str[start - 1] = sepChar;
				}
				
				if(!FastaRecord::IsDefiniteBase(str[end - k]))
				{
					for(; end - start > k && !FastaRecord::IsDefiniteBase(str[end - k]); --end);
					str[end - k] = sepChar;
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

	size_t BlockBuilder::EnumerateBifurcationsSArray(const std::vector<std::string> & data, size_t k_, const std::string & tempDir, std::vector<DeBruijnIndex::ChrBifVector> & ret) const
	{
		ret.assign(2, DeBruijnIndex::ChrBifVector(data.size()));
		Size bifurcationCount = 0;
		std::vector<size_t> cumSize;
		std::string superGenome(1, DeBruijnIndex::SEPARATION_CHAR);
		for(size_t chr = 0; chr < data.size(); chr++)
		{
			cumSize.push_back(superGenome.size());
			superGenome += data[chr];
			superGenome += DeBruijnIndex::SEPARATION_CHAR;
			Flank(superGenome, superGenome.size() - 1 - data[chr].size(), superGenome.size() - 1, k_, DeBruijnIndex::SEPARATION_CHAR);
		}

		for(size_t chr = 0; chr < data.size(); chr++)
		{
			cumSize.push_back(superGenome.size());
			std::string::const_reverse_iterator it1 = data[chr].rbegin();
			std::string::const_reverse_iterator it2 = data[chr].rend();
			superGenome.insert(superGenome.end(), CFancyIterator(it1, FastaRecord::Translate, ' '), CFancyIterator(it2, FastaRecord::Translate, ' '));
			superGenome += DeBruijnIndex::SEPARATION_CHAR;
			Flank(superGenome, superGenome.size() - 1 - data[chr].size(), superGenome.size() - 1, k_, DeBruijnIndex::SEPARATION_CHAR);
		}
		
		std::vector<saidx_t> pos;
		std::vector<saidx_t> lcp;
		CreateOutDirectory(tempDir);
		FilePtr posFile = CalculateLCP(superGenome, lcp, tempDir);
		CharSet prev;
		CharSet next;
		std::vector<size_t> candidateChr;		
		std::vector<FastaRecord::Direction> candidateDir;
		std::vector<DeBruijnIndex::BifurcationInstance> candidate;
		
		for(size_t start = 0; start < superGenome.size(); )
		{
			pos.assign(1, 0);
			posFile->Read(&pos[0], sizeof(pos[0]), 1);
			if(superGenome[pos[0]] == DeBruijnIndex::SEPARATION_CHAR || !FastaRecord::IsDefiniteBase(superGenome[pos[0]]))
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

			if(Bifurcation(prev) || Bifurcation(next))
			{
				candidate.clear();
				candidateDir.clear();
				candidateChr.clear();
				bool terminal = false;
				for(size_t i = 0; i < pos.size(); i++)
				{
					size_t suffix = pos[i];
					size_t chr = std::upper_bound(cumSize.begin(), cumSize.end(), suffix) - cumSize.begin() - 1;
					FastaRecord::Direction strand = chr < data.size() ? FastaRecord::positive : FastaRecord::negative;
					size_t pos = suffix - cumSize[chr];
					chr = chr < data.size() ? chr : chr - data.size();
					if(pos + k_ <= data[chr].size())
					{
						terminal = terminal || superGenome[suffix - 1] == DeBruijnIndex::SEPARATION_CHAR || superGenome[suffix + k_] == DeBruijnIndex::SEPARATION_CHAR;
						candidate.push_back(DeBruijnIndex::BifurcationInstance(bifurcationCount, pos, pos));
						candidateDir.push_back(strand);
						candidateChr.push_back(chr);
					}
				}

				if(candidate.size() > 1 || terminal)
				{
					bifurcationCount++;
					for(size_t i = 0; i < candidate.size(); i++)
					{
						ret[candidateDir[i] == FastaRecord::positive ? 0 : 1][candidateChr[i]].push_back(candidate[i]);
					}
				}
			}
			
			start += pos.size();
		}
		
		for(size_t dir = 0; dir < 2; dir++)
		{
			for(size_t chr = 0; chr < ret[dir].size(); chr++)
			{
				std::sort(ret[dir][chr].begin(), ret[dir][chr].end());
			}
		}

		return bifurcationCount;
	}
}