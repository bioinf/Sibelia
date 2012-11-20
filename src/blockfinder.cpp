//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "blockfinder.h"

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

		bool Bifurcation(const CharSet & set)
		{
			return set.Size() > 1 || set.In(BlockFinder::SEPARATION_CHAR);
		}

		void Flank(std::string & str, size_t start, size_t end, size_t k, char sepChar)
		{
			if(end - start > k)
			{
				if(str[start] == 'n')
				{
					for(; start < str.size() && str[start] == 'n'; ++start);
					str[start - 1] = sepChar;
				}

				if(str[end - k] == 'n')
				{
					for(; end - start > k && str[end - k] == 'n'; --end);
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

		void ShrinkToEmpty(std::vector<saidx_t> & toShrink)
		{
			std::vector<saidx_t> empty;
			toShrink.clear();
			toShrink.swap(empty);
		}

		struct Handler
		{
		public:			
			void AddHandle(FILE * handle)
			{
				toClose.push_back(handle);
			}

			~Handler()
			{
				std::for_each(toClose.begin(), toClose.end(), fclose);
			}

			void Clear()
			{
				toClose.clear();
			}

		private:
			std::vector<FILE*> toClose;
		};

		FILE * TempFile()
		{
			FILE * ret = tmpfile();
			if(ret == 0)
			{
				throw std::runtime_error("Can't create a temporary file");
			}

			return ret;
		}

		void Write(const void * ptr, size_t size, size_t count, FILE * handle)
		{
			if(fwrite(ptr, size, count, handle) != count)
			{
				throw std::runtime_error("Error while writing to a temporary file");
			}

		}

		void Read(void * ptr, size_t size, size_t count, FILE * handle)
		{
			if(fread(ptr, size, count, handle) != count)
			{
				throw std::runtime_error("Error while reading from a temporary file");
			}
		}

		FILE* CreateFileWithSA(const std::string & superGenome)
		{
			FILE * posFile = TempFile();
			std::vector<saidx_t> pos(superGenome.size());
			divsufsort(reinterpret_cast<const sauchar_t*>(superGenome.c_str()), &pos[0], static_cast<saidx_t>(pos.size()));
			Write(&pos[0], sizeof(pos[0]), pos.size(), posFile);
			rewind(posFile);
			return posFile;
		}

		void FindPhi(std::vector<saidx_t> & phi, FILE * posFile)
		{
			saidx_t pos[2];
			Read(pos, sizeof(pos[0]), 1, posFile);
			for(size_t i = 1; i < phi.size(); i++)
			{
				Read(pos + 1, sizeof(pos[0]), 1, posFile);
				phi[pos[1]] = pos[0];
				pos[0] = pos[1];
			}

			rewind(posFile);
		}
		

		void CalculateLCP(const std::string & superGenome, std::vector<saidx_t> & lcp, FILE* & posFile)
		{
			FILE * lcpFile;
			Handler handler;
			handler.AddHandle(lcpFile = TempFile());
			handler.AddHandle(posFile = CreateFileWithSA(superGenome));
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
					Read(&pos, sizeof(pos), 1, posFile);
					saidx_t lcp = phi[pos];
					Write(&lcp, sizeof(lcp), 1, lcpFile);
				}
			}
			
			rewind(lcpFile);
			lcp.assign(superGenome.size(), 0);
			Read(&lcp[0], sizeof(lcp[0]), lcp.size(), lcpFile);
			fclose(lcpFile);
			rewind(posFile);
			handler.Clear();
		}		
	}

	const char BlockFinder::SEPARATION_CHAR = '#';

	BlockFinder::BlockFinder(const std::vector<FASTARecord> & chrList): chrList_(chrList)
	{
		originalPos_.resize(chrList.size());
		for(size_t i = 0; i < originalPos_.size(); i++)
		{
			originalPos_[i].resize(chrList[i].sequence.size());
			for(size_t j = 0; j < originalPos_[i].size(); j++)
			{
				originalPos_[i][j] = static_cast<Pos>(j);
			}
		}
	}

	size_t BlockFinder::EnumerateBifurcationsSArray(size_t k, std::vector<BifurcationInstance> & positiveBif, std::vector<BifurcationInstance> & negativeBif) const
	{
		positiveBif.clear();
		negativeBif.clear();
		Size bifurcationCount = 0;
		std::vector<size_t> cumSize;
		std::string superGenome(1, SEPARATION_CHAR);
		for(size_t chr = 0; chr < chrList_.size(); chr++)
		{
			cumSize.push_back(superGenome.size());
			superGenome += chrList_[chr].sequence;
			superGenome += SEPARATION_CHAR;
			Flank(superGenome, superGenome.size() - 1 - chrList_[chr].sequence.size(), superGenome.size() - 1, k, SEPARATION_CHAR);
		}

		for(size_t chr = 0; chr < chrList_.size(); chr++)
		{
			cumSize.push_back(superGenome.size());
			std::string::const_reverse_iterator it1 = chrList_[chr].sequence.rbegin();
			std::string::const_reverse_iterator it2 = chrList_[chr].sequence.rend();
			superGenome.insert(superGenome.end(), CFancyIterator(it1, DNASequence::Translate, ' '), CFancyIterator(it2, DNASequence::Translate, ' '));
			superGenome += SEPARATION_CHAR;
			Flank(superGenome, superGenome.size() - 1 - chrList_[chr].sequence.size(), superGenome.size() - 1, k, SEPARATION_CHAR);
		}
		
		FILE * posFile;
		std::vector<saidx_t> pos;
		std::vector<saidx_t> lcp;
		CalculateLCP(superGenome, lcp, posFile);
		
		CharSet prev;
		CharSet next;
		std::vector<BifurcationInstance> * ret[] = {&positiveBif, &negativeBif};
		std::vector<std::pair<DNASequence::Direction, BifurcationInstance> > candidate;
		for(size_t start = 0; start < superGenome.size(); )
		{
			pos.assign(1, 0);
			Read(&pos[0], sizeof(pos[0]), 1, posFile);
			if(superGenome[pos[0]] == SEPARATION_CHAR || superGenome[pos[0]] == 'n')
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
					Read(&pos.back(), sizeof(pos[0]), 1, posFile);
					assert(lcp[start + match] == StupidLCP(superGenome, pos[match], pos[match - 1]));
				}

				if(pos[match] > 0)
				{
					prev.Add(superGenome[pos[match] - 1]);
				}

				if(pos[match] + k < superGenome.size())
				{
					next.Add(superGenome[pos[match] + k]);
				}
			}
			while(++match + start < superGenome.size() && lcp[match + start] >= k);

			if(Bifurcation(prev) || Bifurcation(next))
			{
				candidate.clear();
				bool terminal = false;
				for(size_t i = 0; i < pos.size(); i++)
				{
					size_t suffix = pos[i];
					size_t chr = std::upper_bound(cumSize.begin(), cumSize.end(), suffix) - cumSize.begin() - 1;
					DNASequence::Direction strand = chr < chrList_.size() ? DNASequence::positive : DNASequence::negative;
					size_t pos = suffix - cumSize[chr];
					chr = chr < chrList_.size() ? chr : chr - chrList_.size();
					if(pos + k <= chrList_[chr].sequence.size())
					{
						terminal = terminal || superGenome[suffix - 1] == SEPARATION_CHAR || superGenome[suffix + k] == SEPARATION_CHAR;
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
		
		fclose(posFile);
		std::sort(positiveBif.begin(), positiveBif.end());
		std::sort(negativeBif.begin(), negativeBif.end());
		return bifurcationCount;
	}

	void BlockFinder::ConstructBifStorage(const DNASequence & sequence, const std::vector<std::vector<BifurcationInstance> > & bifurcation, BifurcationStorage & bifStorage) const
	{
		for(size_t strand = 0; strand < 2; strand++)
		{
			size_t nowBif = 0;
			DNASequence::Direction dir = static_cast<DNASequence::Direction>(strand);
			for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
			{
				size_t pos = 0;
				StrandIterator end = sequence.End(dir, chr);
				for(DNASequence::StrandIterator it = sequence.Begin(dir, chr); it != end; ++it, ++pos)
				{
					if(nowBif < bifurcation[strand].size() && chr == bifurcation[strand][nowBif].chr && pos == bifurcation[strand][nowBif].pos)
					{
						bifStorage.AddPoint(it, bifurcation[strand][nowBif++].bifId);
					}
				}
			}
		}

	}

	void BlockFinder::PerformGraphSimplifications(size_t k, size_t minBranchSize, size_t maxIterations, ProgressCallBack f)
	{
	#ifdef NEW_ENUMERATION		
		std::vector<std::vector<BifurcationInstance> > bifurcation(2);		
		size_t maxId = EnumerateBifurcationsSArray(k, bifurcation[0], bifurcation[1]);
		DNASequence sequence(chrList_, originalPos_, true);
		{
			BifurcationStorage bifStorage(maxId);
			ConstructBifStorage(sequence, bifurcation, bifStorage);
		#else
			BifurcationStorage bifStorage;
			DNASequence sequence(chrList_, originalPos_);
			EnumerateBifurcationsHash(sequence, bifStorage, k);
		#endif

		#ifdef _DEBUG
			bifStorage.FormDictionary(idMap, k);
		#endif
			SimplifyGraph(sequence, bifStorage, k, minBranchSize, maxIterations, f);
		}

		for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
		{
			originalPos_[chr].clear();
			chrList_[chr].sequence.clear();
			StrandIterator end = sequence.PositiveEnd(chr);
			for(StrandIterator it = sequence.PositiveBegin(chr); it != end; ++it)
			{
				chrList_[chr].sequence.push_back(*it);
				originalPos_[chr].push_back(static_cast<Pos>(it.GetOriginalPosition()));
			}
		}
	}	
}