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

		void MatchBifurcationStorage(const DNASequence & sequence, const BifurcationStorage & now, const BifurcationStorage & old)
		{			
			const size_t NO_BIF = BifurcationStorage::NO_BIFURCATION;
			for(size_t strand = 0; strand < 2; strand++)
			{
				DNASequence::Direction dir = static_cast<DNASequence::Direction>(strand);
				for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
				{
					for(StrandIterator it = sequence.Begin(dir, chr); it != sequence.End(dir, chr); ++it)
					{
						std::vector<StrandIterator> nowPos;
						std::vector<StrandIterator> oldPos;
						size_t oldBifId = old.GetBifurcation(it);
						size_t nowBifId = now.GetBifurcation(it);
						assert((oldBifId == NO_BIF && nowBifId == NO_BIF) || (oldBifId != NO_BIF && nowBifId != NO_BIF));						
						now.ListPositions(nowBifId, std::back_inserter(nowPos));
						old.ListPositions(oldBifId, std::back_inserter(oldPos));
						assert(nowPos.size() == oldPos.size());
						for(size_t i = 0; i < nowPos.size(); i++)
						{
							assert(std::find(oldPos.begin(), oldPos.end(), nowPos[i]) != oldPos.end());
						}
					}
				}
			}
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

	void BlockFinder::NotifyBefore(NotificationData data, StrandIterator begin, StrandIterator end)
	{
		nowInvalid_ = 0;
		invalid_.push_back(std::vector<size_t>());
		data.bifStorage->NotifyBefore(begin, end);		
		for(StrandIterator it = begin; it != end; ++it)
		{			
			IteratorIndexMap::iterator index = data.iteratorIndex->find(it);
			if(index != data.iteratorIndex->end())
			{
				invalid_.back().push_back(index->second);
				RemoveRestricted(*data.restricted, it, index->second, data.k);
				data.iteratorIndex->erase(it);
			}
			else
			{
				invalid_.back().push_back(UNUSED);
			}
		}
	}

	void BlockFinder::NotifyAfter(NotificationData data, StrandIterator begin, StrandIterator end)
	{
		size_t pos = 0;
		data.bifStorage->NotifyAfter(begin, end);
		for(StrandIterator it = begin; it != end; ++it, ++pos)
		{
			if(invalid_[nowInvalid_][pos] != UNUSED)
			{				
				(*data.startKMer)[invalid_[nowInvalid_][pos]] = it;
				AddRestricted(*data.restricted, (*data.startKMer)[invalid_[nowInvalid_][pos]], invalid_[nowInvalid_][pos], data.k);
				(*data.iteratorIndex)[it] = invalid_[nowInvalid_][pos];
			}
		}

		if(++nowInvalid_ == invalid_.size())
		{
			invalid_.clear();
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

		std::vector<Size> lcp;
		std::vector<saidx_t> pos(superGenome.size());
		std::vector<saidx_t> order(superGenome.size());
		std::vector<BifurcationInstance> * ret[] = {&positiveBif, &negativeBif};
		divsufsort(reinterpret_cast<const sauchar_t*>(superGenome.c_str()), &order[0], static_cast<saidx_t>(order.size()));
		for(size_t i = 0; i < order.size(); i++)
		{
			pos[order[i]] = static_cast<saidx_t>(i);
		}

		CharSet prev;
		CharSet next;
		GetHeight(superGenome, order, pos, lcp);
		std::vector<std::pair<DNASequence::Direction, BifurcationInstance> > candidate;
		for(size_t start = 0; start < superGenome.size(); )
		{
			if(superGenome[order[start]] == '$' || superGenome[order[start]] == 'n')
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

				if(order[end] + k < superGenome.size())
				{
					next.Add(superGenome[order[end] + k]);
				}
			}
			while(++end < superGenome.size() && lcp[end] >= k);

			if(Bifurcation(prev) || Bifurcation(next))
			{
				candidate.clear();
				bool terminal = false;
				for(size_t i = start; i < end; i++)
				{
					size_t suffix = order[i];
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
			
			start = end;
		}
		
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
		BifurcationStorage bifStorage;
		std::vector<std::vector<BifurcationInstance> > bifurcation(2);		
		EnumerateBifurcationsSArray(k, bifurcation[0], bifurcation[1]);
		DNASequence sequence(chrList_, originalPos_);
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