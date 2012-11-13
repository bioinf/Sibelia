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
		const size_t PROGRESS_STRIDE = 50;

		struct BifurcationData
		{
		public:
			static const char NO_CHAR;
			static const size_t FORWARD;
			static const size_t BACKWARD;
			typedef BifurcationStorage::BifurcationId BifurcationId;

			BifurcationData(BifurcationId id = BifurcationStorage::NO_BIFURCATION): id_(id), forward_(NO_CHAR), backward_(NO_CHAR) {}
			bool Update(StrandIterator it, size_t direction)
			{
				char BifurcationData::*field[2] = 
				{
					&BifurcationData::forward_,
					&BifurcationData::backward_
				};

				if(id_ == BifurcationStorage::NO_BIFURCATION && it.AtValidPosition())
				{
					if(this->*field[direction] == NO_CHAR)
					{
						this->*field[direction]= *it;
					}
					else if(this->*field[direction] != *it)
					{
						return true;
					}
				}

				return false;
			}
			
			void SetId(BifurcationId newId)
			{
				id_ = newId;
			}

			size_t GetId() const
			{
				return id_;
			}

		private:
			BifurcationData::BifurcationId id_;
			char forward_;
			char backward_;
		};

		const char BifurcationData::NO_CHAR = -1;
		const size_t BifurcationData::FORWARD = 0;
		const size_t BifurcationData::BACKWARD = 1;

		StrandIterator ForwardValidKMer(StrandIterator it, size_t k)
		{
			for(; it.AtValidPosition() && *it == DNASequence::UNKNOWN_BASE; ++it);
			if(!ProperKMer(it, k))
			{
				for(; !it.AtValidPosition(); ++it);
			}

			return it;
		}

		StrandIterator BackwardValidKMer(StrandIterator it, size_t k)
		{
			for(; it.AtValidPosition() && *it == DNASequence::UNKNOWN_BASE; --it);			
			return it;
		}
	}
	
	
#ifdef _DEBUG	
	size_t BlockFinder::GetMustBeBifurcation(StrandIterator jt, size_t k)
	{		
		std::string buf(std::string(jt, AdvanceForward(jt, k)));
		KMerBifMap::iterator kt = idMap.find(buf);
		return kt == idMap.end() ? BifurcationStorage::NO_BIFURCATION : kt->second;
	}

	void BlockFinder::Test(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k)
	{	
		for(size_t strand = 0; strand < 2; strand++)
		{
			for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
			{
				StrandIterator begin = sequence.Begin((DNASequence::Direction)strand, chr);
				StrandIterator end = sequence.End((DNASequence::Direction)strand, chr);
				SlidingWindow<StrandIterator> window(begin, end, k);
				for(; window.Valid(); window.Move())
				{
					StrandIterator jt = window.GetBegin();
					size_t pos = sequence.GlobalIndex(jt);
					size_t actualBifurcation = bifStorage.GetBifurcation(jt);
					size_t mustBeBifurcation = GetMustBeBifurcation(jt, k);
					assert(actualBifurcation == mustBeBifurcation);
				}	
			}
		}
	}
#endif

	size_t BlockFinder::EnumerateBifurcationsHash(const DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, ProgressCallBack callBack) const
	{
		bifStorage.Clear();
		BifurcationData::BifurcationId bifurcationCount = 0;
		typedef boost::unordered_map<uint64_t, BifurcationData> BifurcationMap;
		BifurcationMap bifurcation(sequence.TotalSize());
		size_t threshold = 0;
		std::vector<Bool> permit(sequence.ChrNumber(), 1);
		for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
		{
			size_t length = std::distance(sequence.PositiveBegin(chr), sequence.PositiveEnd(chr));
			threshold += length;
			if(length <= k)
			{
				permit[chr] = 0;
				continue;
			}

			StrandIterator border[] = 
			{
				ForwardValidKMer(sequence.PositiveBegin(chr), k),
				ForwardValidKMer(sequence.NegativeBegin(chr), k),
				BackwardValidKMer(AdvanceBackward(sequence.PositiveEnd(chr), k), k),
				BackwardValidKMer(AdvanceBackward(sequence.NegativeEnd(chr), k), k)
			};

			KMerHashFunction<DNASequence::StrandIterator> hashF(k);
			for(size_t i = 0; i < 4; i++)
			{
				if(border[i].AtValidPosition())
				{
					uint64_t hash = hashF(border[i]);
					BifurcationMap::iterator jt = bifurcation.find(hash);
					if(jt == bifurcation.end())
					{
						jt = bifurcation.insert(std::make_pair(hash, BifurcationData())).first;
						jt->second.SetId(bifurcationCount++);
					}
				}
			}
		}

		size_t count = 0;
		size_t totalProgress = 0;
		threshold = (threshold * 4) / PROGRESS_STRIDE;

		if(!callBack.empty())
		{
			callBack(totalProgress, start);
		}

		for(size_t strand = 0; strand < 2; strand++)
		{
			for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
			{
				if(!permit[chr])
				{
					continue;
				}

				StrandIterator begin = sequence.Begin((DNASequence::Direction)strand, chr);
				StrandIterator end = sequence.End((DNASequence::Direction)strand, chr);
				SlidingWindow<StrandIterator> window(begin, end, k);
				for(; window.Valid(); window.Move())
				{
					StrandIterator it = window.GetBegin();
					uint64_t hash = window.GetValue();
					if(*it != DNASequence::UNKNOWN_BASE)
					{
						BifurcationMap::iterator jt = bifurcation.find(hash);
						if(jt == bifurcation.end())
						{
							jt = bifurcation.insert(std::make_pair(hash, BifurcationData())).first;
						}

						if(jt != bifurcation.end() && (
							jt->second.Update(window.GetEnd(), BifurcationData::FORWARD) ||
							jt->second.Update(--window.GetBegin(), BifurcationData::BACKWARD)))
						{
							jt->second.SetId(bifurcationCount++);
						}
					}

					if(++count >= threshold && !callBack.empty())
					{
						count = 0;
						totalProgress = std::min(totalProgress + 1, PROGRESS_STRIDE);
						callBack(totalProgress, run);
					}
				}
			}
		}

		for(size_t strand = 0; strand < 2; strand++)
		{
			for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
			{
				if(!permit[chr])
				{
					continue;
				}

				StrandIterator begin = sequence.Begin((DNASequence::Direction)strand, chr);
				StrandIterator end = sequence.End((DNASequence::Direction)strand, chr);
				SlidingWindow<StrandIterator> window(begin, end, k);
				for(; window.Valid(); window.Move())
				{
					BifurcationMap::iterator jt = bifurcation.find(window.GetValue());
					if(jt != bifurcation.end() && jt->second.GetId() != BifurcationStorage::NO_BIFURCATION)
					{
						bifStorage.AddPoint(window.GetBegin(), jt->second.GetId());
					}

					if(++count >= threshold && !callBack.empty())
					{
						count = 0;
						totalProgress = std::min(totalProgress + 1, PROGRESS_STRIDE);
						callBack(totalProgress, run);
					}
				}	
			}
		}

		if(!callBack.empty())
		{
			callBack(PROGRESS_STRIDE, end);
		}

		return bifurcationCount;
	}

	
	size_t BlockFinder::SimplifyGraph(DNASequence & sequence, BifurcationStorage & bifStorage, size_t k, size_t minBranchSize, size_t maxIterations, ProgressCallBack callBack)
	{
		size_t count = 0;
		size_t totalBulges = 0;
		size_t iterations = 0;
		size_t totalProgress = 0;
		bool anyChanges = true;
		if(!callBack.empty())
		{
			callBack(totalProgress, start);
		}

		size_t threshold = (bifStorage.GetMaxId() * maxIterations) / PROGRESS_STRIDE;
		do
		{
			iterations++;
			for(size_t id = 0; id < bifStorage.GetMaxId(); id++)
			{			
				totalBulges += RemoveBulges(sequence, bifStorage, k, minBranchSize, id);
				if(++count >= threshold && !callBack.empty())
				{
					count = 0;
					totalProgress = std::min(totalProgress + 1, PROGRESS_STRIDE);
					callBack(totalProgress, run);
				}
			}
		}
		while((totalBulges > 0) && iterations < maxIterations);

		if(!callBack.empty())
		{
			callBack(PROGRESS_STRIDE, end);
		}
		
		return totalBulges;
	}
}
