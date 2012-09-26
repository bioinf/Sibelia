#include "graphalgorithm.h"

namespace SyntenyBuilder
{
	namespace
	{
		struct BifurcationData
		{
		public:
			typedef BifurcationStorage::BifurcationId BifurcationId;
			static const char NO_CHAR;	

			BifurcationData(BifurcationId id = BifurcationStorage::NO_BIFURCATION): id_(id), forward_(NO_CHAR), backward_(NO_CHAR) {}
			bool UpdateForward(const DNASequence & sequence, StrandIterator it)
			{
				if(id_ == BifurcationStorage::NO_BIFURCATION && Valid(++it, sequence))
				{
					if(forward_ == NO_CHAR)
					{
						forward_ = *it;
					}
					else if(forward_ != *it)
					{
						return true;
					}
				}

				return false;
			}

			bool UpdateBackward(const DNASequence & sequence, StrandIterator it)
			{
				if(id_ == BifurcationStorage::NO_BIFURCATION && !AtBegin(it, sequence))
				{
					if(backward_ == NO_CHAR)
					{
						backward_ = *--it;
					}
					else if(backward_ != *--it)
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
	}
	
	
#ifdef _DEBUG
	typedef boost::unordered_map<std::string, size_t> KMerBifMap;
	KMerBifMap idMap;	
	void GraphAlgorithm::Test(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k)
	{
		SlidingWindow<StrandIterator> window[] = 
		{
			SlidingWindow<StrandIterator>(sequence.PositiveBegin(), sequence.PositiveEnd(), k),
			SlidingWindow<StrandIterator>(sequence.NegativeBegin(), sequence.NegativeEnd(), k)
		};
			
		for(size_t strand = 0; strand < 2; strand++)
		{
			size_t pos = 0;
			for(; window[strand].Valid(); window[strand].Move(), ++pos)
			{
				StrandIterator jt = window[strand].GetBegin();
				size_t actualBifurcation = bifStorage.GetBifurcation(jt);
				std::string buf(std::string(jt, AdvanceForward(jt, k)));
				KMerBifMap::iterator kt = idMap.find(buf);
				size_t mustbeBifurcation = kt == idMap.end() ? BifurcationStorage::NO_BIFURCATION : kt->second;
				assert(actualBifurcation == mustbeBifurcation);
			}
		}
	}
#endif

	size_t GraphAlgorithm::EnumerateBifurcations(const DNASequence & sequence, BifurcationStorage & bifStorage, size_t k)
	{
		bifStorage.Clear();
		std::cerr << DELIMITER << std::endl;
		std::cerr << "Finding all bifurcations in the graph..." << std::endl;
		
		const size_t MOD = 1000000;
		BifurcationData::BifurcationId bifurcationCount = 0;
		typedef boost::unordered_map<size_t, BifurcationData> BifurcationMap;
		BifurcationMap bifurcation(sequence.Size());

		StrandIterator border[] = 
		{
			sequence.PositiveBegin(),
			sequence.NegativeBegin(),
			AdvanceBackward(sequence.PositiveEnd(), k),
			AdvanceBackward(sequence.NegativeEnd(), k),	
			sequence.PositiveEnd(),
			sequence.NegativeEnd()
		};

		KMerHashFunction hashF(k);
		for(size_t i = 0; i < 4; i++)
		{
			size_t hash = hashF(border[i]);
			BifurcationMap::iterator jt = bifurcation.find(hash);
			if(jt == bifurcation.end())
			{
				jt = bifurcation.insert(std::make_pair(hash, BifurcationData())).first;
				jt->second.SetId(bifurcationCount++);
			}
		}

		for(size_t i = 0; i < 2; i++)
		{
			SlidingWindow<StrandIterator> window(border[i], border[i + 4], k);
			for(size_t count = 0; window.Valid(); window.Move(), count++)
			{
				if(count % MOD == 0)
				{
					std::cerr << "Pos = " << count << std::endl;
				}

				StrandIterator it = window.GetBegin();
				size_t hash = window.GetValue();
				if(*it != DNASequence::UNKNOWN_BASE)
				{
					BifurcationMap::iterator jt = bifurcation.find(hash);
					if(i == 0 && jt == bifurcation.end())
					{
						jt = bifurcation.insert(std::make_pair(hash, BifurcationData())).first;
					}

					if(jt != bifurcation.end() && (jt->second.UpdateForward(sequence, --window.GetEnd()) || jt->second.UpdateBackward(sequence, it)))
					{
						jt->second.SetId(bifurcationCount++);
					}
				}
			}
		}

		for(size_t i = 0; i < 2; i++)
		{
			SlidingWindow<StrandIterator> window = SlidingWindow<StrandIterator>(border[i], border[i + 4], k);
			for(size_t count = 0; window.Valid(); window.Move(), count++)
			{
				if(count % MOD == 0)
				{
					std::cerr << "Pos = " << count << std::endl;
				}

				BifurcationMap::iterator jt = bifurcation.find(window.GetValue());
				if(jt != bifurcation.end() && jt->second.GetId() != BifurcationStorage::NO_BIFURCATION)
				{
					bifStorage.AddPoint(window.GetBegin(), jt->second.GetId());
				}			
			}	
		}				

	#ifdef _DEBUG	
		idMap.clear();
		std::copy(sequence.PositiveBegin(), sequence.PositiveEnd(), std::ostream_iterator<char>(std::cerr));
		std::cerr << std::endl << DELIMITER << std::endl << "Bifurcations: " << std::endl;
		for(size_t i = 0; i < 2; i++)
		{
			for(StrandIterator it = border[i]; it != border[i + 4]; ++it)
			{				
				size_t bifId = bifStorage.GetBifurcation(it);
				if(bifId != BifurcationStorage::NO_BIFURCATION)
				{
					std::string buf(it, AdvanceForward(it, k));
					if(idMap.count(buf) == 0)
					{
						idMap[buf] = bifId;
						std::cerr << "Id = " << bifId << std::endl << "Body = " << buf << std::endl;
					}
				}
			}
		}	

		bifStorage.Dump(sequence, k, std::cerr);
		Test(sequence, bifStorage, k);
	#endif

		return bifurcationCount;
	}

	void GraphAlgorithm::SimplifyGraph(DNASequence & sequence, 
		BifurcationStorage & bifStorage, size_t k, size_t minBranchSize)
	{
		size_t prevBulges = 0;
		size_t totalBulges = 0;
		size_t iterations = 0;
		const size_t MOD = 1000;
		bool anyChanges = true;
		size_t bifurcationCount = bifStorage.GetMaxId();
		std::cerr << "Total bifurcations: " << bifurcationCount << std::endl;
		do
		{
			totalBulges = 0;
			std::cerr << "Iteration: " << iterations++ << std::endl;
			std::cerr << "Removing bulges..." << std::endl;
			for(size_t id = 0; id < bifurcationCount; id++)
			{
				if(id % MOD == 0)
				{
					std::cout << "id = " << id << std::endl;
				}

				totalBulges += RemoveBulges(sequence, bifStorage, k, minBranchSize, id);
			}

			std::cerr << "Total bulges: " << totalBulges << std::endl;		
		}
		while((totalBulges > 0) && iterations < 4);
	}
}
