#include "graphalgorithm.h"

namespace SyntenyBuilder
{
	namespace
	{
	#ifdef _DEBUG
		typedef boost::unordered_map<std::string, size_t> KMerBifMap;
		KMerBifMap idMap;		
	#endif

		struct BifurcationData
		{
		public:
			typedef BifurcationStorage::BifurcationId BifurcationId;
			static const BifurcationId NO_ID;
			static const char NO_CHAR;	

			BifurcationData(BifurcationId id = NO_ID): id_(id), forward_(NO_CHAR), backward_(NO_CHAR) {}
			bool UpdateForward(char nowForward)
			{
				if(id_ == NO_ID && nowForward != NO_CHAR)
				{
					if(forward_ == NO_CHAR)
					{
						forward_ = nowForward;
					}
					else if(forward_ != nowForward)
					{
						return true;
					}
				}

				return false;
			}

			bool UpdateBackward(char nowBackward)
			{
				if(id_ == NO_ID && nowBackward != NO_CHAR)
				{
					if(backward_ == NO_CHAR)
					{
						backward_ = nowBackward;
					}
					else if(backward_ != nowBackward)
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

		const BifurcationData::BifurcationId BifurcationData::NO_ID = -1;
		const char BifurcationData::NO_CHAR = -1;
	}
	
	
#ifdef _DEBUG
	void GraphAlgorithm::Test(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k)
	{
		SlidingWindow<StrandIterator> window[] = 
		{
			SlidingWindow<StrandIterator>(sequence.PositiveBegin(), sequence.PositiveEnd(), k),
			SlidingWindow<StrandIterator>(sequence.NegativeBegin(), sequence.NegativeEnd(), k)
		};
			
		for(size_t strand = 0; strand < 2; strand++)
		{
			for(; window[strand].Valid(); window[strand].Move())
			{
				StrandIterator jt = window[strand].GetBegin();
				size_t actualBifurcation = bifStorage.GetBifurcation(jt);
				KMerBifMap::iterator kt = idMap.find(std::string(jt, AdvanceForward(jt, k)));
				size_t mustbeBifurcation = kt == idMap.end() ? BifurcationStorage::NO_BIFURCATION :
						kt->second;
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
		KMerHashFunction hashF(k);
		SlidingWindow<StrandIterator> window(++sequence.PositiveBegin(), --sequence.PositiveEnd(), k); 
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

		StrandIterator addBorder[] = 
		{
			sequence.PositiveBegin(),
			AdvanceBackward(sequence.PositiveEnd(), k)
		};

		for(size_t i = 0; i < 2; i++)
		{
			size_t hash = hashF(addBorder[i]);
			bifurcation.insert(std::make_pair(hash, BifurcationData(bifurcationCount++)));
			bifurcation[hash].UpdateForward(*AdvanceForward(addBorder[i], k));
		}

		for(size_t count = 0; window.Valid(); window.Move(), count++)
		{
			if(count % MOD == 0)
			{
				std::cerr << "Pos = " << count << std::endl;
			}

			StrandIterator it = window.GetBegin();
			size_t hash = window.GetValue();
			if(*it != 'n')
			{
				BifurcationMap::iterator jt = bifurcation.find(hash);
				if(jt == bifurcation.end())
				{
					jt = bifurcation.insert(std::make_pair(hash, BifurcationData())).first;
				}

				if(jt->second.UpdateForward(*window.GetEnd()) || jt->second.UpdateBackward(*(--it)))
				{
					jt->second.SetId(bifurcationCount++);
				}
			}
		}

		window = SlidingWindow<StrandIterator>(++sequence.NegativeBegin(), --sequence.NegativeEnd(), k); 
		for(size_t count = 0; window.Valid(); window.Move(), count++)
		{
			StrandIterator it = window.GetBegin();
			size_t hash = window.GetValue();
			BifurcationMap::iterator jt = bifurcation.find(hash);

			if(jt != bifurcation.end())
			{
				if(jt->second.UpdateForward(*window.GetEnd()) || jt->second.UpdateBackward(*(--it)))
				{
					jt->second.SetId(bifurcationCount++);
				}
			}			
		}

		for(size_t i = 0; i < 2; i++)
		{
			window = SlidingWindow<StrandIterator>(border[i], border[i + 4], k);
			for(size_t count = 0; window.Valid(); window.Move(), count++)
			{
				if(count % MOD == 0)
				{
					std::cerr << "Pos = " << count << std::endl;
				}

				StrandIterator it = window.GetBegin();
				size_t hash = window.GetValue();
				BifurcationMap::iterator jt = bifurcation.find(hash);
				if(jt != bifurcation.end() && jt->second.GetId() != BifurcationData::NO_ID)
				{
					bifStorage.AddPoint(it, jt->second.GetId());
				}			
			}	
		}

	#ifdef _DEBUG	
		std::cerr << DELIMITER << std::endl << "Bifurcations: " << std::endl;
		for(BifurcationMap::iterator it = bifurcation.begin(); it != bifurcation.end(); ++it)
		{
			if(it->second.GetId() != BifurcationData::NO_ID)
			{
				idMap[std::string(it->first, AdvanceForward(it->first, k))] = it->second.GetId();
				std::cerr << "Id = " << it->second.GetId() << std::endl << "Body = ";
				CopyN(it->first, k, std::ostream_iterator<char>(std::cerr));
				std::cerr << std::endl;				
			}
		}
	#endif

		return bifurcationCount;
	}
/*
	void GraphAlgorithm::FindGraphBulges(const DNASequence & sequence, size_t k)
	{
		size_t totalBulges = 0;
		size_t totalWhirls = 0;
		const size_t MOD = 1000;
		BifurcationStorage bifStorage;
		size_t bifurcationCount = GraphAlgorithm::EnumerateBifurcations(sequence, k, bifStorage);
		std::cerr << "Total bifurcations: " << bifurcationCount << std::endl;
		std::cerr << "Finding bulges..." << std::endl;
		for(size_t id = 0; id < bifurcationCount; id++)
		{
			if(id % MOD == 0)
			{
				std::cout << "id = " << id << std::endl;
			}

			FindBulges(sequence, bifStorage, k, id);
		}
	}
	*/

	void GraphAlgorithm::SimplifyGraph(DNASequence & sequence, 
		BifurcationStorage & bifStorage, size_t k, size_t minBranchSize)
	{
		size_t prevBulges = 0;
		size_t totalBulges = 0;
		size_t totalWhirls = 0;
		size_t iterations = 0;
		const size_t MOD = 1000;
		bool anyChanges = true;
		size_t bifurcationCount = bifStorage.GetMaxId();
		std::cerr << "Total bifurcations: " << bifurcationCount << std::endl;
		do
		{
			totalBulges = 0;
			size_t counter = 0;
			//std::cerr << "Removing whirls..." << std::endl;
			//totalWhirls = RemoveWhirls(bifStorage, sequence, k, minBranchSize);

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

			//std::cerr << "Total whirls: " << totalWhirls << std::endl;
			std::cerr << "Total bulges: " << totalBulges << std::endl;		
		}
		while((totalBulges > 0) && iterations < 4);
	}
}
