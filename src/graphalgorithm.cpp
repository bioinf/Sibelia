#include "graphalgorithm.h"

namespace SyntenyBuilder
{
	namespace
	{
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
	}
	
	
#ifdef _DEBUG
	typedef boost::unordered_map<std::string, size_t> KMerBifMap;
	KMerBifMap idMap;	
	void GraphAlgorithm::Test(const DNASequence & sequence, const BifurcationStorage & bifStorage, size_t k)
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
					size_t actualBifurcation = bifStorage.GetBifurcation(jt);
					std::string buf(std::string(jt, AdvanceForward(jt, k)));
					KMerBifMap::iterator kt = idMap.find(buf);
					size_t mustbeBifurcation = kt == idMap.end() ? BifurcationStorage::NO_BIFURCATION : kt->second;
					assert(actualBifurcation == mustbeBifurcation);
				}	
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
		BifurcationMap bifurcation(sequence.TotalSize());
		for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
		{
			StrandIterator border[] = 
			{
				sequence.PositiveBegin(chr),
				sequence.NegativeBegin(chr),
				AdvanceBackward(sequence.PositiveEnd(chr), k),
				AdvanceBackward(sequence.NegativeEnd(chr), k),	
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
		}

		for(size_t strand = 0; strand < 2; strand++)
		{
			for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
			{
				StrandIterator begin = sequence.Begin((DNASequence::Direction)strand, chr);
				StrandIterator end = sequence.End((DNASequence::Direction)strand, chr);
				SlidingWindow<StrandIterator> window(begin, end, k);
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
				}
			}
		}

		for(size_t strand = 0; strand < 2; strand++)
		{
			for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
			{
				StrandIterator begin = sequence.Begin((DNASequence::Direction)strand, chr);
				StrandIterator end = sequence.End((DNASequence::Direction)strand, chr);
				SlidingWindow<StrandIterator> window(begin, end, k);
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
		}

	#ifdef _DEBUG	
		idMap.clear();
		PrintRaw(sequence, std::cerr);
		std::cerr << DELIMITER << std::endl << "Bifurcations: " << std::endl;
		for(size_t strand = 0; strand < 2; strand++)
		{
			for(size_t chr = 0; chr < sequence.ChrNumber(); chr++)
			{
				StrandIterator begin = sequence.Begin((DNASequence::Direction)strand, chr);
				StrandIterator end = sequence.End((DNASequence::Direction)strand, chr);
				SlidingWindow<StrandIterator> window(begin, end, k);
				for(; window.Valid(); window.Move())
				{
					BifurcationMap::iterator jt = bifurcation.find(window.GetValue());
					std::string buf(window.GetBegin(), AdvanceForward(window.GetBegin(), k));
					if(jt != bifurcation.end() && jt->second.GetId() != BifurcationStorage::NO_BIFURCATION && idMap.count(buf) == 0)
					{
						idMap[buf] = jt->second.GetId();
						std::cerr << "Id = " << jt->second.GetId() << std::endl << "Body = " << buf << std::endl;
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
