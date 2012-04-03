#ifndef _K_MER_MULTI_SET_H_

#include "common.h"

namespace SyntenyBuilder
{
	//Class for storing k-mers of a sequence. It looks quite monstrous and complicated,
	//but it works. Such structure is needed because there is no google::sparse_hash_multiset
	//and other implementations (such as boost::unordered) of multisets are memory-inefficient.
	class KMerMultiSet
	{
	public:
		KMerMultiSet(int k, std::string::const_iterator sequence, char deletedCharacter);
		void Clear();
		void EraseShift(int shift);
		void InsertShift(int shift);		
		size_t FindEquivalentShift(int shift, std::vector<int> & ret) const;
		size_t FindEquivalentShift(std::string::const_iterator kMer, std::vector<int> & ret) const;

		class KMerEqualTo
		{
		public:
			KMerEqualTo(int k, std::string::const_iterator sequence, char deletedCharacter, const std::string::const_iterator * auxilary = 0):
				k_(k), deletedCharacter_(deletedCharacter), sequence_(sequence), auxilary_(auxilary) {}
			bool operator()(int shift1, int shift2) const;
		private:
			int k_;
			char deletedCharacter_;
			std::string::const_iterator sequence_;
			const std::string::const_iterator * auxilary_;
			bool Mismatch(std::string::const_iterator it1, std::string::const_iterator it2) const;
		};

		class KMerHashFunction
		{
		public:
			KMerHashFunction(int k, std::string::const_iterator sequence, char deletedCharacter, const std::string::const_iterator * auxilary = 0):
				k_(k), deletedCharacter_(deletedCharacter) ,sequence_(sequence), auxilary_(auxilary) {}
			size_t operator()(int shift) const;
		private:
			int k_;
			char deletedCharacter_;
			std::string::const_iterator sequence_;
			const std::string::const_iterator * auxilary_;			
			static const size_t HASH_BASE = 57;

			size_t CalculateHash(std::string::const_iterator it) const;
		};

	private:
		DISALLOW_COPY_AND_ASSIGN(KMerMultiSet);		
		
		static const int DELETED_KEY = -1;
		static const int AUXILARY_KEY = -2;

		int k_;		
		std::string::const_iterator sequence_;
		mutable std::string::const_iterator auxilary_;		
		typedef google::sparse_hash_set<int, KMerHashFunction, KMerEqualTo> ShiftSet; 
		typedef google::sparse_hash_map<int, std::vector<int>, KMerHashFunction, KMerEqualTo> ShiftMultiSet; 
		//This set is dedicated to store K-mers that occur only once (most of k-mer in genome occur rare)
		ShiftSet singleKMer;
		//This set is dedicated to store K-mers that occur many times
		ShiftMultiSet multipleKMer;
	};
}


#define _K_MER_MULTI_SET_H_
#endif _K_MER_MULTI_SET_H_