#include "kmermultiset.h"

namespace SyntenyBuilder
{
	bool KMerMultiSet::KMerEqualTo::Mismatch(std::string::const_iterator it1, std::string::const_iterator it2) const
	{
		for(int i = 0; i < k_; i++)
		{
			for(; *it1 == deletedCharacter_; it1++);
			for(; *it2 == deletedCharacter_; it2++);
			if(*it1++ != *it2++)
			{
				return true;
			}
		}

		return false;
	}

	bool KMerMultiSet::KMerEqualTo::operator()(int shift1, int shift2) const
	{
		if(shift1 == DELETED_KEY || shift2 == DELETED_KEY)
		{
			return shift1 == shift2;
		}

		if(shift1 == AUXILARY_KEY || shift2 == AUXILARY_KEY)	//We work with some "auxilary" sequence
		{
			shift1 = std::max(shift1, shift2);
			return !Mismatch(sequence_ + shift1, *auxilary_);	//This means that we work with some subseqeunce of the "native" sequence_
		}
		else													
		{
			return !Mismatch(sequence_ + shift1, sequence_ + shift2);
		}

		return false;
	}

	size_t KMerMultiSet::KMerHashFunction::CalculateHash(std::string::const_iterator it) const
	{
		size_t base = 1;
		size_t hash = 0;
		for(int i = 0; i < k_; it++)
		{
			if(*it != deletedCharacter_)			
			{
				i++;
				hash += *it * base;
				base *= HASH_BASE;
			}
		}		

		return hash;
	}

	size_t KMerMultiSet::KMerHashFunction::operator()(int shift) const
	{
		if(shift != AUXILARY_KEY)			//This means that we work with some subseqeunce of the "native" sequence_
		{
			return CalculateHash(sequence_ + shift);
		}
		else								//We work with some "auxilary" sequence
		{
			return CalculateHash(*auxilary_); 
		}

		return 0;
	}

	KMerMultiSet::KMerMultiSet(int k, std::string::const_iterator sequence, char deletedCharacter):
		k_(k), sequence_(sequence),
		singleKMer(0, KMerHashFunction(k, sequence, deletedCharacter, &auxilary_), KMerEqualTo(k, sequence, deletedCharacter, &auxilary_)),
		multipleKMer(0, KMerHashFunction(k, sequence, deletedCharacter, &auxilary_), KMerEqualTo(k, sequence, deletedCharacter, &auxilary_))
	{
		singleKMer.set_deleted_key(DELETED_KEY);
		multipleKMer.set_deleted_key(DELETED_KEY);
	}

	void KMerMultiSet::Clear()
	{
		singleKMer.clear();
		multipleKMer.clear();
	}

	void KMerMultiSet::InsertShift(int shift)
	{
		ShiftSet::iterator jt = singleKMer.find(shift);
		if(jt != singleKMer.end())								 //Test if k-mer occured before as single k-mer
		{
			if(shift != *jt)									 //If it is not the same as stored, store it as a multiple k-mer
			{
				multipleKMer.insert(make_pair(*jt, std::vector<int>(1, shift)));
				singleKMer.erase(shift);				
			}
		}
		else				
		{
			ShiftMultiSet::iterator it = multipleKMer.find(shift);
			if(it != multipleKMer.end())						//Test if k-mer occured before as multiple k-mer
			{
				if(shift != it->first)
				{
					it->second.push_back(shift);
				}
			}
			else												//No multiple or single occurence before, store as a single k-mer
			{
				singleKMer.insert(shift);
			}
		}
	}

	size_t KMerMultiSet::FindEquivalentShift(std::string::const_iterator kMer, std::vector<int> & ret) const
	{
		ret.clear();
		auxilary_ = kMer;
		ShiftSet::const_iterator it = singleKMer.find(AUXILARY_KEY);
		if(it != singleKMer.end())
		{
			ret.push_back(*it);
		}
		else
		{
			ShiftMultiSet::const_iterator jt = multipleKMer.find(AUXILARY_KEY);
			if(jt != multipleKMer.end())
			{
				ret.reserve(jt->second.size() + 1);
				ret.assign(jt->second.begin(), jt->second.end());
				ret.push_back(jt->first);
			}
		}

		return ret.size();
	}

	size_t KMerMultiSet::FindEquivalentShift(int shift, std::vector<int> & ret) const
	{
		ret.clear();
		ShiftSet::const_iterator it = singleKMer.find(shift);
		if(it != singleKMer.end())
		{
			ret.push_back(*it);
		}
		else
		{
			ShiftMultiSet::const_iterator jt = multipleKMer.find(shift);
			if(jt != multipleKMer.end())
			{
				ret.reserve(jt->second.size() + 1);
				ret.assign(jt->second.begin(), jt->second.end());
				ret.push_back(jt->first);
			}
		}

		return ret.size();
	}

	void KMerMultiSet::EraseShift(int shift)
	{
		ShiftSet::iterator it = singleKMer.find(shift);
		if(it != singleKMer.end())
		{
			singleKMer.erase(it);
		}
		else
		{
			ShiftMultiSet::iterator jt = multipleKMer.find(shift);
			if(jt != multipleKMer.end())
			{
				if(shift != jt->first)
				{
					jt->second.erase(std::remove(jt->second.begin(), jt->second.end(), shift), jt->second.end());
				}
				else
				{
					std::vector<int> temp(jt->second.begin(), jt->second.end());
					multipleKMer.erase(jt);
					if(temp.size() > 0)
					{
						multipleKMer.insert(make_pair(temp[0], std::vector<int>(temp.begin() + 1, temp.end())));
					}
				}
			}
		}
	}
}