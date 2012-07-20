#ifndef _SPARSE_MULTI_SET_H_
#define _SPARSE_MULTI_SET_H_

#include "common.h"
#pragma warning(disable:4355)

#undef _DENSE_MAP_

namespace SyntenyBuilder
{
	template<class ActualStore,
			 class Transformer,
			 class StoreHashFunction,
			 class StoreEqualTo>
		class IndexMultiSet
		{
		public:
			IndexMultiSet(size_t size, Transformer transformer, StoreHashFunction hashFunction, StoreEqualTo equalTo):
				transformer_(transformer), hashFunction_(hashFunction), equalTo_(equalTo),
				single_(size, HashFunction(this), EqualTo(this)), multiple_(size, HashFunction(this), EqualTo(this))
			{
			#ifdef _DENSE_MAP_
				single_.set_empty_key(EMPTY_KEY);
				multiple_.set_empty_key(EMPTY_KEY);
			#endif
				single_.set_deleted_key(DELETED_KEY);
				multiple_.set_deleted_key(DELETED_KEY);
			}

			void Clear()
			{
				single_.clear();
				multiple_.clear();
			}

			void Erase(size_t index)
			{
				typename Set::iterator it = single_.find(index);
				if(it != single_.end())
				{
					single_.erase(it);
				}
				else
				{
					typename MultiSet::iterator jt = multiple_.find(index);
					if(jt != multiple_.end())
					{
						if(index != jt->first)
						{
							jt->second.erase(std::remove(jt->second.begin(), jt->second.end(), index), jt->second.end());
						}
						else
						{
							std::vector<size_t> temp(jt->second.begin(), jt->second.end());
							multiple_.erase(jt);
							if(temp.size() > 0)
							{
								multiple_.insert(make_pair(temp[0], std::vector<size_t>(temp.begin() + 1, temp.end())));
							}
						}
					}
				}
			}

			void Insert(size_t index)
			{
				typename Set::iterator jt = single_.find(index);
				if(jt != single_.end())									 //Test if k-mer occured before as single index
				{
					if(index != *jt)									 //If it is not the same as stored, store it as a multiple index
					{
						multiple_.insert(make_pair(*jt, std::vector<size_t>(1, index)));
						single_.erase(index);				
					}
				}
				else				
				{
					typename MultiSet::iterator it = multiple_.find(index);
					if(it != multiple_.end())							//Test if k-mer occured before as multiple index
					{
						if(index != it->first && std::find(it->second.begin(), it->second.end(), index) == it->second.end())
						{
							it->second.push_back(index);
						}
					}
					else												//No multiple or single occurence before, store as a single index
					{
						single_.insert(index);
					}
				}
			}

			size_t Count(size_t index) const
			{
				typename Set::const_iterator it = single_.find(index);
				if(it != single_.end())
				{
					return 1;
				}
				else
				{
					typename MultiSet::const_iterator jt = multiple_.find(index);
					if(jt != multiple_.end())
					{
						return jt->second.size() + 1;
					}
				}

				return 0;
			}

			size_t Count(const ActualStore & auxilary) const
			{
				auxilary_ = &auxilary;
				typename Set::const_iterator it = single_.find(AUXILARY_KEY);
				if(it != single_.end())
				{
					return 1;
				}
				else
				{
					typename MultiSet::const_iterator jt = multiple_.find(AUXILARY_KEY);
					if(jt != multiple_.end())
					{
						return jt->second.size() + 1;
					}
				}

				return 0;
			}

			template<class Iterator>
				size_t Find(size_t index, Iterator ret) const
				{
					size_t count = 0;
					typename Set::const_iterator it = single_.find(index);
					if(it != single_.end())
					{
						*ret++ = *it;
						count = 1;
					}
					else
					{
						typename MultiSet::const_iterator jt = multiple_.find(index);
						if(jt != multiple_.end())
						{
							*ret++ = jt->first;
							std::copy(jt->second.begin(), jt->second.end(), ret);
							count = jt->second.size() + 1;
						}
					}

					return count;
				}
			
			template<class Iterator>
				size_t Find(const ActualStore & auxilary, Iterator ret) const
				{
					size_t count = 0;
					auxilary_ = &auxilary;
					typename Set::const_iterator it = single_.find(AUXILARY_KEY);
					if(it != single_.end())
					{
						*ret++ = *it;
					}
					else
					{
						typename MultiSet::const_iterator jt = multiple_.find(AUXILARY_KEY);
						if(jt != multiple_.end())
						{
							*ret++ = jt->first;;
							std::copy(jt->second.begin(), jt->second.end(), ret);
							count = jt->second.size() + 1;
						}
					}

					return count;
				}

			template<class Iterator>
				void Dump(Iterator ret) const
				{
					for(typename Set::const_iterator it = single_.begin(); it != single_.end(); 
						*ret++ = transformer_(*it++));
					for(typename MultiSet::const_iterator it = multiple_.begin(); it != multiple_.end();
						*ret++ = transformer_(it++->first))
					{
						std::transform(it->second.begin(), it->second.end(), ret, transformer_);
					}
				}

		private:
			DISALLOW_COPY_AND_ASSIGN(IndexMultiSet);
			class HashFunction
			{
			public:
				HashFunction(const IndexMultiSet * set): set_(set) {}
				size_t operator ()(size_t index) const
				{
					if(index == EMPTY_KEY)
					{
						return 0;
					}

					if(index == AUXILARY_KEY)
					{
						return set_->hashFunction_(*set_->auxilary_);
					}

					return set_->hashFunction_(set_->transformer_(index));
				}
			private:
				const IndexMultiSet * set_;
			};

			class EqualTo
			{
			public:
				EqualTo(const IndexMultiSet * set): set_(set) {}
				bool operator () (size_t index1, size_t index2) const
				{
					if(index1 == DELETED_KEY || index2 == DELETED_KEY || index1 == EMPTY_KEY || index2 == EMPTY_KEY)
					{
						return index1 == index2;
					}

					ActualStore item1 = index1 == AUXILARY_KEY ? *set_->auxilary_ : set_->transformer_(index1);
					ActualStore item2 = index2 == AUXILARY_KEY ? *set_->auxilary_ : set_->transformer_(index2);
					return set_->equalTo_(item1, item2);
				}
			private:
				const IndexMultiSet * set_;
			};

		#ifdef _DENSE_MAP_
			typedef google::dense_hash_set<size_t, HashFunction, EqualTo> Set; 
			typedef google::dense_hash_map<size_t, std::vector<size_t>, HashFunction, EqualTo> MultiSet; 
		#else
			typedef google::sparse_hash_set<size_t, HashFunction, EqualTo> Set; 
			typedef google::sparse_hash_map<size_t, std::vector<size_t>, HashFunction, EqualTo> MultiSet; 
		#endif
			size_t size_;
			mutable const ActualStore * auxilary_;
			static const size_t DELETED_KEY;
			static const size_t AUXILARY_KEY;
			static const size_t EMPTY_KEY;
			Transformer transformer_;
			StoreHashFunction hashFunction_;
			StoreEqualTo equalTo_;

			Set single_;
			MultiSet multiple_;
		};

	template<class T1, class T2, class T3, class T4>
		const size_t IndexMultiSet<T1, T2, T3, T4>::DELETED_KEY = -1;	
	template<class T1, class T2, class T3, class T4>
		const size_t IndexMultiSet<T1, T2, T3, T4>::AUXILARY_KEY = -2;	
	template<class T1, class T2, class T3, class T4>
		const size_t IndexMultiSet<T1, T2, T3, T4>::EMPTY_KEY = -3;	
}


#endif
