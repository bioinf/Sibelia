//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _UNROLLED_LIST_H_
#define _UNROLLED_LIST_H_

#include <vector>
#include <list>
#include <iterator>
#include <cassert>
#include <iostream>
#include <boost/function.hpp>

namespace SyntenyFinder
{
	template<class T, class A, size_t NODE_SIZE>
	class unrolled_list
	{
	private:
		struct chunk
		{
			chunk(const T& erased_value);
			T data[NODE_SIZE];
			A meta[NODE_SIZE];
			T erased_value;
			size_t 	count;
			bool	is_end;
		};

	public:
		typedef uint32_t chunk_size;

		class iterator: public std::iterator<std::bidirectional_iterator_tag, T>
		{
		public:
			iterator();
			iterator(const iterator & it);
			A& meta() const;
			T& operator * () const;
			T* operator -> () const;
			iterator& operator ++ ();
			iterator operator ++ (int);
			iterator& operator -- ();
			iterator operator -- (int);
			bool operator == (const iterator & comp) const;
			bool operator != (const iterator & comp) const;
			iterator& operator = (const iterator & to_copy);
			chunk_size& get_padding_int();
			const chunk_size& get_padding_int() const;			
		private:
			friend class unrolled_list;
			typename std::list<chunk>::iterator m_ListPos;
			chunk_size              			m_ArrayPos;
			chunk_size							m_PaddingInt;
		};

		typedef std::reverse_iterator<typename unrolled_list::iterator> reverse_iterator;

	//	class const_iterator: public std::iterator<std::bidirectional_iterator_tag, T>
	//	{
	//	public:
	//		const_iterator();
	//		const_iterator(const const_iterator & it);
	//		const_iterator(const_reverse_iterator);
	//		const T& operator * () const;
	//		const_iterator& operator ++ ();
	//		const_iterator operator ++ (int);
	//		const_iterator& operator -- ();
	//		const_iterator operator -- (int);
	//		bool operator == (const const_iterator & comp) const;
	//		bool operator != (const const_iterator & comp) const;
	//		const_iterator& operator = (const const_iterator & to_copy);
	//	};
	//
	//	class const_reverse_iterator: public std::iterator<std::bidirectional_iterator_tag, T>
	//	{
	//	public:
	//		const_reverse_iterator();
	//		const_reverse_iterator(const const_reverse_iterator & it);
	//		const_reverse_iterator(const_iterator);
	//		const T& operator * () const;
	//		const_reverse_iterator& operator ++ ();
	//		const_reverse_iterator operator ++ (int);
	//		const_reverse_iterator& operator -- ();
	//		const_reverse_iterator operator -- (int);
	//		bool operator == (const const_reverse_iterator & comp) const;
	//		bool operator != (const const_reverse_iterator & comp) const;
	//		const_reverse_iterator& operator = (const const_reverse_iterator & to_copy);
	//	};

		unrolled_list();
		unrolled_list(const T& erased_value);
		unrolled_list(const unrolled_list& other);
		unrolled_list& operator= (const unrolled_list& other);

		iterator begin();
		iterator end();
		//const_iterator begin() const;
		//const_iterator end() const;
		reverse_iterator rbegin();
		reverse_iterator rend();
		//const_reverse_iterator rbegin() const;
		//const_reverse_iterator rend() const;

		typedef boost::function<void (iterator, iterator)> notify_func;

		size_t  size()  const;
		bool    empty() const;

		void 	set_erased_value(const T& erased_value);

		iterator 			erase(iterator start, iterator end);
		iterator 			erase(iterator position);
		reverse_iterator 	erase(reverse_iterator start, reverse_iterator end);
		reverse_iterator 	erase(reverse_iterator position);


		template<class out_it>
		iterator 			insert(iterator target, out_it source_begin, out_it source_end,
									notify_func before = 0, notify_func after = 0);

		template<class out_it>
		reverse_iterator 	insert(reverse_iterator target, out_it source_begin, out_it source_end,
									notify_func before = 0, notify_func after = 0);

		iterator 			insert(iterator pos, const T & value);
		reverse_iterator 	insert(reverse_iterator pos, const T & value);

		void 				push_back(const T & value);

		void				debugPrintList();
		size_t				debugGetNodesCount(){return m_Data.size();}
	private:
		typedef typename std::list<chunk>::iterator type_iter;
		typedef typename std::list<chunk>::reverse_iterator type_rev_iter;

		static iterator 	create_iterator(type_iter listPos, size_t arrayPos);
		void				lazyUpdateBeginEnd();
		void				updateEndMark();

		std::list<chunk>    m_Data;
		size_t              m_Size;
		T					m_ErasedValue;
		bool				m_ErasedValueSet;

		bool				m_BeginEndDirty;
		iterator			m_Begin;
		iterator			m_End;
		type_iter			m_LastChunk;
		bool				m_LastChunkSet;
	};

	template<class T, class A, size_t NODE_SIZE>
	unrolled_list<T, A, NODE_SIZE>::chunk::chunk(const T& _erased_value):
		erased_value(_erased_value),
		count(0),
		is_end(false)
	{
		std::fill(data, data + NODE_SIZE, erased_value);
	}

	////////////////
	//begin iterator
	template<class T, class A, size_t NODE_SIZE>
	unrolled_list<T, A, NODE_SIZE>::iterator::iterator()
	{}

	template<class T, class A, size_t NODE_SIZE>
	unrolled_list<T, A, NODE_SIZE>::iterator::iterator(const iterator & it):
		m_ListPos(it.m_ListPos),
		m_ArrayPos(it.m_ArrayPos),
		m_PaddingInt(it.m_PaddingInt)
	{}

	template<class T, class A, size_t NODE_SIZE>
	T& unrolled_list<T, A, NODE_SIZE>::iterator::operator * () const
	{
		return m_ListPos->data[m_ArrayPos];
	}

	template<class T, class A, size_t NODE_SIZE>
	T* unrolled_list<T, A, NODE_SIZE>::iterator::operator -> () const
	{
		return &m_ListPos->data[m_ArrayPos];
	}

	template<class T, class A, size_t NODE_SIZE>
	typename unrolled_list<T, A, NODE_SIZE>::iterator&
	unrolled_list<T, A, NODE_SIZE>::iterator::operator ++ ()
	{
		iterator prev = *this;
		for (;;)
		{
			if (m_ArrayPos < NODE_SIZE - 1)
			{
				++m_ArrayPos;
			}
			else
			{
				m_ArrayPos = 0;
				if (!m_ListPos->is_end)
				{
					++m_ListPos;
				}
				else
				{
					*this = prev;
					++m_ArrayPos;
					if (m_ArrayPos == NODE_SIZE)
					{
						++m_ListPos;
						m_ArrayPos = 0;
					}
					break;
				}
			}
			if (**this != m_ListPos->erased_value) break;
		}
		return *this;
	}

	template<class T, class A, size_t NODE_SIZE>
	typename unrolled_list<T, A, NODE_SIZE>::iterator
	unrolled_list<T, A, NODE_SIZE>::iterator::operator ++ (int)
	{
		iterator tmp = *this;
		++(*this);
		return tmp;
	}

	template<class T, class A, size_t NODE_SIZE>
	typename unrolled_list<T, A, NODE_SIZE>::chunk_size&
	unrolled_list<T, A, NODE_SIZE>::iterator::get_padding_int ()
	{
		return m_PaddingInt;
	}

	template<class T, class A, size_t NODE_SIZE>
	A& unrolled_list<T, A, NODE_SIZE>::iterator::meta() const
	{
		return m_ListPos->meta[m_ArrayPos];
	}
	
	template<class T, class A, size_t NODE_SIZE>
	const typename unrolled_list<T, A, NODE_SIZE>::chunk_size&
	unrolled_list<T, A, NODE_SIZE>::iterator::get_padding_int () const
	{
		return m_PaddingInt;
	}

	template<class T, class A, size_t NODE_SIZE>
	typename unrolled_list<T, A, NODE_SIZE>::iterator&
	unrolled_list<T, A, NODE_SIZE>::iterator::operator -- ()
	{
		for (;;)
		{
			if (m_ArrayPos > 0)
			{
				--m_ArrayPos;
			}
			else
			{
				m_ArrayPos = NODE_SIZE - 1;
				--m_ListPos;
			}

			if (**this != m_ListPos->erased_value) break;
		}
		return *this;
	}

	template<class T, class A, size_t NODE_SIZE>
	typename unrolled_list<T, A, NODE_SIZE>::iterator
	unrolled_list<T, A, NODE_SIZE>::iterator::operator -- (int)
	{
		iterator tmp = *this;
		--(*this);
		return tmp;
	}

	template<class T, class A, size_t NODE_SIZE>
	bool unrolled_list<T, A, NODE_SIZE>::iterator::operator == (const iterator & comp) const
	{
		return (m_ListPos == comp.m_ListPos) && (m_ArrayPos == comp.m_ArrayPos);
	}

	template<class T, class A, size_t NODE_SIZE>
	bool unrolled_list<T, A, NODE_SIZE>::iterator::operator != (const iterator & comp) const
	{
		return !(*this == comp);
	}

	template<class T, class A, size_t NODE_SIZE>
	typename unrolled_list<T, A, NODE_SIZE>::iterator&
	unrolled_list<T, A, NODE_SIZE>::iterator::operator = (const iterator & toCopy)
	{
		m_ListPos = toCopy.m_ListPos;
		m_ArrayPos = toCopy.m_ArrayPos;
		m_PaddingInt = toCopy.m_PaddingInt;
		return *this;
	}
	//end iterator
	//////////////

	template<class T, class A, size_t NODE_SIZE>
	void unrolled_list<T, A, NODE_SIZE>::debugPrintList()
	{
		for (type_iter itList = m_Data.begin(); itList != m_Data.end(); ++itList)
		{
			std::cout << "|";
			for (size_t pos = 0; pos < NODE_SIZE; ++pos)
			{
				std::cout << itList->data[pos] << " ";
			}
			if (itList->is_end) std::cout << "e ";
			std::cout << "|" << std::endl;
		}
		std::cout << std::endl;
	}

	template<class T, class A, size_t NODE_SIZE>
	typename unrolled_list<T, A, NODE_SIZE>::iterator
	unrolled_list<T, A, NODE_SIZE>::create_iterator(type_iter listPos, size_t arrayPos)
	{
		iterator toReturn;
		toReturn.m_ListPos = listPos;
		toReturn.m_ArrayPos = static_cast<chunk_size>(arrayPos);
		return toReturn;
	}

	template<class T, class A, size_t NODE_SIZE>
	size_t unrolled_list<T, A, NODE_SIZE>::size() const
	{
		return m_Size;
	}

	template<class T, class A, size_t NODE_SIZE>
	bool unrolled_list<T, A, NODE_SIZE>::empty() const
	{
		return m_Size == 0;
	}

	template<class T, class A, size_t NODE_SIZE>
	unrolled_list<T, A, NODE_SIZE>::unrolled_list():
		m_Size(0),
		m_ErasedValueSet(false),
		m_BeginEndDirty(false),
		m_Begin(create_iterator(m_Data.begin(), 0)),
		m_End(create_iterator(m_Data.end(), 0)),
		m_LastChunkSet(false)
	{
	}


	template<class T, class A, size_t NODE_SIZE>
	unrolled_list<T, A, NODE_SIZE>::unrolled_list(const T& erased_value):
		m_Size(0),
		m_ErasedValue(erased_value),
		m_ErasedValueSet(true),
		m_BeginEndDirty(false),
		m_Begin(create_iterator(m_Data.begin(), 0)),
		m_End(create_iterator(m_Data.end(), 0)),
		m_LastChunkSet(false)
	{
	}

	template<class T, class A, size_t NODE_SIZE>
	unrolled_list<T, A, NODE_SIZE>::unrolled_list(const unrolled_list& other):
		m_Size(other.m_Size),
		m_Data(other.m_Data),
		m_ErasedValue(other.m_ErasedValue),
		m_ErasedValueSet(other.m_ErasedValueSet),
		m_BeginEndDirty(other.m_BeginEndDirty),
		m_LastChunk(other.m_LastChunk),
		m_LastChunkSet(other.m_LastChunkSet),
		m_Begin(other.m_Begin),
		m_End(other.m_End)
	{
	}

	template<class T, class A, size_t NODE_SIZE>
	unrolled_list<T, A, NODE_SIZE>& unrolled_list<T, A, NODE_SIZE>::operator=(const unrolled_list& other)
	{
		m_Size = other.m_Size;
		m_Data = other.m_Data;
		m_ErasedValue = other.m_ErasedValue;
		m_ErasedValueSet = other.m_ErasedValueSet;
		m_BeginEndDirty = other.m_BeginEndDirty;
		m_LastChunk = other.m_LastChunk;
		m_LastChunkSet = other.m_LastChunkSet;
		m_Begin = other.m_Begin;
		m_End = other.m_End;
	}

	template<class T, class A, size_t NODE_SIZE>
	void unrolled_list<T, A, NODE_SIZE>::set_erased_value(const T& erased_value)
	{
		m_ErasedValue = erased_value;
		m_ErasedValueSet = true;
	}

	template<class T, class A, size_t NODE_SIZE>
	void unrolled_list<T, A, NODE_SIZE>::lazyUpdateBeginEnd()
	{
		if (!m_BeginEndDirty) return;

		m_Begin = this->create_iterator(m_Data.begin(), 0);
		if (m_Begin.m_ListPos != m_Data.end())
		{
			while (*m_Begin == m_ErasedValue)
			{
				++m_Begin.m_ArrayPos;
				assert(m_Begin.m_ArrayPos < NODE_SIZE);
			}
		}

		m_End = this->create_iterator(m_Data.end(), 0);
		if (!m_Data.empty())
		{
			m_End.m_ArrayPos = NODE_SIZE;
			--m_End.m_ListPos;
			--m_End;
			++m_End;
		}

		m_BeginEndDirty = false;
	}

	template<class T, class A, size_t NODE_SIZE>
	void unrolled_list<T, A, NODE_SIZE>::updateEndMark()
	{
		if (m_LastChunkSet) m_LastChunk->is_end = false;
		type_iter prevToEnd = --m_Data.end();
		prevToEnd->is_end = true;
		m_LastChunk = prevToEnd;
		m_LastChunkSet = true;
	}

	template<class T, class A, size_t NODE_SIZE>
	typename unrolled_list<T, A, NODE_SIZE>::iterator
	unrolled_list<T, A, NODE_SIZE>::begin()
	{
		this->lazyUpdateBeginEnd();
		return m_Begin;
	}

	template<class T, class A, size_t NODE_SIZE>
	typename unrolled_list<T, A, NODE_SIZE>::iterator
	unrolled_list<T, A, NODE_SIZE>::end()
	{
		this->lazyUpdateBeginEnd();
		return m_End;
	}

	template<class T, class A, size_t NODE_SIZE>
	typename unrolled_list<T, A, NODE_SIZE>::reverse_iterator
	unrolled_list<T, A, NODE_SIZE>::rbegin()
	{
		return reverse_iterator(this->end());
	}

	template<class T, class A, size_t NODE_SIZE>
	typename unrolled_list<T, A, NODE_SIZE>::reverse_iterator
	unrolled_list<T, A, NODE_SIZE>::rend()
	{
		return reverse_iterator(this->begin());
	}

	template<class T, class A, size_t NODE_SIZE>
	typename unrolled_list<T, A, NODE_SIZE>::iterator
	unrolled_list<T, A, NODE_SIZE>::erase(iterator start, iterator end)
	{
		assert(m_ErasedValueSet);

		while(start != end)
		{
			*start = m_ErasedValue;
			assert(start.m_ListPos->count > 0);
			assert(start.m_ListPos->count <= NODE_SIZE);
			--start.m_ListPos->count;
			if (start.m_ListPos->count == 0)
			{
				iterator old = start;
				++old;
				m_Data.erase(start.m_ListPos);
				start = old;
			}
			else
			{
				++start;
			}
			--m_Size;
		}

		if (!m_Data.empty())
		{
			(--m_Data.end())->is_end = true;
			m_LastChunk = --m_Data.end();
		}

		m_BeginEndDirty = true;
		return end;
	}

	template<class T, class A, size_t NODE_SIZE>
	typename unrolled_list<T, A, NODE_SIZE>::iterator
	unrolled_list<T, A, NODE_SIZE>::erase(iterator position)
	{
		iterator next = position;
		++next;
		return this->erase(position, next);
	}


	template<class T, class A, size_t NODE_SIZE>
	typename unrolled_list<T, A, NODE_SIZE>::reverse_iterator
	unrolled_list<T, A, NODE_SIZE>::erase(reverse_iterator start, reverse_iterator end)
	{
		assert(m_ErasedValueSet);
		iterator for_start = (++start).base();
		iterator for_end = (++end).base();
		this->erase(for_end, for_start);

		return end;
	}

	template<class T, class A, size_t NODE_SIZE>
	typename unrolled_list<T, A, NODE_SIZE>::reverse_iterator
	unrolled_list<T, A, NODE_SIZE>::erase(reverse_iterator position)
	{
		reverse_iterator next = position;
		++next;
		return this->erase(position, next);
	}

	template<class T, class A, size_t NODE_SIZE>
	typename unrolled_list<T, A, NODE_SIZE>::iterator
	unrolled_list<T, A, NODE_SIZE>::insert(iterator pos, const T & value)
	{
		const T * ptr = &value;
		return this->insert(pos, ptr, ptr + 1);
	}

	template<class T, class A, size_t NODE_SIZE>
	typename unrolled_list<T, A, NODE_SIZE>::reverse_iterator
	unrolled_list<T, A, NODE_SIZE>::insert(reverse_iterator pos, const T & value)
	{
		const T* ptr = &value;
		return this->insert(pos, ptr, ptr + 1);
	}

	template<class T, class A, size_t NODE_SIZE>
	void unrolled_list<T, A, NODE_SIZE>::push_back(const T & value)
	{
		this->insert(this->end(), value);
	}

	template<class T, class A, size_t NODE_SIZE>
	template <class out_it>
	typename unrolled_list<T, A, NODE_SIZE>::iterator
	unrolled_list<T, A, NODE_SIZE>::insert(iterator target, out_it source_begin, out_it source_end,
										notify_func notify_before, notify_func notify_after)
	{
		assert(m_ErasedValueSet);
		assert(source_begin != source_end);

		//look behind
		if (target != this->begin())
		{
			--target;
			++target.m_ArrayPos;
			if (target.m_ArrayPos == NODE_SIZE)
			{
				++target.m_ListPos;
				target.m_ArrayPos = 0;
			}
		}

		type_iter itList = target.m_ListPos;
		size_t arrayPos = target.m_ArrayPos;
		iterator to_return;
		iterator inv_after_begin;
		iterator inv_after_end;
		bool invalidated = false;
		bool return_set = false;

		while (source_begin != source_end)
		{
			//insert new chunk if need
			if (itList == m_Data.end())
			{
				itList = m_Data.insert(itList, chunk(m_ErasedValue));
				this->updateEndMark();
			}

			//cant insert, make new chunks
			if (itList->data[arrayPos] != m_ErasedValue)
			{
				//just insert new chunk
				if (arrayPos == 0)
				{
					itList = m_Data.insert(itList, chunk(m_ErasedValue));
					this->updateEndMark();
				}
				//insert new chunk and move previous elements to it
				else
				{
					//before notification
					if (notify_before)
					{
						iterator inv_before_begin = this->create_iterator(itList, arrayPos);
						iterator inv_before_end = this->create_iterator(itList, 0);

						for (size_t index = NODE_SIZE - 1; index >= arrayPos; --index)
						{
							if (itList->data[index] != m_ErasedValue)
							{
								inv_before_end.m_ArrayPos = static_cast<chunk_size>(index);
								break;
							}
						}
						++inv_before_end;

						m_BeginEndDirty = true;
						notify_before(inv_before_begin, inv_before_end);
					}

					//add new chunk
					type_iter nextNode = itList;
					++nextNode;
					type_iter newChunk = m_Data.insert(nextNode, chunk(m_ErasedValue));
					this->updateEndMark();

					size_t idFrom = arrayPos;
					size_t idTo = 0;

					while (idFrom < NODE_SIZE)
					{
						if (itList->data[idFrom] != m_ErasedValue)
						{
							newChunk->data[idTo] = itList->data[idFrom];
							newChunk->meta[idTo] = itList->meta[idFrom];
							++newChunk->count;
							itList->data[idFrom] = m_ErasedValue;
							--itList->count;
							++idTo;
						}
						++idFrom;
					}

					//for after notifications
					if (notify_after)
					{
						assert(!invalidated);
						inv_after_begin = this->create_iterator(newChunk, 0);
						inv_after_end = this->create_iterator(newChunk, idTo - 1);
						invalidated = true;
					}
				}
			}

			//inserting new element over erased value
			itList->data[arrayPos] = *source_begin;
			++itList->count;
			++source_begin;

			if (!return_set)
			{
				to_return = this->create_iterator(itList, arrayPos);
				return_set = true;
			}

			++m_Size;
			++arrayPos;

			if (arrayPos == NODE_SIZE)
			{
				++itList;
				arrayPos = 0;
			}
		}

		m_BeginEndDirty = true;
		if (invalidated)
		{
			++inv_after_end;
			notify_after(inv_after_begin, inv_after_end);
		}

		return to_return;
	}

	template<class T, class A, size_t NODE_SIZE>
	template <class out_it>
	typename unrolled_list<T, A, NODE_SIZE>::reverse_iterator
	unrolled_list<T, A, NODE_SIZE>::insert(reverse_iterator target, out_it source_begin, out_it source_end,
										notify_func notify_before, notify_func notify_after)
	{
		assert(m_ErasedValueSet);
		assert(source_begin != source_end);

		//look behind
		bool modified = false;
		if (target != this->rbegin())
		{
			--target;
			modified = true;
		}

		type_rev_iter itList = --std::reverse_iterator<type_iter> ((--target.base()).m_ListPos);
		size_t arrayPos = (--target.base()).m_ArrayPos;

		if (modified)
		{
			if (arrayPos > 0)
			{
				--arrayPos;
			}
			else
			{
				++itList;
				arrayPos = NODE_SIZE - 1;
			}
		}

		reverse_iterator to_return;
		iterator inv_after_begin;
		iterator inv_after_end;
		bool invalidated = false;
		bool return_set = false;

		while (source_begin != source_end)
		{
			if (itList == m_Data.rend())
			{
				itList = std::reverse_iterator<type_iter> ( m_Data.insert(itList.base(), chunk(m_ErasedValue)) );
				--itList;
				this->updateEndMark();
			}

			if (itList->data[arrayPos] != m_ErasedValue)
			{
				if (arrayPos == NODE_SIZE - 1)
				{
					itList = std::reverse_iterator<type_iter> ( m_Data.insert(itList.base(), chunk(m_ErasedValue)) );
					this->updateEndMark();
					--itList;
				}
				else
				{
					//before notification
					if (notify_before)
					{
						iterator inv_before_end = this->create_iterator(--itList.base(), arrayPos);
						iterator inv_before_begin = this->create_iterator(--itList.base(), 0);

						for (size_t index = 0; index <= arrayPos; ++index)
						{
							if (itList->data[index] != m_ErasedValue)
							{
								inv_before_begin.m_ArrayPos = static_cast<chunk_size>(index);
								break;
							}
						}
						++inv_before_end;

						m_BeginEndDirty = true;
						notify_before(inv_before_begin, inv_before_end);
					}

					//adding new chunk
					type_iter nextNode = --itList.base();
					++itList;	//prevent invalidation if itList.base() == end()
					type_iter newChunk = m_Data.insert(nextNode, chunk(m_ErasedValue));
					--itList;
					this->updateEndMark();

					size_t idFrom = arrayPos;
					size_t idTo = NODE_SIZE - 1;
					for (;;)
					{
						if (itList->data[idFrom] != m_ErasedValue)
						{
							newChunk->data[idTo] = itList->data[idFrom];
							newChunk->meta[idTo] = itList->meta[idFrom];
							++newChunk->count;
							itList->data[idFrom] = m_ErasedValue;
							--itList->count;
						}
						if (idFrom > 0)
						{
							--idFrom;
							--idTo;
						}
						else
						{
							break;
						}
					}

					//after notifications
					if (notify_after)
					{
						assert(!invalidated);
						inv_after_end = this->create_iterator(newChunk, NODE_SIZE - 1);
						inv_after_begin = this->create_iterator(newChunk, idTo);
						invalidated = true;
					}
				}
			}

			//insert over erased value
			itList->data[arrayPos] = *source_begin;
			++itList->count;
			++source_begin;

			if (!return_set)
			{
				iterator forward = this->create_iterator(--itList.base(), arrayPos);
				to_return = reverse_iterator(++forward);
				return_set = true;
			}

			++m_Size;
			if (arrayPos > 0)
			{
				--arrayPos;
			}
			else
			{
				++itList;
				arrayPos = NODE_SIZE - 1;
			}
		}

		m_BeginEndDirty = true;
		if (invalidated)
		{
			++inv_after_end;
			notify_after(inv_after_begin, inv_after_end);
		}
		return to_return;
	}
}
#endif
