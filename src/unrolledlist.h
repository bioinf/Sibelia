#include <vector>
#include <list>
#include <iterator>
#include <cassert>
#include <boost/function.hpp>

template<class T, size_t NODE_SIZE>
class unrolled_list
{
private:
	struct chunk
	{
		chunk(const T& erased_value);
		T data[NODE_SIZE];
		size_t count;
	};

public:
	class reverse_iterator;
	class iterator: public std::iterator<std::bidirectional_iterator_tag, T>
	{
	public:
		iterator();
		iterator(const iterator & it);
		T& operator * () const;
		iterator& operator ++ ();
		iterator operator ++ (int);
		iterator& operator -- ();
		iterator operator -- (int);
		bool operator == (const iterator & comp) const;
		bool operator != (const iterator & comp) const;
		iterator& operator = (const iterator & to_copy);
	private:
		friend class unrolled_list;

		typename std::list<chunk>::iterator m_ListPos;
		size_t              m_ArrayPos;
		unrolled_list*		m_ListInstance;
	};

	class reverse_iterator: public std::iterator<std::bidirectional_iterator_tag, T>
	{
	public:
		reverse_iterator();
		reverse_iterator(const reverse_iterator & it);
		explicit reverse_iterator(const iterator & it);
		T& operator * () const;
		reverse_iterator& operator ++ ();
		reverse_iterator operator ++ (int);
		reverse_iterator& operator -- ();
		reverse_iterator operator -- (int);
		bool operator == (const reverse_iterator & comp) const;
		bool operator != (const reverse_iterator & comp) const;
		reverse_iterator& operator = (const reverse_iterator & to_copy);
		iterator base() const;
    private:
    	iterator m_Base;
	};


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

	unrolled_list(const T& erased_value);
	unrolled_list(const unrolled_list&);

	iterator begin();
	iterator end();
	//const_iterator begin() const;
	//const_iterator end() const;
	reverse_iterator rbegin();
	reverse_iterator rend();
	//const_reverse_iterator rbegin() const;
	//const_reverse_iterator rend() const;

	typedef boost::function<bool (iterator)> notify_predicate;
	typedef boost::function<bool (reverse_iterator)> notify_reverse_predicate;

    inline size_t  size()  const   	{return m_Size;}
    inline bool    empty() const  	{return m_Size == 0;}

	void erase(iterator start, iterator end);
	void erase(iterator position);
	void erase(reverse_iterator start, reverse_iterator end);
    void erase(reverse_iterator position);

	template<class out_it>
	void insert(iterator target, out_it source_begin, out_it source_end, notify_predicate pd,
                notify_reverse_predicate, std::vector<iterator> & invalidated,
                std::vector<reverse_iterator> & reverse_invalidated);

	template<class out_it>
	void insert(reverse_iterator target, out_it source_begin, out_it source_end, notify_predicate pd,
                notify_reverse_predicate, std::vector<iterator> & invalidated,
				std::vector<reverse_iterator> & reverse_invalidated);
private:

    typedef typename std::list<chunk>::iterator type_iter;
    typedef typename std::list<chunk>::reverse_iterator type_rev_iter;

	std::list<chunk>    m_Data;
	const T				m_ErasedValue;
	size_t              m_Size;
};

template<class T, size_t NODE_SIZE>
unrolled_list<T, NODE_SIZE>::chunk::chunk(const T& erased_value):
	count(0)
{
	std::fill(data, data + NODE_SIZE, erased_value);
}

////////////////
//begin iterator
template<class T, size_t NODE_SIZE>
unrolled_list<T, NODE_SIZE>::iterator::iterator()
{}

template<class T, size_t NODE_SIZE>
unrolled_list<T, NODE_SIZE>::iterator::iterator(const iterator & it):
    m_ListPos(it.m_ListPos),
    m_ArrayPos(it.m_ArrayPos),
    m_ListInstance(it.m_ListInstance)
{}

template<class T, size_t NODE_SIZE>
T& unrolled_list<T, NODE_SIZE>::iterator::operator * () const
{
	return m_ListPos->data[m_ArrayPos];
}

template<class T, size_t NODE_SIZE>
typename unrolled_list<T, NODE_SIZE>::iterator&
unrolled_list<T, NODE_SIZE>::iterator::operator ++ ()
{
	for (;;)
	{
		if (m_ArrayPos < NODE_SIZE - 1)
		{
		    ++m_ArrayPos;
		}
		else
        {
            m_ArrayPos = 0;
			++m_ListPos;
        }

		if (m_ListPos == m_ListInstance->m_Data.end())
        {
            m_ArrayPos = 0;
            break;
        }
        if (**this != m_ListInstance->m_ErasedValue) break;
	}
	return *this;
}

template<class T, size_t NODE_SIZE>
typename unrolled_list<T, NODE_SIZE>::iterator
unrolled_list<T, NODE_SIZE>::iterator::operator ++ (int)
{
	iterator tmp = *this;
	++(*this);
	return tmp;
}

template<class T, size_t NODE_SIZE>
typename unrolled_list<T, NODE_SIZE>::iterator&
unrolled_list<T, NODE_SIZE>::iterator::operator -- ()
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

		if (**this != m_ListInstance->m_ErasedValue) break;
	}
	return *this;
}

template<class T, size_t NODE_SIZE>
typename unrolled_list<T, NODE_SIZE>::iterator
unrolled_list<T, NODE_SIZE>::iterator::operator -- (int)
{
	iterator tmp = *this;
	--(*this);
	return tmp;
}

template<class T, size_t NODE_SIZE>
bool unrolled_list<T, NODE_SIZE>::iterator::operator == (const iterator & comp) const
{
	return (m_ListPos == comp.m_ListPos) && (m_ArrayPos == comp.m_ArrayPos);
}

template<class T, size_t NODE_SIZE>
bool unrolled_list<T, NODE_SIZE>::iterator::operator != (const iterator & comp) const
{
	return !(*this == comp);
}

template<class T, size_t NODE_SIZE>
typename unrolled_list<T, NODE_SIZE>::iterator&
unrolled_list<T, NODE_SIZE>::iterator::operator = (const iterator & toCopy)
{
	m_ListPos = toCopy.m_ListPos;
	m_ArrayPos = toCopy.m_ArrayPos;
	m_ListInstance = toCopy.m_ListInstance;
	return *this;
}
//end iterator
//////////////

////////////////////////
//begin reverse_iterator
template<class T, size_t NODE_SIZE>
unrolled_list<T, NODE_SIZE>::reverse_iterator::reverse_iterator()
{}

template<class T, size_t NODE_SIZE>
unrolled_list<T, NODE_SIZE>::reverse_iterator::reverse_iterator(const reverse_iterator & it):
	m_Base(it.m_Base)
{}

template<class T, size_t NODE_SIZE>
unrolled_list<T, NODE_SIZE>::reverse_iterator::reverse_iterator(const iterator & it):
	m_Base(it)
{}

template<class T, size_t NODE_SIZE>
T& unrolled_list<T, NODE_SIZE>::reverse_iterator::operator * () const
{
	unrolled_list<T, NODE_SIZE>::iterator temp(m_Base);
	return *(--temp);
}

template<class T, size_t NODE_SIZE>
typename unrolled_list<T, NODE_SIZE>::reverse_iterator&
unrolled_list<T, NODE_SIZE>::reverse_iterator::operator ++ ()
{
	--m_Base;
	return *this;
}

template<class T, size_t NODE_SIZE>
typename unrolled_list<T, NODE_SIZE>::reverse_iterator
unrolled_list<T, NODE_SIZE>::reverse_iterator::operator ++ (int)
{
	reverse_iterator tmp = *this;
	++(*this);
	return tmp;
}

template<class T, size_t NODE_SIZE>
typename unrolled_list<T, NODE_SIZE>::reverse_iterator&
unrolled_list<T, NODE_SIZE>::reverse_iterator::operator -- ()
{
	++m_Base;
	return *this;
}

template<class T, size_t NODE_SIZE>
typename unrolled_list<T, NODE_SIZE>::reverse_iterator
unrolled_list<T, NODE_SIZE>::reverse_iterator::operator -- (int)
{
	reverse_iterator tmp = *this;
	--(*this);
	return tmp;
}

template<class T, size_t NODE_SIZE>
bool unrolled_list<T, NODE_SIZE>::reverse_iterator::operator == (const reverse_iterator & comp) const
{
	return m_Base == comp.m_Base;
}

template<class T, size_t NODE_SIZE>
bool unrolled_list<T, NODE_SIZE>::reverse_iterator::operator != (const reverse_iterator & comp) const
{
	return !(*this == comp);
}

template<class T, size_t NODE_SIZE>
typename unrolled_list<T, NODE_SIZE>::reverse_iterator&
unrolled_list<T, NODE_SIZE>::reverse_iterator::operator = (const reverse_iterator & toCopy)
{
	m_Base = toCopy.m_Base;
	return *this;
}

template<class T, size_t NODE_SIZE>
typename unrolled_list<T, NODE_SIZE>::iterator
unrolled_list<T, NODE_SIZE>::reverse_iterator::base() const
{
	return m_Base;
}

//end reverse_iterator
//////////////////////

template<class T, size_t NODE_SIZE>
typename unrolled_list<T, NODE_SIZE>::iterator
unrolled_list<T, NODE_SIZE>::begin()
{
    iterator toReturn;
    toReturn.m_ListPos = m_Data.begin();
    toReturn.m_ArrayPos = 0;
    toReturn.m_ListInstance = this;

    if (toReturn.m_ListPos != m_Data.end())
    {
        while (*toReturn == m_ErasedValue)
        {
            ++toReturn.m_ArrayPos;
            assert(toReturn.m_ArrayPos < NODE_SIZE);
        }
    }
    return toReturn;
}

template<class T, size_t NODE_SIZE>
typename unrolled_list<T, NODE_SIZE>::iterator
unrolled_list<T, NODE_SIZE>::end()
{
    iterator toReturn;
    toReturn.m_ListPos = m_Data.end();
    toReturn.m_ArrayPos = 0;
    toReturn.m_ListInstance = this;

    return toReturn;
}

template<class T, size_t NODE_SIZE>
typename unrolled_list<T, NODE_SIZE>::reverse_iterator
unrolled_list<T, NODE_SIZE>::rbegin()
{
    return reverse_iterator(this->end());
}

template<class T, size_t NODE_SIZE>
typename unrolled_list<T, NODE_SIZE>::reverse_iterator
unrolled_list<T, NODE_SIZE>::rend()
{
    return reverse_iterator(this->begin());
}

template<class T, size_t NODE_SIZE>
unrolled_list<T, NODE_SIZE>::unrolled_list(const T& erased_value):
	m_ErasedValue(erased_value),
	m_Size(0)
{
}

template<class T, size_t NODE_SIZE>
unrolled_list<T, NODE_SIZE>::unrolled_list(const unrolled_list& other):
	m_Size(other.m_Size),
	m_Data(other.m_Data),
	m_ErasedValue(other.m_ErasedValue)
{
}

template<class T, size_t NODE_SIZE>
void unrolled_list<T, NODE_SIZE>::erase(iterator start, iterator end)
{
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
}

template<class T, size_t NODE_SIZE>
void unrolled_list<T, NODE_SIZE>::erase(iterator position)
{
    *position = m_ErasedValue;
    assert(position.m_ListPos->count > 0);
    assert(position.m_ListPos->count <= NODE_SIZE);
    --position.m_ListPos->count;
    if (position.m_ListPos->count == 0)
    {
        m_Data.erase(position.m_ListPos);
    }
    --m_Size;
}


template<class T, size_t NODE_SIZE>
void unrolled_list<T, NODE_SIZE>::erase(reverse_iterator start, reverse_iterator end)
{
	iterator for_start((++start).base());
	iterator for_end((++end).base());
	this->erase(for_end, for_start);
}

template<class T, size_t NODE_SIZE>
void unrolled_list<T, NODE_SIZE>::erase(reverse_iterator position)
{
	iterator for_pos ((++position).base());
	this->erase(for_pos);
}

template<class T, size_t NODE_SIZE>
template <class out_it>
void unrolled_list<T, NODE_SIZE>::insert(iterator target, out_it source_begin, out_it source_end,
										notify_predicate pd, notify_reverse_predicate rpd,
										std::vector<iterator> & invalidated,
										std::vector<reverse_iterator> & reverse_invalidated)
{
	type_iter itList = target.m_ListPos;
	size_t arrayPos = target.m_ArrayPos;
	bool once = true;
	for (;;)
	{
		if (source_begin == source_end) return;
		if (itList == m_Data.end())
		{
			itList = m_Data.insert(itList, chunk(m_ErasedValue));
		}
		if (itList->data[arrayPos] == m_ErasedValue)
		{
			itList->data[arrayPos] = *source_begin;
			++itList->count;
			++source_begin;
		}
		else
		{
			if (arrayPos == 0)
			{
				itList = m_Data.insert(itList, chunk(m_ErasedValue));
				itList->data[arrayPos] = *source_begin;
				++itList->count;
				++source_begin;
			}
			else
			{
			    assert(once);
			    once = false;
			    type_iter nextNode = itList;
			    ++nextNode;
				type_iter newChunk = m_Data.insert(nextNode, chunk(m_ErasedValue));
				size_t idFrom = arrayPos;
				size_t idTo = 0;
				while (idFrom < NODE_SIZE)
				{
					if (itList->data[idFrom] != m_ErasedValue)
					{
						iterator inv_iter;
						inv_iter.m_ListInstance = this;
						inv_iter.m_ArrayPos = idFrom;
						inv_iter.m_ListPos = itList;
						if (pd(inv_iter)) invalidated.push_back(inv_iter);
						reverse_iterator inv_rev_iter(inv_iter);
						if (rpd(inv_rev_iter)) reverse_invalidated.push_back(inv_rev_iter);

						newChunk->data[idTo] = itList->data[idFrom];
						++newChunk->count;
						itList->data[idFrom] = m_ErasedValue;
						--itList->count;
						++idTo;
					}
                    ++idFrom;
				}
				itList->data[arrayPos] = *source_begin;
				++itList->count;
				++source_begin;
			}
		}
		++m_Size;
		++arrayPos;
		if (arrayPos == NODE_SIZE)
		{
			++itList;
			arrayPos = 0;
		}
	}
}

template<class T, size_t NODE_SIZE>
template <class out_it>
void unrolled_list<T, NODE_SIZE>::insert(reverse_iterator target, out_it source_begin,
										out_it source_end, notify_predicate pd, notify_reverse_predicate rpd,
										std::vector<iterator> & invalidated,
										std::vector<reverse_iterator> & reverse_invalidated)
{
	type_rev_iter itList = --std::reverse_iterator<type_iter> ((--target.base()).m_ListPos);
	size_t arrayPos = (--target.base()).m_ArrayPos;
	bool once = true;
	for (;;)
	{
		if (source_begin == source_end) break;
		if (itList == m_Data.rend())
		{
			itList = std::reverse_iterator<type_iter> ( m_Data.insert(itList.base(), chunk(m_ErasedValue)) );
			--itList;
		}
		if (itList->data[arrayPos] == m_ErasedValue)
		{
			itList->data[arrayPos] = *source_begin;
			++itList->count;
			++source_begin;
		}
		else
		{
			if (arrayPos == NODE_SIZE - 1)
			{
				itList = std::reverse_iterator<type_iter> ( m_Data.insert(itList.base(), chunk(m_ErasedValue)) );
				--itList;
				itList->data[arrayPos] = *source_begin;
				++itList->count;
				++source_begin;
			}
			else
			{
			    assert(once);
			    once = false;
			    type_iter nextNode = --itList.base();
			    ++itList;	//prevent invalidation if itList.base() == end()
				type_iter newChunk = m_Data.insert(nextNode, chunk(m_ErasedValue));
				--itList;
				size_t idFrom = arrayPos;
				size_t idTo = NODE_SIZE - 1;
				for (;;)
				{
					if (itList->data[idFrom] != m_ErasedValue)
					{
						iterator inv_iter;
						inv_iter.m_ListInstance = this;
						inv_iter.m_ArrayPos = idFrom;
						inv_iter.m_ListPos = --itList.base();
						if (pd(inv_iter)) invalidated.push_back(inv_iter);
						reverse_iterator inv_rev_iter(inv_iter);
						if (rpd(inv_rev_iter)) reverse_invalidated.push_back(inv_rev_iter);

						newChunk->data[idTo] = itList->data[idFrom];
						++newChunk->count;
						itList->data[idFrom] = m_ErasedValue;
						--itList->count;
						--idTo;
					}
					if (idFrom > 0)
                    {
                        --idFrom;
                    }
                    else
                    {
                        break;
                    }
				}
				itList->data[arrayPos] = *source_begin;
				++itList->count;
				++source_begin;
			}
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
}

