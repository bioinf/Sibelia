#include <vector>
#include <list>
#include <iterator>
#include <cassert>

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
class unrolled_list
{
private:
	struct chunk
	{
		T* data;
		size_t count;
	};

public:
	class reverse_iterator;
	class iterator: public std::iterator<std::bidirectional_iterator_tag, T>
	{
	public:
		iterator();
		iterator(const iterator & it);
		iterator(const reverse_iterator & it);
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
		std::list<chunk>*   m_ListInstance;
	};

	class reverse_iterator: public std::iterator<std::bidirectional_iterator_tag, T>
	{
	public:
		reverse_iterator();
		reverse_iterator(const reverse_iterator & it);
		reverse_iterator(const iterator & it);
		T& operator * () const;
		reverse_iterator& operator ++ ();
		reverse_iterator operator ++ (int);
		reverse_iterator& operator -- ();
		reverse_iterator operator -- (int);
		bool operator == (const reverse_iterator & comp) const;
		bool operator != (const reverse_iterator & comp) const;
		reverse_iterator& operator = (const reverse_iterator & to_copy);
    private:
        friend class unrolled_list;
		typename std::list<chunk>::reverse_iterator m_ListPos;
		size_t              m_ArrayPos;
		std::list<chunk>*   m_ListInstance;
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

	unrolled_list();
//	unrolled_list(const unrolled_list&);
//	unrolled_list(const T & erased_indicator);
	~unrolled_list();

	iterator begin();
	iterator end();
	//const_iterator begin() const;
	//const_iterator end() const;
	reverse_iterator rbegin();
	reverse_iterator rend();
	//const_reverse_iterator rbegin() const;
	//const_reverse_iterator rend() const;

	//typedef boost::function<bool (iterator)> notify_predicate;
	//typedef boost::function<bool (reverse_iterator)> notify_reverse_predicate;
	typedef bool notify_predicate;
	typedef bool notify_reverse_predicate;

    inline size_t  size()   {return m_Size;}
    inline bool    empty()  {return m_Size == 0;}

	void erase(iterator start, iterator end);
	void erase(iterator position);
	//void erase(reverse_iterator start, reverse_iterator end);
    //void erase(reverse_iterator position);

	template<class out_it>
	void insert(iterator target, out_it source_begin, out_it source_end, notify_predicate pd,
                notify_reverse_predicate, std::vector<iterator> & invalidated,
                std::vector<reverse_iterator> & reverse_invalidated);

	template<class out_it>
	void insert(reverse_iterator target, out_it source_begin, out_it source_end, notify_predicate pd,
                notify_reverse_predicate, std::vector<iterator> & invalidated,
                 std::vector<reverse_iterator> & reverse_invalidated);
private:
	chunk new_chunk();

    typedef typename std::list<chunk>::iterator type_iter;
    typedef typename std::list<chunk>::reverse_iterator type_rev_iter;
	std::list<chunk>    m_Data;
	size_t              m_Size;
};

////////////////
//begin iterator
template<class T, size_t NODE_SIZE, T ERASED_VALUE>
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator::iterator()
{}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator::iterator(const iterator & it):
    m_ListPos(it.m_ListPos),
    m_ArrayPos(it.m_ArrayPos),
    m_ListInstance(it.m_ListInstance)
{}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator::iterator(const reverse_iterator & it):
    m_ListPos(it.m_ListPos),
    m_ArrayPos(it.m_ArrayPos),
    m_ListInstance(it.m_ListInstance)
{}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
T& unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator::operator * () const
{
	return m_ListPos->data[m_ArrayPos];
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
typename unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator&
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator::operator ++ ()
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

		if (m_ListPos == m_ListInstance->end())
        {
            m_ArrayPos = 0;
            break;
        }
        if (**this != ERASED_VALUE) break;
	}
	return *this;
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
typename unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator::operator ++ (int)
{
	iterator tmp = *this;
	++(*this);
	return tmp;
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
typename unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator&
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator::operator -- ()
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

		if (**this != ERASED_VALUE) break;
	}
	return *this;
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
typename unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator::operator -- (int)
{
	iterator tmp = *this;
	--(*this);
	return tmp;
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
bool unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator::operator == (const iterator & comp) const
{
	return (m_ListPos == comp.m_ListPos) && (m_ArrayPos == comp.m_ArrayPos);
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
bool unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator::operator != (const iterator & comp) const
{
	return !(*this == comp);
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
typename unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator&
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator::operator = (const iterator & toCopy)
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
template<class T, size_t NODE_SIZE, T ERASED_VALUE>
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator::reverse_iterator()
{}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator::reverse_iterator(const reverse_iterator & it):
    m_ListPos(it.m_ListPos),
    m_ArrayPos(it.m_ArrayPos),
    m_ListInstance(it.m_ListInstance)
{}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator::reverse_iterator(const iterator & it):
    m_ListPos(it.m_ListPos),
    m_ArrayPos(it.m_ArrayPos),
    m_ListInstance(it.m_ListInstance)
{}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
T& unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator::operator * () const
{
	return m_ListPos->data[m_ArrayPos];
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
typename unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator&
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator::operator ++ ()
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
			++m_ListPos;
        }

		if (m_ListPos == m_ListInstance->rend())
        {
            m_ArrayPos = NODE_SIZE - 1;
            break;
        }
        if (**this != ERASED_VALUE) break;
	}
	return *this;
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
typename unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator::operator ++ (int)
{
	reverse_iterator tmp = *this;
	++(*this);
	return tmp;
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
typename unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator&
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator::operator -- ()
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
			--m_ListPos;
        }

		if (**this != ERASED_VALUE) break;
	}
	return *this;
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
typename unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator::operator -- (int)
{
	reverse_iterator tmp = *this;
	--(*this);
	return tmp;
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
bool unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator::operator == (const reverse_iterator & comp) const
{
	return (m_ListPos == comp.m_ListPos) && (m_ArrayPos == comp.m_ArrayPos);
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
bool unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator::operator != (const reverse_iterator & comp) const
{
	return !(*this == comp);
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
typename unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator&
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator::operator = (const reverse_iterator & toCopy)
{
	m_ListPos = toCopy.m_ListPos;
	m_ArrayPos = toCopy.m_ArrayPos;
	m_ListInstance = toCopy.m_ListInstance;
	return *this;
}
//end reverse_iterator
//////////////////////

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
typename unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::begin()
{
    iterator toReturn;
    toReturn.m_ListPos = m_Data.begin();
    toReturn.m_ArrayPos = 0;
    toReturn.m_ListInstance = &m_Data;

    if (toReturn.m_ListPos != m_Data.end())
    {
        while (*toReturn == ERASED_VALUE)
        {
            ++toReturn.m_ArrayPos;
            assert(toReturn.m_ArrayPos < NODE_SIZE);
        }
    }
    return toReturn;
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
typename unrolled_list<T, NODE_SIZE, ERASED_VALUE>::iterator
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::end()
{
    iterator toReturn;
    toReturn.m_ListPos = m_Data.end();
    toReturn.m_ArrayPos = 0;
    toReturn.m_ListInstance = &m_Data;

    return toReturn;
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
typename unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::rbegin()
{
    reverse_iterator toReturn;
    toReturn.m_ListPos = m_Data.rbegin();
    toReturn.m_ArrayPos = NODE_SIZE - 1;
    toReturn.m_ListInstance = &m_Data;

    if (toReturn.m_ListPos != m_Data.rend())
    {
        while (*toReturn == ERASED_VALUE)
        {
            --toReturn.m_ArrayPos;
            assert(toReturn.m_ArrayPos >= 0);
        }
    }
    return toReturn;
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
typename unrolled_list<T, NODE_SIZE, ERASED_VALUE>::reverse_iterator
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::rend()
{
    reverse_iterator toReturn;
    toReturn.m_ListPos = m_Data.rend();
    toReturn.m_ArrayPos = NODE_SIZE - 1;
    toReturn.m_ListInstance = &m_Data;

    return toReturn;
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
typename unrolled_list<T, NODE_SIZE, ERASED_VALUE>::chunk
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::new_chunk()
{
	chunk ch;
	ch.data = new T[NODE_SIZE];
	ch.count = 0;
	std::fill(ch.data, ch.data + NODE_SIZE, ERASED_VALUE);
	return ch;
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::unrolled_list():
	m_Size(0)
{
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
unrolled_list<T, NODE_SIZE, ERASED_VALUE>::~unrolled_list()
{
    for (type_iter itList = m_Data.begin(); itList != m_Data.end(); ++itList)
    {
        delete[] itList->data;
    }
}


//template<class T, size_t NODE_SIZE, T ERASED_VALUE>
//unrolled_list<T, NODE_SIZE, ERASED_VALUE>::unrolled_list(const T & erased_indicator):
//	m_Size(0)
//{
//}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
void unrolled_list<T, NODE_SIZE, ERASED_VALUE>::erase(iterator start, iterator end)
{
	while(start != end)
	{
		*start = ERASED_VALUE;
		assert(start.m_ListPos->count > 0);
		--start.m_ListPos->count;
		if (start.m_ListPos->count == 0)
		{
			iterator old = start;
			++old;
			delete[] start.m_ListPos->data;
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

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
void unrolled_list<T, NODE_SIZE, ERASED_VALUE>::erase(iterator position)
{
    *position = ERASED_VALUE;
    assert(position.m_ListPos->count > 0);
    --position.m_ListPos->count;
    if (position.m_ListPos->count == 0)
    {
        delete[] position.m_ListPos->data;
        m_Data.erase(position.m_ListPos);
    }
    --m_Size;
}

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
template <class out_it>
void unrolled_list<T, NODE_SIZE, ERASED_VALUE>::insert	(iterator target, out_it source_begin, out_it source_end,
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
			itList = m_Data.insert(itList, this->new_chunk());
		}
		if (itList->data[arrayPos] == ERASED_VALUE)
		{
			itList->data[arrayPos] = *source_begin;
			++itList->count;
			++source_begin;
		}
		else
		{
			if (arrayPos == 0)
			{
				itList = m_Data.insert(itList, this->new_chunk());
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
				type_iter newChunk = m_Data.insert(nextNode, this->new_chunk());
				size_t idFrom = arrayPos;
				size_t idTo = 0;
				while (idFrom < NODE_SIZE)
				{
					if (itList->data[idFrom] != ERASED_VALUE)
					{
						newChunk->data[idTo] = itList->data[idFrom];
						++newChunk->count;
						itList->data[idFrom] = ERASED_VALUE;
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

//#include <iostream>

template<class T, size_t NODE_SIZE, T ERASED_VALUE>
template <class out_it>
void unrolled_list<T, NODE_SIZE, ERASED_VALUE>::insert	(reverse_iterator target, out_it source_begin,
                                                        out_it source_end, notify_predicate pd,
                                                        notify_reverse_predicate rpd,
														std::vector<iterator> & invalidated,
														std::vector<reverse_iterator> & reverse_invalidated)
{
    //std::cout << "-------------\n";
	type_rev_iter itList = target.m_ListPos;
	size_t arrayPos = target.m_ArrayPos;
	bool once = true;
	for (;;)
	{
		if (source_begin == source_end) break;
		if (itList == m_Data.rend())
		{
			itList = std::reverse_iterator<type_iter> ( m_Data.insert(itList.base(), this->new_chunk()) );
			--itList;
		}
		if (itList->data[arrayPos] == ERASED_VALUE)
		{
			itList->data[arrayPos] = *source_begin;
			++itList->count;
			++source_begin;
		}
		else
		{
			if (arrayPos == NODE_SIZE - 1)
			{
				itList = std::reverse_iterator<type_iter> ( m_Data.insert(itList.base(), this->new_chunk()) );
				--itList;
				itList->data[arrayPos] = *source_begin;
				++itList->count;
				++source_begin;
			}
			else
			{
			    assert(once);
			    once = false;
			    type_rev_iter nextNode = itList;
			    ++nextNode;
				type_iter newChunk = m_Data.insert(nextNode.base(), this->new_chunk());
				size_t idFrom = arrayPos;
				size_t idTo = NODE_SIZE - 1;
				for (;;)
				{
					if (itList->data[idFrom] != ERASED_VALUE)
					{
						newChunk->data[idTo] = itList->data[idFrom];
						++newChunk->count;
						itList->data[idFrom] = ERASED_VALUE;
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
		//
		//for(iterator it = this->begin(); it != this->end(); ++it) std::cout << *it << " ";
        //std::cout << std::endl;
        //
	}
	//std::cout << "-------------\n";
}
