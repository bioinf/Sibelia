template<class T, size_t NODE_SIZE>
	class unrolled_list
	{
	public:
		class Iterator: public std::iterator...
		{
		public:
			iterator();
			iterator(const iterator & it);
			iterator(reverse_iterator);
			T& operator * () const;
			iterator& operator ++ ();
			iterator operator ++ (int);
			iterator& operator -- ();
			iterator operator -- (int);
			bool operator == (const iterator & comp) const;
			bool operator != (const iterator & comp) const;
			iterator& operator = (const iterator & to_copy);
		};

		class const_iterator: public std::iterator...
		{
		public:
			const_iterator();
			const_iterator(const const_iterator & it);
			const_iterator(const_reverse_iterator);
			const T& operator * () const;
			const_iterator& operator ++ ();
			const_iterator operator ++ (int);
			const_iterator& operator -- ();
			const_iterator operator -- (int);
			bool operator == (const const_iterator & comp) const;
			bool operator != (const const_iterator & comp) const;
			const_iterator& operator = (const const_iterator & to_copy);
		};

		class reverse_iterator: public std::iterator...
		{
		public:
			reverse_iterator();
			reverse_iterator(const reverse_iterator & it);
			reverse_iterator(reverse_iterator);
			T& operator * () const;
			reverse_iterator& operator ++ ();
			reverse_iterator operator ++ (int);
			reverse_iterator& operator -- ();
			reverse_iterator operator -- (int);
			bool operator == (const reverse_iterator & comp) const;
			bool operator != (const reverse_iterator & comp) const;
			reverse_iterator& operator = (const reverse_iterator & to_copy);
		};

		class const_reverse_iterator: public std::iterator...
		{
		public:
			const_reverse_iterator();
			const_reverse_iterator(const const_reverse_iterator & it);
			const_reverse_iterator(const_iterator);
			const T& operator * () const;
			const_reverse_iterator& operator ++ ();
			const_reverse_iterator operator ++ (int);
			const_reverse_iterator& operator -- ();
			const_reverse_iterator operator -- (int);
			bool operator == (const const_reverse_iterator & comp) const;
			bool operator != (const const_reverse_iterator & comp) const;
			const_reverse_iterator& operator = (const const_reverse_iterator & to_copy);
		};

		unrolled_list();
		unrolled_list(const unrolled_list&);
		unrolled_list(const T & erased_indicator);
		iterator begin();
		const_iterator begin() const;
		reverse_iterator rbegin();
		const_reverse_iterator rbegin() const;
		iterator end();
		const_iterator end() const;
		reverse_iterator rend();
		const_reverse_iterator rend() const;
		typedef boost::function<bool (iterator)> notify_predicate;
		void erase(iterator start, iterator end);
		void erase(reverse_iterator start, reverse_iterator end);
		template<class out_it>
			void insert(iterator target, out_it source_begin, out_it source_end, notify_predicate p, std::vector<iterator> & invalidated);
		template<class out_it>
			void insert(reverse_iterator target, out_it source_begin, out_it source_end, notify_predicate p, std::vector<iterator> & invalidated);
	};