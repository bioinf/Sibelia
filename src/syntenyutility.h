#ifndef _SYNTENY_UTILITY_H_
#define _SYNTENY_UTILITY_H_

#include "common.h"
#include "skipiterator.h"

namespace SyntenyBuilder
{
	extern const char DELETED_CHARACTER;
	extern const std::string alphabet_;				
	extern const std::string complement_;	

	typedef char Bool;
	typedef SkipIterator<std::string::iterator> StringIterator;
	typedef SkipConstIterator<std::string::const_iterator> StringConstIterator;
	typedef SkipIterator<std::string::reverse_iterator> StringReverseIterator;
	typedef SkipConstIterator<std::string::const_reverse_iterator> StringConstReverseIterator;
			
	struct VisitData
	{			
		int classId;			
		int startVertex;
		bool direct;
		int step;
		
		VisitData() {}			
		VisitData(int classId, int startVertex, bool direct, int step):
			classId(classId), startVertex(startVertex), direct(direct), step(step) { }
			
		bool operator == (const VisitData & vis) const
		{
			return classId == vis.classId && direct == vis.direct;
		}

		bool operator != (const VisitData & vis) const
		{
			return !(*this == vis);
		}
	};

	template<class InputIterator, class OutputIterator>
		void BuildReverseComplementary(InputIterator it, size_t size, OutputIterator out)
		{
			size_t count = 0;
			OutputIterator start = out;
			for(; count < size; count++)
			{
				*out++ = complement_[*it++];
			}
			
			std::reverse(start, out);
		}

	class RevCompIterator: public std::iterator<std::forward_iterator_tag, char, void> 
	{
	public:
		RevCompIterator(const StringConstReverseIterator & it): it_(it)
		{
				
		}

		char operator *() const
		{
			return complement_[*it_];
		}

		RevCompIterator operator ++ ()
		{			
			return RevCompIterator(++it_);
		}

		RevCompIterator operator ++ (int)
		{
			StringConstReverseIterator ret(it_++);
			return ret;
		}

	private:
		StringConstReverseIterator it_;
	};
}

#endif 