#ifndef _SYNTENY_UTILITY_H_
#define _SYNTENY_UTILITY_H_

#include "common.h"

namespace SyntenyBuilder
{
	extern const char DELETED_CHARACTER;
	extern const std::string alphabet_;				
	extern const std::string complement_;	

	class IsInvalid
	{
	public:
		IsInvalid(const std::vector<char> & invalid): invalid_(&invalid) {}
		bool operator () (int position)
		{
			return (*invalid_)[position] == 1;
		}
	private:
		const std::vector<char> * invalid_;
	};
	
	class SetInvalid
	{
	public:
		SetInvalid(std::vector<char> & invalid, int shift): shift_(shift), invalid_(&invalid) {}
		void operator () (int position)
		{
			(*invalid_)[position + shift_] = true;
		}
	private:
		int shift_;
		std::vector<char> * invalid_;
	};
			
	//Structure that contains multiple copies of equivalent edges.
	struct EdgeList
	{		
	public:
		std::vector<int> direct;
		std::vector<int> revComp;
		EdgeList() { }
		EdgeList(const std::vector<int> & direct, const std::vector<int> & revComp):
				direct(direct), revComp(revComp) {}
		size_t Size() const
		{
			return direct.size() + revComp.size();
		}

		void Clear()
		{
			direct.clear();
			revComp.clear();
		}

		size_t Filter(const std::vector<char> & invalid)
		{
			direct.erase(std::remove_if(direct.begin(), direct.end(), IsInvalid(invalid)), direct.end());
			revComp.erase(std::remove_if(revComp.begin(), revComp.end(), IsInvalid(invalid)), revComp.end());
			return Size();
		}	

		void Invalidate(std::vector<char> & invalid, int shift) const
		{
			std::for_each(direct.begin(), direct.end(), SetInvalid(invalid, shift));
			std::for_each(revComp.begin(), revComp.end(), SetInvalid(invalid, shift));
		}
	};

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

	template<class OutputIterator>
		void BuildReverseComplementary(std::string::const_iterator src, size_t size, OutputIterator out)
		{
			size_t count = 0;
			OutputIterator start = out;
			for(std::string::const_iterator iter = src; count < size; iter++)
			{
				if(*iter != DELETED_CHARACTER)
				{
					count++;
					*out++ = complement_[*iter];
				}
			}
			
			std::reverse(start, out);
		}

	template<class OutputIterator>
		void CopySequence(std::string::const_iterator src, size_t size, OutputIterator out)
		{
			size_t count = 0;
			OutputIterator start = out;std::cout << DELETED_CHARACTER;
			for(std::string::const_iterator iter = src; count < size; iter++)
			{
				if(*iter != DELETED_CHARACTER)
				{
					count++;
					*out++ = *iter;
				}
			}
		}
}

#endif 