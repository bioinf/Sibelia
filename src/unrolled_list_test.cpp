#include "unrolledlist.h"
#include <iostream>
#include <algorithm>

typedef SyntenyFinder::unrolled_list<int, 20> MyList;

typedef SyntenyFinder::unrolled_list<int, 100> MyList100;

void f (MyList::iterator, MyList::iterator) {}
//bool rf (MyList::reverse_iterator) {return true;}

void f100 (MyList100::iterator, MyList100::iterator) {}
//bool rf100 (MyList100::reverse_iterator) {return true;}

template<class T>
void advance(T& c, int n)
{
	for (int i = 0; i < n; ++i)
	{
		++c;
	}
}

void randomTest();
void memoryTest();

int main()
{
	randomTest();

	std::cin.get();

	return 0;
}

void memoryTest()
{
	MyList100 list(-1);
	std::vector<MyList100::iterator> inv;
	std::vector<MyList100::reverse_iterator> rinv;

	MyList100::notify_func func = f100;
	//MyList100::notify_reverse_predicate rpred = rf100;

	for (int i = 0; i < 400 / 4 * 1024; ++i)
    {
        int * tmp = new int[1024];
        list.insert(list.begin(), tmp, tmp + 1024, func, func);
        delete[] tmp;
    }

    //std::cin.get();
}

void randomTest()
{
	const int TEST_SIZE = 10000;
	const int TEST_UNIT = 1000;


	MyList list(-1);
	std::vector<MyList::iterator> inv;
	std::vector<MyList::reverse_iterator> rinv;

	MyList::notify_func func = f;
	//MyList::notify_reverse_predicate rpred = rf;


	//int* toInsert = new int[TEST_SIZE];
	std::list<int> testList;
	for (int i = 0; i < TEST_SIZE; ++i)
    {
    	int random = rand() % 1000;
    	list.push_back(random);
    	testList.push_back(random);
    	//toInsert[i] = rand() % 1000;
    }
    //list.insert(list.begin(), toInsert, toInsert + TEST_SIZE, func, func);
    //testList.insert(testList.begin(), toInsert, toInsert + TEST_SIZE);
    //delete[] toInsert;

	MyList::iterator itList = list.begin();
	std::list<int>::iterator itTest = testList.begin();
    for (; itList != list.end() && itTest != testList.end(); ++ itList, ++itTest) assert(*itList == *itTest);

	srand(time(0));
    int listSize = TEST_SIZE;
    int iteration = 0;
    bool action = false;
    for (;;)
	{
		action = false;

		//point erase
		if (listSize > TEST_UNIT * 2)
		{
			int pos = rand() % (listSize - 2);

			MyList::iterator itBegin = list.begin(); advance(itBegin, pos);
			itBegin = list.erase(itBegin);
			list.erase(itBegin);

			std::list<int>::iterator testBegin = testList.begin(); advance(testBegin, pos);
			testBegin = testList.erase(testBegin);
			testList.erase(testBegin);

			listSize -= 2;
			action = true;
		}

		//point insert
		if (listSize < TEST_SIZE * 2)
		{
			int pos = rand() % (listSize - 2);

			int value = rand() % 1000;
			MyList::iterator itBegin = list.begin(); advance(itBegin, pos);
			itBegin = list.insert(itBegin, value);
			list.insert(itBegin, value);

			std::list<int>::iterator testBegin = testList.begin(); advance(testBegin, pos);
			testBegin = testList.insert(testBegin, value);
			testList.insert(testBegin, value);

			listSize += 2;
			action = true;
		}

		//forward erase
		if (listSize > TEST_UNIT * 2)
		{
			//how much erase
			int toErase = rand() % TEST_UNIT;
			int pos = rand() % (listSize - toErase - 1);
			MyList::iterator itBegin = list.begin(); advance(itBegin, pos);
			MyList::iterator itEnd = list.begin(); advance(itEnd, pos + toErase);
			std::list<int>::iterator testBegin = testList.begin(); advance(testBegin, pos);
			std::list<int>::iterator testEnd = testList.begin(); advance(testEnd, pos + toErase);

			list.erase(itBegin, itEnd);
			testList.erase(testBegin, testEnd);
			listSize -= toErase;
			action = true;
		}
		//reverse erase
		if (listSize > TEST_UNIT * 2)
		{
			//how much erase
			int toErase = rand() % TEST_UNIT;
			int pos = rand() % (listSize - toErase - 1);
			MyList::reverse_iterator itBegin = list.rbegin(); advance(itBegin, pos);
			MyList::reverse_iterator itEnd = list.rbegin(); advance(itEnd, pos + toErase);
			std::list<int>::reverse_iterator testBegin = testList.rbegin(); advance(testBegin, pos);
			std::list<int>::reverse_iterator testEnd = testList.rbegin(); advance(testEnd, pos + toErase);

			list.erase(itBegin, itEnd);
			testList.erase((++testEnd).base(), (++testBegin).base());
			listSize -= toErase;
			action = true;
		}

		//forward insert
		if (listSize < TEST_SIZE * 2)
		{
			int sizeInsert = rand() % TEST_UNIT;
			int* toInsert = new int[sizeInsert];

			for (int i = 0; i < sizeInsert; ++i)
			{
				toInsert[i] = rand() % 1000;
			}

			int pos = rand() % listSize;
			MyList::iterator itBegin = list.begin(); advance(itBegin, pos);
			std::list<int>::iterator testBegin = testList.begin(); advance(testBegin, pos);

			list.insert(itBegin, toInsert, toInsert + sizeInsert, func, func);
			testList.insert(testBegin, toInsert, toInsert + sizeInsert);
			listSize += sizeInsert;

			delete[] toInsert;
			action = true;
		}

		//reverse insert
		if (listSize < TEST_SIZE * 2)
		{
			int sizeInsert = rand() % TEST_UNIT;
			int* toInsert = new int[sizeInsert];

			for (int i = 0; i < sizeInsert; ++i)
			{
				toInsert[i] = rand() % 1000;
			}

			//for(MyList::iterator it = list.begin(); it != list.end(); ++it) std::cout << *it << " ";
			//std::cout << std::endl;

			int pos = rand() % listSize;
			MyList::reverse_iterator itBegin = list.rbegin(); advance(itBegin, pos);
			std::list<int>::reverse_iterator testBegin = testList.rbegin(); advance(testBegin, pos);

			list.insert(itBegin, toInsert, toInsert + sizeInsert, func, func);
			std::reverse(toInsert, toInsert + sizeInsert);
			testList.insert((testBegin).base(), toInsert, toInsert + sizeInsert);
			listSize += sizeInsert;

			delete[] toInsert;
			action = true;
		}


		//////////////////////////////////
		//checking
		int iter = 0;
		MyList::iterator itList = list.begin();
		std::list<int>::iterator itTest = testList.begin();
		while (itTest != testList.end())
		{
			assert(itList != list.end());
			if (*itList != *itTest)
			{
				for(MyList::iterator it = list.begin(); it != list.end(); ++it) std::cout << *it << " ";
				std::cout << "\n====================\n";
				for(std::list<int>::iterator it = testList.begin(); it != testList.end(); ++it) std::cout << *it << " ";
				std::cout << std::endl;
				assert(0);
			}
			++itList;
			++itTest;
			++iter;
		}
		//////////////////////////////
		std::cout 	<< "iteration: " << iteration++ << " list size: " << list.size() << std::endl;

		if (!action) break;

		inv.clear();
		rinv.clear();

		//if (iteration > 500000) break;
	}
}
