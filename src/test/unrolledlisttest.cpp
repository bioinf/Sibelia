#include "unrolledlisttest.h"

void randomTest();
void memoryTest();

typedef SyntenyFinder::unrolled_list<int, 5> UList5;
typedef SyntenyFinder::unrolled_list<int, 20> UList20;
typedef SyntenyFinder::unrolled_list<int, 100> UList100;

void dummy20 (UList20::iterator, UList20::iterator) {}
void dummy100 (UList100::iterator, UList100::iterator) {}

typedef std::pair<UList5::iterator, int> KeyValue;
typedef std::vector<KeyValue> CustomMap;

bool Cmp(const KeyValue & a, const KeyValue & b)
{
	return a.first == b.first;
}

void UpdateIteratorsBefore(CustomMap & validIterator, std::vector<size_t> & pos, UList5::iterator start, UList5::iterator end)
{
	for(; start != end; ++start)
	{
		pos.push_back(std::find_if(validIterator.begin(), validIterator.end(), boost::bind(Cmp, KeyValue(start, -1), _1)) - validIterator.begin());
	}
}

void UpdateIteratorsAfter(CustomMap & validIterator, std::vector<size_t> & pos, UList5::iterator start, UList5::iterator end)
{
	for(size_t i = 0; start != end; ++start, ++i)
	{
		validIterator[pos[i]].first = start;
	}
}

void UnrolledListConsistencyTest()
{
	UList5 store;
	store.set_erased_value(-1);
	std::vector<int> buf;
	for(int i = 1; i <= 20; i++)
	{
		buf.push_back(i);
	}
	
	store.insert(store.end(), buf.begin(), buf.end());
	int counter = 1;
	//Only valid iterators here!	
	CustomMap validIterator;
	for(UList5::iterator it = store.begin(); it != store.end(); ++it)
	{
		validIterator.push_back(std::make_pair(it, counter++));
	}

	buf.assign(10, 0);
	std::vector<size_t> invalidatedPos;

	//Perform an insertion
	UList5::iterator pos = ++store.begin();	
	store.insert(pos,
		buf.begin(),
		buf.end(),
		boost::bind(UpdateIteratorsBefore, boost::ref(validIterator), boost::ref(invalidatedPos), _1, _2),
		boost::bind(UpdateIteratorsAfter, boost::ref(validIterator), boost::ref(invalidatedPos), _1, _2));	

	for(UList5::iterator it = store.begin(); it != store.end(); ++it)
	{
		std::cout << *it << ' ';
	}

	//Test consistency
	for(size_t i = 0; i < validIterator.size(); i++)
	{
		assert(*validIterator[i].first == validIterator[i].second);
	}
}

void TestUnrolledList()
{
	randomTest();
}

void memoryTest()
{
	UList100 list(-1);
	std::vector<UList100::iterator> inv;
	std::vector<UList100::reverse_iterator> rinv;

	UList100::notify_func func = dummy100;

	for (int i = 0; i < 400 / 4 * 1024; ++i)
    {
        int * tmp = new int[1024];
        list.insert(list.begin(), tmp, tmp + 1024, func, func);
        delete[] tmp;
    }

    std::cin.get();
}

void checkValid(UList20& list, std::list<int>& testList)
{
	UList20::iterator itList = list.begin();
	std::list<int>::iterator itTest = testList.begin();
	while (itTest != testList.end())
	{
		if (*itList != *itTest || itList == list.end())
		{
			for(UList20::iterator it = list.begin(); it != list.end(); ++it) std::cout << *it << " ";
			std::cout << "\n====================\n";
			for(std::list<int>::iterator it = testList.begin(); it != testList.end(); ++it) std::cout << *it << " ";
			std::cout << std::endl;
			assert(0);
		}
		++itList;
		++itTest;
	}
}

void randomTest()
{
	const int TEST_SIZE = 1000;
	const int TEST_UNIT = 50;


	UList20 list(-1);
	std::vector<UList20::iterator> inv;
	std::vector<UList20::reverse_iterator> rinv;

	UList20::notify_func func = dummy20;

	std::list<int> testList;
	for (int i = 0; i < TEST_SIZE; ++i)
    {
    	int random = rand() % 1000;
    	list.push_back(random);
    	testList.push_back(random);
    }

	UList20::iterator itList = list.begin();
	std::list<int>::iterator itTest = testList.begin();
    for (; itList != list.end() && itTest != testList.end(); ++ itList, ++itTest) assert(*itList == *itTest);

	srand((unsigned)time(0));
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

			UList20::iterator itBegin = list.begin(); advance(itBegin, pos);
			itBegin = list.erase(itBegin);
			list.erase(itBegin);

			std::list<int>::iterator testBegin = testList.begin(); advance(testBegin, pos);
			testBegin = testList.erase(testBegin);
			testList.erase(testBegin);

			listSize -= 2;
			action = true;
		}
		checkValid(list, testList);

		//point insert
		if (listSize < TEST_SIZE * 2)
		{
			int pos = rand() % (listSize - 2);

			int value = rand() % 1000;
			UList20::iterator itBegin = list.begin(); advance(itBegin, pos);
			itBegin = list.insert(itBegin, value);
			list.insert(itBegin, value);

			std::list<int>::iterator testBegin = testList.begin(); advance(testBegin, pos);
			testBegin = testList.insert(testBegin, value);
			testList.insert(testBegin, value);

			listSize += 2;
			action = true;
		}
		checkValid(list, testList);

		//forward erase
		if (listSize > TEST_UNIT * 2)
		{
			//how much erase
			int toErase = rand() % TEST_UNIT;
			int pos = rand() % (listSize - toErase - 1);
			UList20::iterator itBegin = list.begin(); advance(itBegin, pos);
			UList20::iterator itEnd = list.begin(); advance(itEnd, pos + toErase);
			std::list<int>::iterator testBegin = testList.begin(); advance(testBegin, pos);
			std::list<int>::iterator testEnd = testList.begin(); advance(testEnd, pos + toErase);

			list.erase(itBegin, itEnd);
			testList.erase(testBegin, testEnd);
			listSize -= toErase;
			action = true;
		}
		checkValid(list, testList);

		//reverse erase
		if (listSize > TEST_UNIT * 2)
		{
			//how much erase
			int toErase = rand() % TEST_UNIT;
			int pos = rand() % (listSize - toErase - 1);
			UList20::reverse_iterator itBegin = list.rbegin(); advance(itBegin, pos);
			UList20::reverse_iterator itEnd = list.rbegin(); advance(itEnd, pos + toErase);
			std::list<int>::reverse_iterator testBegin = testList.rbegin(); advance(testBegin, pos);
			std::list<int>::reverse_iterator testEnd = testList.rbegin(); advance(testEnd, pos + toErase);

			list.erase(itBegin, itEnd);
			testList.erase((++testEnd).base(), (++testBegin).base());
			listSize -= toErase;
			action = true;
		}
		checkValid(list, testList);

		//forward insert
		if (listSize < TEST_SIZE * 2)
		{
			int sizeInsert = rand() % TEST_UNIT + 1;
			int* toInsert = new int[sizeInsert];

			for (int i = 0; i < sizeInsert; ++i)
			{
				toInsert[i] = rand() % 1000;
			}

			int pos = rand() % listSize;
			UList20::iterator itBegin = list.begin(); advance(itBegin, pos);
			std::list<int>::iterator testBegin = testList.begin(); advance(testBegin, pos);

			list.insert(itBegin, toInsert, toInsert + sizeInsert, func, func);
			testList.insert(testBegin, toInsert, toInsert + sizeInsert);
			listSize += sizeInsert;

			delete[] toInsert;
			action = true;
		}
		checkValid(list, testList);

		//reverse insert
		if (listSize < TEST_SIZE * 2)
		{
			int sizeInsert = rand() % TEST_UNIT + 1;
			int* toInsert = new int[sizeInsert];

			for (int i = 0; i < sizeInsert; ++i)
			{
				toInsert[i] = rand() % 1000;
			}

			int pos = rand() % listSize;
			UList20::reverse_iterator itBegin = list.rbegin(); advance(itBegin, pos);
			std::list<int>::reverse_iterator testBegin = testList.rbegin(); advance(testBegin, pos);

			list.insert(itBegin, toInsert, toInsert + sizeInsert, func, func);
			std::reverse(toInsert, toInsert + sizeInsert);
			testList.insert((testBegin).base(), toInsert, toInsert + sizeInsert);
			listSize += sizeInsert;

			delete[] toInsert;
			action = true;
		}

		checkValid(list, testList);
		std::cout 	<< "iteration: " << iteration++ << " list size: " << list.size()
					<< " with " << list.debugGetNodesCount() << " nodes" << std::endl;

		if (!action) break;
	}
}
