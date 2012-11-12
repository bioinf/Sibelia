#ifndef UNROLLED_LIST_TEST_H
#define UNROLLED_LIST_TEST_H

#include "../common.h"
#include "../unrolledlist.h"

void UnrolledListConsistencyTest();
void memoryTest();
void randomTest();

inline void testUnrolledList()
{
	UnrolledListConsistencyTest();
	randomTest();
}

#endif