//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef UNROLLED_LIST_TEST_H
#define UNROLLED_LIST_TEST_H

#include "../common.h"
#include "../unrolledlist.h"

void InsertTest();
void UnrolledListConsistencyTest();
void memoryTest();
void randomTest();

inline void TestUnrolledList()
{
	InsertTest();
	UnrolledListConsistencyTest();
	randomTest();
}

#endif