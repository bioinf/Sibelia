#ifndef STRUCTS
#define STRUCTS

//general defines
#include <stdio.h>
#include <vector>
#include <map>
#include <set>
#include <stdlib.h>
#include <iostream>
#include <limits.h>
#include <list>
#include <string.h>

using namespace std;

#define RIGHT 0
#define LEFT 1
#define UNRELATED 2

#define NEGINF LLONG_MIN

#define UPSTRANDBITS  3
#define DOWNSTRANDBITS 3
#define RELPOSBITS 3


#define UPSTRANDSHIFT 0
#define DOWNSTRANDSHIFT UPSTRANDBITS
#define RELPOSSHIFT UPSTRANDBITS + DOWNSTRANDBITS
#define TOTALSHIFT UPSTRANDBITS + DOWNSTRANDBITS + RELPOSBITS

#define POSITIVE 1
#define NEGATIVE 0
#define CUTOFF 0

#define TRUE 1
#define FALSE 0

#define INF LLONG_MAX
#define MIN LLONG_MIN
#define NAMESIZE 100


struct ltstr {
	bool operator() (const  char* s1, const char* s2) const {
		return strcmp(s1,s2) < 0;
	}
};


typedef map<const char*,long long int ,ltstr> Name;


typedef struct Fragment {
	long long int seq1Start,seq2Start,seq1End,seq2End;
	char strand;
	float score;
	float totalScore;
	struct Fragment *back;
	char deleted;
	char seq1Name[NAMESIZE];
	Name::iterator nameIter;
	char seq2Name[NAMESIZE];
	long long int base;
	long long int getSeq2End(long long int reflectFlag){ return this->seq2End*((reflectFlag == TRUE)?(-1): 1);};
	long long int getSeq2Start(long long int reflectFlag){return this->seq2Start*((reflectFlag == TRUE)?(-1): 1);};
} Fragment;


typedef struct HitLocationList {
	long long int seq1start;
	long long int seq2start;
	long long int seq1end;
	long long int seq2end;
	float score;
	char strand;
	struct HitLocationList *next;
	struct HitLocationList *bkptr;
	float scoreSoFar;
	char seq1Name[NAMESIZE];
	char seq2Name[NAMESIZE];
} hll;



typedef struct Point {
	long long int seq1,seq2;
	Fragment *frag;
} Point;

#endif
