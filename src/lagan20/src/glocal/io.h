#ifndef IO
#define IO

#include<stdio.h>
#include<vector>
#include<map>
#include<stdlib.h>
#include<ctype.h>
#include<structs.h>


long long int printChain(Fragment *current);
long long int readInput(char * fileName);
void printAllFragments( long long int numFragments);
void createPointLists(long long int numFragments);
void printPointLists(long long int numFragments);
void printFragment ( Fragment * curfrag );
void findAllNames(long long int numFragments);
void storeIterators(long long int numFragments);
void decideContigBase();

#endif
