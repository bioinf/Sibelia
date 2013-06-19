#ifndef __MULTIAL_H
#define __MULTIAL_H


#include <stdio.h>

#define NUC_FILE "nucmatrix.txt"
#define NUC_FILE_SIZE 6

#define MAX_SEQ 63
#define CNTS_LEN 8
#define CNTS_A 0
#define CNTS_T 1
#define CNTS_C 2
#define CNTS_G 3
#define CNTS_CB 4
#define CNTS_GS 5
#define CNTS_GC 6
#define CNTS_GE 7


typedef struct HitLocationList {
  int seq1start;
  int seq2start;
  int seq1end;
  int seq2end;
  float score;
  struct HitLocationList *next;
  struct HitLocationList *bkptr;
  float scoreSoFar;
  char dirty;
} hll;

typedef struct hllpointer {
  int number;
  char isstart;
  hll* myhll;
} hptr;

typedef struct Sequence {
  char* lets;
  int numlets, numsiglets;
  char* name;
  char* rptr;
  char* filename;
  int leftbound, rightbound;
  int index;
} seq;

typedef struct align_res {
  int num;
  int index;
  int score;
  int algnlen;
  int numseq;
  seq* seqs[MAX_SEQ];
  long long int* algn;
  char* cnts[CNTS_LEN];
  hll* hlls[MAX_SEQ];
  int dirty;
  struct align_res* nextalign;
} align;


seq* mkConsensus(align* ali);
align* mkSimAlign(seq* seq1);
align* makeAlign(align* ali1, align* ali2, hll* anchors, align **uni);
align* removeSeq(align* ali, int seqnum);
void swapHLL(hll* arg);
hll* remapHLLs(hll* anchs, int which, align* aln, int seqnum);
hll* mergeHLLs(hll* anchs1, int wh1, hll* anchs2, int wh2);
hll* getAnchsFromAlign(align* current, int seqnum, int cutoff);
int getSeqNum(align* ali, seq* trgt);
int printTextAlign(FILE *, align* myalign);
int printFASTAAlign(FILE *, align* myalign);
void printSeqsNames(align *a);
void buildcache();

void freeHLLs(hll *myHLL);
void freeSequence(seq *mySequence);
void freeAlign(align *myAlign);

void setScores(int gapperseqV, int overlapV, int glwidthV);

extern char* alpha;

extern int s1start;
extern int s1end;
extern int s2start;
extern int s2end;
//int match;
//int mismatch;
extern int gapstart;
extern int gapend;
extern int gapcont;
extern int gapperseq;
extern int overlap;
extern int glwidth;
extern char dobin;
extern char* nucmatrixfile;

extern float factor, offset;
extern int logs[MAX_SEQ*MAX_SEQ];

extern FILE* outfile;

#endif


















