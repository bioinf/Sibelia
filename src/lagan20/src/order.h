#ifndef ORDER__H
#define ORDER__H

#include "fchaos.h"

typedef struct align_res {
  int score;
  int algnlen;
  char* algn;
  struct align_res *nextalign;
  int nextloc;
  char dirty;
} align;


//align* makeAlign(dmat* mydm, char* seq1, char* seq2);
int printAlign(char* seq1, char* seq2, align* myalign);
void freeAlign(align* t);
int printBinAlign(char* seq1, char* seq2, align* myalign);
int printTextAlign(char* seq1, char* seq2, align* myalign);

#endif







