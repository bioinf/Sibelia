#ifndef __FCHAOS_H
#define __FCHAOS_H

typedef struct GapFreeChunkList {
  int offset;
  int length;
  int score;
  struct GapFreeChunkList *next;
} gfc;

typedef struct HitLocationList {
  int seq1start;
  int seq2start;
  int seq1end;
  int seq2end;
  float score;
  gfc* first;
  gfc* last;
  struct HitLocationList *next;
  char dirty;
} hll;




typedef struct Sequence {
  char* lets;
  int numlets, numsiglets;
  int leftbound, rightbound;
  char* name;
  char* rptr;
} seq;



hll* fchaos(int argc, char** argv);
int mergeOverlap(hll* h1, hll* h2, seq* seq1, seq* seq2);

#endif
