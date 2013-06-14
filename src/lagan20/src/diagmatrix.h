#ifndef __DIAGMATRIX_H
#define __DIAGMATRIX_H

#ifdef MULTIAL__FLAG
#include "multial.h"
#else
#include "order.h"
#endif

#define Mmask 0x3
#define Nmask 0x4
#define Omask 0x8
#define NACT 3

typedef struct AlignElement {
  long int M;
  long int N; 
  long int O;
} alel;

typedef struct diagmatrix {
  int d1;
  int d2;
  int* diagindex;  /* this points to where in myelems a certain diagonal starts*/
  int* diagstart;   /* the elem on which the "cross-section" starts*/
  int* diagend;   /* the elem on which the "cross-section" ends */
  int* isneck;   /* if so, give size of next block, 0 ow */
  int numelems;
  int elemsize;
  char* myptrs;
  alel* myelems[NACT];  /* NACT(3) diags active at a time */
  int currdiag;   /*current diagonal */
  int rangelow;
  int currneck;
  align** myneck[2][3]; /* The past 2 necks, 3 ptrs for each */
  int neckdiag[2]; /* For each the size of its 2 diagonals */
} dmat;


dmat* makeDM(int d1, int d2);
void freeDM(dmat* trgt);
void DMinitDiag(dmat* trgt, int* starts, int* ends);
alel* DMgetElem(dmat* trgt, int x, int y);
alel* DMgetElem2(dmat* trgt, int x, int y, alel* prev);
char DMgetPtr(dmat* trgt, int x, int y);
void DMsetPtr(dmat* trgt, char ptr, int x, int y);
align* DMgetNeck(dmat* trgt, int x, int y, int which);
void DMsetNeck(dmat* trgt, align* myal, int x, int y, int which);
alel* DMgetDiagStart(dmat* trgt, int dn, int* size, int* startx, int* starty);
void DMsetElem(dmat* trgt, alel* elem, int x, int y, char ptr);
char DMnextDiag(dmat* trgt);
int DMnextNecks(dmat* trgt, int diag);

#endif
