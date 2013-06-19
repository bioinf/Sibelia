#ifndef __DIAGMATRIX_C
#define __DIAGMATRIX_C

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "diagmatrix.h"

#define MAX2(x,y)   ( (x) >= (y) ? (x) : (y) )
#define MIN2(x,y)   ( (x) <= (y) ? (x) : (y) )

alel dummy;

#ifdef MULTIAL__FLAG
extern int *freed, freedsize, freedcap;
extern align *freedptr;
#endif

dmat* makeDM(int d1, int d2) {
  dmat* trgt = (dmat*)malloc(sizeof(dmat));
  int i;
  trgt->d1 = d1;
  trgt->d2 = d2;
 
  trgt->diagindex = (int*) calloc(d1+d2+1, sizeof(int));
  trgt->diagstart = (int*) calloc(d1+d2+1, sizeof(int));
  trgt->diagend = (int*) calloc(d1+d2+1, sizeof(int));
  trgt->isneck = (int*) calloc(d1+d2+1, sizeof(int));
  for (i=0; i < d1+ d2+1; i++) {
    trgt->diagindex[i] = trgt->diagstart[i] = trgt->diagend[i] = -1;
    trgt->isneck[i] = 0;
  }
  trgt->numelems = 0;
  trgt->currdiag = 0;
  trgt->currneck = 0;
  dummy.M = dummy.N = dummy.O = INT_MIN+(1<<28);
  return trgt;
}


void freeDM(dmat* trgt) {
  
  int olddiag = trgt->neckdiag[trgt->currneck%2];
  int prevsize = (olddiag>0)?trgt->diagend[olddiag]-
    trgt->diagstart[olddiag]+1 + trgt->diagend[olddiag-1]-
    trgt->diagstart[olddiag-1]+1 : 0;
  int i, j;

  //  printf("next neck\n");

  for (i=0; i < prevsize; i++) {
    for (j=0; j<3; j++) {
      //      freeAlign(trgt->myneck[trgt->currneck%2][j][i]);
    }
  }

  for (i=0; i< NACT; i++) {
    free (trgt->myelems[i]);
  }
  free(trgt->myptrs);
  free(trgt->diagindex);
  free(trgt->diagstart);
  free(trgt->diagend);
  free(trgt->isneck);
  free(trgt);
  
}

void DMinitDiag(dmat* trgt, int* starts, int* ends) {
  int i, sav = 0;
  long long int j = 0, ts = 0;
  int k = ends[1]-starts[1]+1, ko=-1, kf;
  int ctr=0, cond=0;

  for (i=1; i < trgt->d1+trgt->d2; i++) {
    trgt->diagindex[i] = j;
    trgt->diagstart[i] = starts[i];
    trgt->diagend[i] = ends[i]; 
    kf = (i == trgt->d1+trgt->d2-1)? -1 : ends[i+1]-starts[i+1]+1;

    j += k;
    cond = (k < kf) || (k <= kf && ctr >= 1000 && k <= 200);
    if ((ko >= k) && cond) {
      ctr = 0;
      //      printf("neck %d\n",i);
      
      if (sav) {
	trgt->isneck[sav] = j;
      } 
      else {
	trgt->myptrs = (char*) calloc (j/2+1, sizeof(char));
      }
      ts += j;
      j = k + ko;
      sav = i;
    }
    ctr++;
    ko = k;
    k = kf;
  }
  trgt->diagindex[i] = j;
  trgt->diagstart[i] = starts[i];
  trgt->diagend[i] = ends[i];
  if (sav) 
    trgt->isneck[sav] = j;
  else
    trgt->myptrs = (char*) calloc (j/2+1, sizeof(char));
  trgt->numelems = j;  
  trgt->currdiag = 0;
  ts += j;
  for (i=0; i < NACT; i++)
    trgt->myelems[i] = 0;
  for (i=0; i < 2; i++) {
    for (j=0; j<3; j++)
    trgt->myneck[i][j] = 0;
    trgt->neckdiag[i] = -1;
  }
  fprintf(stderr,"Total size = %lld * 10^6\n", ts/1000000);
}

alel* DMgetDiagStart(dmat* trgt, int dn, int* size, int* startx, int* starty) {

  alel* res = trgt->myelems[dn%NACT];
  *size = trgt->diagend[dn] - trgt->diagstart[dn]+1;

  if (dn < trgt->d2) {
    *startx = trgt->diagstart[dn]+1;
    *starty = dn - trgt->diagstart[dn];
  }
  else {
    *startx = dn - trgt->d2 + trgt->diagstart[dn]+1;
    *starty = trgt->d2 - trgt->diagstart[dn];
  }
  return res;
}

char DMgetPtr(dmat* trgt, int x, int y) {
  int dn = x+y-1;
  int elem = (dn < trgt->d2)? (x-1): trgt->d2-y;
  int res, loc;
  if (dn <= 0 || dn >= trgt->d1+trgt->d2 ||
      elem < trgt->diagstart[dn] || elem > trgt->diagend[dn]){
      
    return -1;
  }  
  loc = trgt->diagindex[dn] + elem-trgt->diagstart[dn];
  res= trgt->myptrs[loc >> 1];
  if (!(loc & 1))
    res = res >> 4;
  return res & 0xf;
}

void DMsetPtr(dmat* trgt, char ptr, int x, int y) {
  int dn = x+y-1, loc; 
  char res;
  int elem = (dn < trgt->d2)? (x-1): trgt->d2-y;

  if (dn <= 0 || dn >= trgt->d1+trgt->d2 ||
      elem < trgt->diagstart[dn] || elem > trgt->diagend[dn]){      
    fprintf(stderr,"range error!!!\n");
    return;
  }

  dn = trgt->diagindex[dn] + elem-trgt->diagstart[dn];
  if (dn & 1)
    trgt->myptrs[dn >> 1] = (char)(trgt->myptrs[dn >> 1] & 0xf0) | (char)(ptr & 0x0f);
  else
    trgt->myptrs[dn >> 1] = (char)(trgt->myptrs[dn >> 1] & 0x0f) | (char)(ptr << 4);
  
}

alel* DMgetElem(dmat* trgt, int x, int y) {
  register int dn = x+y-1;
  register int elem = (dn < trgt->d2)? (x-1): trgt->d2-y;

  if (dn <= 0 || dn >= trgt->d1+trgt->d2 ||
      elem < trgt->diagstart[dn] || elem > trgt->diagend[dn]){      
    return &dummy;
  }
  return (trgt->myelems[dn % NACT] + elem-trgt->diagstart[dn]);
}

alel* DMgetElem2(dmat* trgt, int x, int y, alel* prev) {
  register int dn = x+y-1;
  register int elem = (dn < trgt->d2)? (x-1): trgt->d2-y;

  if (dn <= 0 || dn >= trgt->d1+trgt->d2 ||
      elem < trgt->diagstart[dn] || elem > trgt->diagend[dn]){      
    return &dummy;
  }

  if (prev != &dummy)
    return prev + 1;
  return (trgt->myelems[dn % NACT] + elem-trgt->diagstart[dn]);
}

void DMsetElem(dmat* trgt, alel* tbi, int x, int y, char ptr) {
  int dn = x+y-1;
  int elem = (dn < trgt->d2)? x: trgt->d2-y;
  if (elem < trgt->diagstart[dn] || elem > trgt->diagend[dn]) {
    fprintf(stderr,"Dummy\n");
    return;
  }
  *(trgt->myelems[dn%NACT]+elem-trgt->diagstart[dn]) = *tbi;
  trgt->myptrs[trgt->diagindex[dn] + elem-trgt->diagstart[dn]]=ptr;
}

char DMnextDiag(dmat* trgt) {
  char* newptrs;
  int i;

  int size = trgt->diagend[trgt->currdiag+1] - trgt->diagstart[trgt->currdiag+1] + 1;
  free(trgt->myelems[(trgt->currdiag+1)%NACT]);
  trgt->myelems[(trgt->currdiag+1)%NACT] = (alel*) calloc(size, sizeof(alel));

  if (trgt->isneck[trgt->currdiag]) {
    //    printf("new pointers!\n");
    newptrs = (char*) calloc ((trgt->isneck[trgt->currdiag]+1)/2+1, sizeof(char)); 
    for (i=0; i< (trgt->isneck[trgt->currdiag]+1)/2+1; i++)
      newptrs[i] = -1;
    free(trgt->myptrs);
    trgt->myptrs = newptrs;
    trgt->diagindex[trgt->currdiag-1] = 0;
    trgt->diagindex[trgt->currdiag] = (trgt->diagend[trgt->currdiag-1] -
				       trgt->diagstart[trgt->currdiag-1] + 1);
  }

  return trgt->isneck[++trgt->currdiag] != 0;
}

int DMnextNecks(dmat* trgt, int diag) {
  int size = trgt->diagend[diag]-trgt->diagstart[diag]+1 +
    trgt->diagend[diag-1]-trgt->diagstart[diag-1]+1;
  
  int olddiag = trgt->neckdiag[trgt->currneck%2];
  int prevsize = (olddiag>0)?trgt->diagend[olddiag]-trgt->diagstart[olddiag]+1 +
    trgt->diagend[olddiag-1]-trgt->diagstart[olddiag-1]+1 : 0;
  int i, j, t1;
  int norm=0;
  int minn = 0;
  //  printf("next neck\n");

  for (i=0; i < prevsize; i++) {
    for (j=0; j<3; j++) {
      if ((trgt->myneck[trgt->currneck%2][j])[i] && 
	  !(trgt->myneck[trgt->currneck%2][j])[i]->dirty){
	freeAlign(trgt->myneck[trgt->currneck%2][j][i]);
	trgt->myneck[trgt->currneck%2][j][i] = 0;
      }
      /*      else if ((trgt->myneck[trgt->currneck%2][j])[i] && 
	       (trgt->myneck[trgt->currneck%2][j])[i]->dirty &&
	       !(trgt->myneck[trgt->currneck%2][j])[i]->nextalign) {
	       fprintf(stderr, "WARN: diag = %d(%d:%d) \n", diag, olddiag, 
	       (trgt->myneck[trgt->currneck%2][j])[i]->algnlen); 
	       }
      */
    }
  }
  for (j=0; j<3; j++) {
    free (trgt->myneck[trgt->currneck%2][j]);
    trgt->myneck[trgt->currneck%2][j] = (align**) calloc (size, sizeof (align*));
    trgt->neckdiag[trgt->currneck%2] = diag;
    for (i=0; i< size; i++) 
      (trgt->myneck[trgt->currneck%2][j])[i] = 0;
  }
  

  size = trgt->diagend[trgt->currdiag] - trgt->diagstart[trgt->currdiag]+1;
  //  fprintf(stderr, "size = %d\n ", size);
  minn  = norm = trgt->myelems[(trgt->currdiag)%NACT][0].M;
  for (j=1; j<size; j++) {
    norm = MAX2 (trgt->myelems[(trgt->currdiag)%NACT][j].M , norm);
    minn = MIN2 (trgt->myelems[(trgt->currdiag)%NACT][j].M , minn);
  } 
  //  fprintf(stderr, "currdiag = %d norm = %d minn = %d\n", trgt->currdiag, norm, minn);
  for (i=0; i < NACT; i++) {
    size = trgt->diagend[trgt->currdiag-i] - trgt->diagstart[trgt->currdiag-i]+1;
    for (j=0; j<size; j++) {
      t1 = trgt->myelems[(trgt->currdiag-i)%NACT][j].M - norm;
      trgt->myelems[(trgt->currdiag-i)%NACT][j].M = (norm > 0)?
	MIN2(trgt->myelems[(trgt->currdiag-i)%NACT][j].M, t1):
	MAX2(trgt->myelems[(trgt->currdiag-i)%NACT][j].M, t1);

      t1 = trgt->myelems[(trgt->currdiag-i)%NACT][j].N - norm;
      trgt->myelems[(trgt->currdiag-i)%NACT][j].N = (norm > 0)?
	MIN2(trgt->myelems[(trgt->currdiag-i)%NACT][j].N, t1):
	MAX2(trgt->myelems[(trgt->currdiag-i)%NACT][j].M, t1);
      t1 = trgt->myelems[(trgt->currdiag-i)%NACT][j].O - norm;
      trgt->myelems[(trgt->currdiag-i)%NACT][j].O = (norm > 0)?
	MIN2(trgt->myelems[(trgt->currdiag-i)%NACT][j].O, t1):
	MAX2(trgt->myelems[(trgt->currdiag-i)%NACT][j].M, t1);
    }
  }

  trgt->currneck++;
  return norm;
}


align* DMgetNeck(dmat* trgt, int x, int y, int which) {
  int dn = x + y - 1;
  int elem = (dn < trgt->d2)? (x-1): trgt->d2-y;
  int fd;

  if (dn <= 0 || dn >= trgt->d1+trgt->d2) {
    return 0;
  }
  if (elem < trgt->diagstart[dn] || elem > trgt->diagend[dn]){
    return 0;
  }
  if (trgt->neckdiag[trgt->currneck%2] == dn) {
    return *(trgt->myneck[trgt->currneck%2][which] + elem-trgt->diagstart[dn]);    
  }
  else if (trgt->neckdiag[trgt->currneck%2] == dn+1) {
    fd = trgt->diagend[dn+1]-trgt->diagstart[dn+1]+1;
    return *(trgt->myneck[trgt->currneck%2][which] + elem-trgt->diagstart[dn] + fd);
  }
  else { fprintf(stderr, "Some dumb error: %d/%d %d %d\n", dn, trgt->d1+trgt->d2-1, trgt->neckdiag[(trgt->currneck-1)%2], trgt->currneck); return 0; }
}

void DMsetNeck(dmat* trgt, align* myal, int x, int y, int which) {
  int dn = x + y - 1;
  int elem = (dn < trgt->d2)? (x-1): trgt->d2-y;
  int fd;

  if (dn <= 0 || dn >= trgt->d1+trgt->d2) {
    fprintf(stderr, "setNeck failed at %d, %d\n", x,y);
    return;
  }
  if (elem < trgt->diagstart[dn] || elem > trgt->diagend[dn]){
    fprintf(stderr, "setNeck failed2 at %d, %d\n", x,y);
    return;
  }
  if (trgt->neckdiag[(trgt->currneck-1)%2] == dn) {
    *(trgt->myneck[(trgt->currneck-1)%2][which] + elem-trgt->diagstart[dn]) = myal;    
  }
  else if (trgt->neckdiag[(trgt->currneck-1)%2] == dn+1) {
    fd = trgt->diagend[dn+1]-trgt->diagstart[dn+1]+1;
    *(trgt->myneck[(trgt->currneck-1)%2][which] + elem-trgt->diagstart[dn] + fd)=myal;
  }
  else { fprintf(stderr, "Some dumb error2: %d %d %d\n", dn, trgt->neckdiag[(trgt->currneck)%2], trgt->currneck); }
}

#endif
