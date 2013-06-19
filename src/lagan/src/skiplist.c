#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include "skiplist.h" 
#include <time.h>
#include <assert.h>


char init = 0;

void printSLE(sle* tbp) {
  printf("  %d   %x\n", tbp->index, tbp->myelem);
}

int makeLevel() {
  unsigned int r = rand();
  int i = 1;
  while ((r&1) && (i<MAX_LISTS)) {
    i++;
    r = r >> 1;
  }
  /*  printf("lev = %d\n", i);*/
  return i;
}

void initLib() {
  init = 1;
  srand(time(0));
}

/* makes a new skip list*/
sklst* makeSkLst() {
  int i;
  sklst* res = (sklst*) malloc (sizeof(sklst));
  if (!init) {
    fprintf(stderr, "Skip Lists not initialized\n");
    exit(2);
  }
  res->sentinel = mksle(MAX_LISTS, INT_MIN, 0);
  res->maxlevel = 1;
  return res;
}

/*deletes an old skip list */
void delSkLst(sklst* trgt) {
  sle *next, *tbd = trgt->sentinel;
  while(tbd) {
    next = tbd->next[0];
    delSLE(tbd);
    tbd = next;
  }
}

void chklst2(sklst* trgt) {
  sle* tt = trgt->sentinel;
  sle* tt2 = tt->next[0];
  while (tt2) {
    assert(tt->index <= tt2->index);
    assert(tt == tt2->prev[0]);
    tt = tt->next[0];
    tt2 = tt2->next[0];
  }
}

void chklst(sklst* trgt) {
  sle* tt = trgt->sentinel;
  sle* tt2 = tt->next[0];
  while (tt2) {
    assert(tt->index <= tt2->index);
    assert(tt == tt2->prev[0]);
    tt = tt->next[0];
    tt2 = tt2->next[0];
  }
}

sle* SLinsertAfter(sklst* trgt, sle* prev, int index, void* elem) {
  int i;
  sle *tbe;
  int lc = makeLevel();
  if (lc > trgt->maxlevel) {
    trgt->maxlevel = lc;
  }
  tbe = mksle(lc, index, elem);
  for (i = 0; i < tbe->linkcnt; i++) {
    tbe->prev[i] = prev; 
    if (prev->next[i]) {
      prev->next[i]->prev[i] = tbe;
    }
    tbe->next[i] = prev->next[i];
    prev->next[i] = tbe;
    while (prev && i >= prev->linkcnt-1) 
      prev = prev->prev[i];

  }
  return tbe;
}

/*inserts the elem with the index */
sle* SLinsert(sklst* trgt, int index, void* elem) {
  sle* prev = SLfind(trgt, index), *tbe;
  return SLinsertAfter(trgt, prev, index, elem);
}

/*removes & destroys this element */
void SLremove(sklst* trgt, sle* tbr) {
  int i;
  if (trgt)
  for (i = 0; i < tbr->linkcnt; i++) {
    if (tbr->prev[i])
      tbr->prev[i]->next[i] = tbr->next[i];
    if (tbr->next[i])
      tbr->next[i]->prev[i] = tbr->prev[i];
  }
  delSLE(tbr);
}


/* I could just keep a pointer to last, but since I'll rarely 
   use it I'll find it this way instead.. */

sle* SLgetLast(sklst* trgt) {
  int i;
  sle* currpivot = trgt->sentinel;
  i = trgt->maxlevel-1;
  for ( ; i >= 0; i--) {
    while (currpivot->next[i]) {
      currpivot = currpivot->next[i];
    }
  }
  return currpivot;

}

/* Same as the method below, but good for searching for things 
   near the beginning. it uses an up-down method */

sle* SLlowFind(sklst* trgt, int index) {
  int i;
  sle* currpivot = trgt->sentinel;
  i = 0;
  for ( ; i < trgt->maxlevel-1; i++) {
    if (!currpivot->next[i] || currpivot->next[i]->index > index)
      break;
    currpivot = currpivot->next[i];
  }

  for ( ; i >= 0; i--) {

    while (currpivot->index < index) {
      if (!currpivot->next[i]) {
	goto cont;
      }
      currpivot = currpivot->next[i];
    }
    currpivot = currpivot->prev[i];
  cont: {}
  }
  return currpivot;
}

/*gets the elem with the next lowest index. 0 if none */
sle* SLfind(sklst* trgt, int index) {
  int i;
  sle* currpivot = trgt->sentinel;
  i = trgt->maxlevel-1;
  for ( ; i >= 0; i--) {

    while (currpivot->index < index) {
      if (!currpivot->next[i]) {
	goto cont;
      }
      currpivot = currpivot->next[i];
    }
    currpivot = currpivot->prev[i];
  cont: {}
  }
  return currpivot;
  
}

sle* mksle(int linkcnt, int index, void* myelem) {
  int i;
  sle* res = (sle*)malloc (sizeof(sle));
  res->next = (sle**) malloc(linkcnt*sizeof(sle*));
  res->prev = (sle**) malloc(linkcnt*sizeof(sle*));
  res->linkcnt = linkcnt;
  res->index = index;
  res->myelem = myelem;
  for (i = 0; i < linkcnt; i++) {
    res->next[i] = 0;
    res->prev[i] = 0;
  } 
  return res;
}

void delSLE(sle* tbd) {
  free(tbd->next);
  free(tbd->prev);
  free(tbd);
}










