#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "skiplist.h"
#include "thrtrie.h"
#include <assert.h>
int triealphasize=0;
int nnodes=0;


#define DEBUG 1
#define JQ_SIZE 1024
#include "mempage.c"

TJob* jobqueue=0;
int jqsize = 1;
int numjobs = 0;

void makeAlpha(char* alpha) {
  int i;
  int isin = 0;
  for (i=0; i < 256; i++)
    indeces[i] = -1;
  i = 0;
  while (*alpha) {
    if (!isin && *alpha == '[') 
      isin = 1;
    else if (isin && *alpha == ']') {
      isin = 0;
      i++;
    }
    else if (isin) 
      indeces[*alpha] = i;

    else indeces [*alpha] = i++;
    alpha++;
  }
  triealphasize = i;
}

int lookup(char c) {
  return indeces[c];
}


TNode* makeTrie(int height, char* alphabet) {
  TNode* root;
  initMP(0);
  makeAlpha(alphabet);
  if (!jobqueue)
    jobqueue = (TJob*) malloc(sizeof(TJob));
  root = makeNode(height);
  return root;
}

void junker (TNode** m){ 
  
}

int tccc = 0;

void freeTrie (TNode* trgt) {
  /*
  int i;
  if (trgt->height) {
    for (i = 0; i < triealphasize; i++)
      if (trgt->kids.ptrs[i])
	freeTrie(trgt->kids.ptrs[i]);
    junker (trgt->kids.ptrs);
  }
  else
    free(trgt->kids.locator.locs);
  free (trgt);
  */
  MPallfree();
}

TNode* makeNode(int height) {
  TNode* tn = (TNode*) MPmalloc(sizeof(TNode));
  int i;
  tn->height=height;
  if (height) {
    tn->kids.ptrs = (TNode**) MPmalloc(sizeof(TNode*)*triealphasize);
    for (i=0; i < triealphasize; i++) 
      tn->kids.ptrs[i]=0;
  }
  else {
    tn->kids.locator.numlocs = 0;
    tn->kids.locator.locs = (int*)MPmalloc(sizeof(int)*2);
    tn->kids.locator.locssize = 2;
  }
  return tn;
}

int insertLoc (int word, locs* locator) {
  locator->locs[locator->numlocs++] = word;
  if (locator->numlocs >= locator->locssize) {
    locator->locs = (int*) MPrealloc (locator->locs,  sizeof(int)*locator->locssize,
				      sizeof(int)*locator->locssize*2);
    locator->locssize *= 2;
  }
  return 0;
}


int insertWordHelp(TNode* currnode, char* word, char* strbeg, int height,int wordlen) {
  int letter;
  if (height == 0)
    return insertLoc((int)(word-strbeg), &(currnode->kids.locator));
  else {
    letter = lookup(word[wordlen-height]);
    if (letter < 0)
      return 1;
    if (!currnode->kids.ptrs[letter]) {
      currnode->kids.ptrs[letter] = makeNode(height-1);
    }
    return insertWordHelp(currnode->kids.ptrs[letter], word, strbeg, height-1, wordlen);
  }
  return 42;
}

int insertWord(TNode* currnode, char* word, char* strbeg) {
  return insertWordHelp(currnode, word, strbeg, currnode->height, currnode->height);
}

LList* appendLList(LList* a , LList* b) {
  if (!a)
    return b;
  if (!b)
    return a;
  b->next = appendLList(a, b->next);
  return b;
}

/*no longer works */
 /* make iterative??? */
/*
LList* lookupZZZWord(TNode* currnode, char* word, int ndegen) {
  int letter,i;
  LList *temp, *help, *res=0;
  int height = currnode->height;
  if (!currnode || ndegen < 0)
    return 0;
  if (!currnode->height) {
    res = (LList*) malloc (sizeof(LList));
    res->myloc = &currnode->kids.locator;
    res->degleft = 0;
    res->next = 0;
    return res;
  }
  letter = lookup(word[currnode->height-1]);
  if (letter >=0 && currnode->kids.ptrs[letter]) {
    temp = lookupZZZWord(currnode->kids.ptrs[letter], word, ndegen);
    res = appendLList(res, temp);
  }
  for (i=0; i < triealphasize; i++) {
    if (ndegen > 0 && i != letter) {
      if (currnode->kids.ptrs[i]) {
	temp = lookupZZZWord(currnode->kids.ptrs[i], word, ndegen-1);
	help = temp;
	while (help != 0) {
	  help->degloc[help->degleft++] = currnode->height;
	  help = help->next;
	}
	res = appendLList(res, temp);
      }
    }
  }
  return res;
  }*/

void insertString(TNode* root, char* word) {
  char* begin = word;
  int i, j, wordlen = root->height, letprev, letcurr;
  TNode* prev, *curr;
  insertWord(root, word, begin); 
  word++;
  root->backptr = root;
  while (*word) {
    curr = prev = root;
    insertWord(root, word, begin); 
    for (i=0; i < wordlen; i++) {
      letprev = lookup(word[i-1]);
      letcurr = lookup(word[i]);
      if (letprev >= 0)
	prev = prev->kids.ptrs[letprev];
      else break;
      prev->backptr = curr;
      if (letcurr >= 0)
	curr = curr->kids.ptrs[letcurr];
      else break;
    }
    word++;
  }
  letcurr = lookup(*(word-1));
  if (letcurr >=0)
    root->kids.ptrs[letcurr]->backptr = root;
}

void addjob(TNode* tn, char *thisdeg, char dirty, int oldindex) {
  int i;
  jobqueue[numjobs].mynode = tn;
  jobqueue[numjobs].dirty = dirty;
  if (oldindex >= 0) {
    jobqueue[numjobs].numdeg = jobqueue[oldindex].numdeg;
    for (i = 0; i < jobqueue[oldindex].numdeg; i++)
      jobqueue[numjobs].degloc[i] = jobqueue[oldindex].degloc[i];
  }
  else {
    jobqueue[numjobs].numdeg = 0;
  }
  if (thisdeg>0) {
    jobqueue[numjobs].degloc[jobqueue[numjobs].numdeg++] = thisdeg;
  }
  numjobs++;
  if (jqsize == numjobs)
    jobqueue = (TJob*)realloc(jobqueue, sizeof(TJob)*(jqsize *=2));

}

void cleanJobQueue() {
  numjobs = 0;
}


void remjob(int i) {
  jobqueue[i]= jobqueue[--numjobs];
}

LList* makeLList(TJob* tj, char* word, int offset) {
  LList* res;
  int i;
  TNode* currnode = tj->mynode;
  res = (LList*) malloc (sizeof(LList));
  res->myloc = &(currnode->kids.locator);
  res->degleft = tj->numdeg;

  for (i = 0; i < tj->numdeg; i++)
    res->degloc[i] = (char *)(word - tj->degloc[i]);
  res->next = 0;
  return res;
}

LList* getNextWords (TNode* currnode, char* word, int ndegen) {
  int i, j;
  int height = currnode->height;
  int letter = lookup(*word);
  int mynjobs;
  char mydirty;
  char myflags;
  char first = 0;
  LList* res=0, *temp;

  // -1 --> 0 (second param)
  if (letter >= 0 && numjobs == 0) /*new string*/
    addjob(currnode, 0, 0, -1);
  mydirty = jobqueue[0].dirty;
  mynjobs = numjobs; /* need my own copy so that I don't go over inserted things */
  for (i = 0; i < mynjobs; i++) {
    myflags = - 1 - (1 << triealphasize)+1;
    first = 0;
    //    printf("jqdl = %d, w = %d, mnh = %d\n", jobqueue[i].degloc[0],(int)word, jobqueue[i].mynode->height);
    if (jobqueue[i].numdeg > 0 && ((char *) jobqueue[i].degloc[0] < word - (height -jobqueue[i].mynode->height))) {
      remjob(i);
      if (jobqueue[i].dirty == mydirty) {
	mynjobs--;
	i--;     
      }
      continue;
    }
    do {
      if (!jobqueue[i].mynode) {
	remjob(i);
	if (jobqueue[i].dirty == mydirty) {
	  mynjobs--;
	  i--;     /* need this if the guy I moved in the old place is in my pass */
	}
	break;
      }
      if (jobqueue[i].mynode->height == 0 || first) {
	jobqueue[i].mynode = jobqueue[i].mynode->backptr;
      }
      first = 1;
      if (ndegen - jobqueue[i].numdeg > 0) {
	for (j = 0; j < triealphasize; j++) {
	  if (!(myflags & (1<< j)) && jobqueue[i].mynode->kids.ptrs[j]) {
	    // changed -1 --> 0
	    addjob(jobqueue[i].mynode->kids.ptrs[j], (j==letter)?0:word, !mydirty,i);
	    if (jobqueue[i].mynode->height == 1) {
	      temp = makeLList(&jobqueue[numjobs-1], word, j);
	      temp->next = res;
	      res = temp;
	    }
	    myflags = myflags | (1 << j); 
	  }
	}
      }
      
      else {
	if (letter >= 0 && jobqueue[i].mynode->kids.ptrs[letter]) {
	  jobqueue[i].mynode = jobqueue[i].mynode->kids.ptrs[letter];
	  jobqueue[i].dirty = !mydirty;
	  if (jobqueue[i].mynode->height == 0) {
	    temp = makeLList(&jobqueue[i], word, letter);
	    temp->next = res;
	    res = temp;
	  }
	  myflags = -1;
	}
      }
      if (myflags == -1) {
	break;
      }
    } while(jobqueue[i].mynode != jobqueue[i].mynode->backptr);
    if (jobqueue[i].dirty == mydirty) {
      remjob(i);
      if (jobqueue[i].dirty == mydirty) {
	mynjobs--;
	i--;     /* need this if the guy I moved in the old place is in my pass */
      }
    }
  }
  return res;
}






