#include "fchaos.h"
#define MAX_DEGEN 2


int indeces[256];

typedef struct PrevHits {
  int* inds1;
  int* inds2;
  int numind;
} phits;

typedef struct Locator {
  int* locs;
  int numlocs;
  int locssize;
} locs;

typedef struct LocatorList {
  locs* myloc;
  int degleft;
  char* degloc[MAX_DEGEN];
  struct LocatorList* next;

  /* Stuff below is for chaining */
  int location;
  char* toberemoved;
  float* scores;
  int* seq1startpnt;
  int* seq2startpnt;
  int* seq1endpnt;
  int* seq2endpnt;
  phits* myhits;
  sle** mysles;
} LList;

typedef struct TrieNode {
  union children {
    struct TrieNode** ptrs;
    locs locator;
  } kids;
  struct TrieNode* backptr;   /* added for threading */
  int height;
} TNode;

typedef struct TrieJob {
  TNode* mynode;
  int numdeg;
  char *degloc[MAX_DEGEN];
  char dirty;
} TJob;

LList* appendLList(LList* a , LList* b);
LList* savenfreeLList (LList* tbf, seq* seq1, seq* seq2);
TNode* makeTrie(int height, char* alphabet); 
void freeTrie (TNode* root);
TNode* makeNode(int height);
int insertWord(TNode* root, char* word, char* strbeg);
LList* lookupWord(TNode* currnode, char* word, int ndegen);

/* above this line are things for all tries */

/*this is for threaded stuff */
void cleanJobQueue();
LList* getNextWords(TNode* root, char* word, int ndegen);
void insertString(TNode* root, char* tbi);

