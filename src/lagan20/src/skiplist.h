#define MAX_LISTS 32

typedef struct skiplistelem {
  struct skiplistelem** next;
  struct skiplistelem** prev;
  int linkcnt;
  int index;
  void* myelem;
} sle;

typedef struct skiplist {
  sle* sentinel;
  int maxlevel;
} sklst;


void initLib();
sklst* makeSkLst();
void chklst(sklst* trgt);
void delSkLst(sklst* trgt);
sle* SLinsertAfter(sklst* trgt, sle* prev, int index, void* elem);
sle* SLinsert(sklst* trgt, int index, void* elem);
sle* SLgetLast(sklst* trgt);
void SLremove(sklst* trgt, sle* tbr);
sle* SLfind(sklst* trgt, int index);
sle* SLlowFind(sklst* trgt, int index);
sle* mksle(int linkcnt, int index, void* myelem);
void delSLE(sle* tbd);

