#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "skiplist.h"

typedef struct GapFreeChunkList {
  int x;
  int y;
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
  struct HitLocationList *next;
  struct HitLocationList *bkptr;
  gfc* first;
  gfc* last;
  float scoreSoFar;
} hll;

typedef struct hllpointer {
  int number;
  char isstart;
  hll* myhll;
} hptr;

char seq1name[255];
char seq2name[255];

float gapopen =0, gapcont=0;
int gapfreechunks = 0;
hll* parseCHAOS(FILE* infile, int* numhits);
hll* findBestChain(hptr* myarr, int arrsize);
void doOutput(hll* mylist);
hll* sortList(hll* mylist);


static int hptrcomp (const void *p1, const void *p2) {
  int i = ((hptr*)p1)->number;
  int j = ((hptr*)p2)->number;
  int it = ((hptr*)p1)->isstart;
  int jt = ((hptr*)p2)->isstart;
  if (i > j)
    return (1);
  if (i < j)
    return (-1);
  if (it)
    return -1;
  else 
    return 1;
}

int main(int argc, char** argv){
  FILE* inf;
  hll* mylist, *temp, *best;
  int numhits, i=0;
  hptr* myptrs;
  
  if (argc < 1 || argc > 6) {
    printf("usage: anchors [filename] [-gap # #]\n");
    printf("For -gap the first # is the gap open penalty, the second the gap continue");
    return 1;
  }
  i = 2;
  if (argc == 1 || strchr(argv[1], '-')) {
    i = 1;
    inf = stdin;
  }
  else if (!(inf = fopen(argv[1],"r"))) {
    printf("couldn't open input file\n");
    return 2;
  }
  while  (i < argc) {
    if (!strcmp(argv[i], "-gap")) {
      sscanf(argv[i+1],"%f",&gapopen);
      sscanf(argv[i+2],"%f",&gapcont);
      i += 3;
    }
    else if (!strcmp(argv[i], "-gfc")) {
      gapfreechunks = 1;
      i += 1;
    }
  }
  initLib();

  mylist = parseCHAOS(inf, &numhits);
  if (!numhits)
    return 0;
  myptrs = (hptr*) malloc (sizeof(hptr) * numhits *2);
  i = 0;
  for (temp = mylist; temp; temp = temp->next) {
    myptrs[i].number  = temp->seq1start;
    myptrs[i].isstart = 1;
    myptrs[i].myhll = temp;
    myptrs[i+1].number  = temp->seq1end;
    myptrs[i+1].isstart = 0;
    myptrs[i+1].myhll = temp;
    i = i+2;
  }
  qsort(myptrs, numhits*2, sizeof(hptr), hptrcomp);
  best = findBestChain(myptrs, numhits*2);
  doOutput(best);
  return 0;
}

int whRulez(hll* one, hll* two) {
  float gapdiff = ((float)(two->seq2end - one->seq2end)) * gapcont;
  return two->scoreSoFar-one->scoreSoFar-gapdiff > 0;
}

float gapPen(hll* next, hll* prev) {
  float j= ((float)(next->seq2start-prev->seq2end))*gapcont + gapopen;
  //  printf("%d (%f)*(%f) %f gap\n", next->seq2start-prev->seq2end, ((float)(next->seq2start-prev->seq2end)),gapcont,j);
  return j;
}

hll* findBestChain(hptr* array, int arrsize) {
  sklst* skipper = makeSkLst();
  sle* help, *bestptr;
  float best = -1;
  int i;
  for (i = 0; i < arrsize; i++) {
    if (array[i].isstart) {
      help = SLfind(skipper, array[i].myhll->seq2start);
      if (help->myelem && 
	  (gapPen(array[i].myhll, ((hll*)help->myelem)) + ((hll*)help->myelem)->scoreSoFar) > 0) {
	array[i].myhll->bkptr = help->myelem;
	array[i].myhll->scoreSoFar = ((hll*)help->myelem)->scoreSoFar + array[i].myhll->score + gapPen(array[i].myhll, ((hll*)help->myelem));
      }
      else {
	array[i].myhll->bkptr = 0;
	array[i].myhll->scoreSoFar = array[i].myhll->score;
      }
    }
    else {
      help = SLfind(skipper, array[i].myhll->seq2end);

      if (help->myelem && whRulez(array[i].myhll,((hll*)help->myelem)))
	continue;
      SLinsertAfter(skipper, help, array[i].myhll->seq2end, array[i].myhll);
      help = help->next[0];

      while (help->next[0] && 
	     !whRulez(((hll*)help->myelem), ((hll*)help->next[0]->myelem)))
	SLremove(skipper, help->next[0]);
    }
  }
  help = skipper->sentinel->next[0];
  while (help) {
    if (((hll*)help->myelem)->scoreSoFar > best) {
      best = ((hll*)help->myelem)->scoreSoFar;
      bestptr = help;
    } 
    help = help->next[0];
  }

  return (hll*)bestptr->myelem;
}

void doOutput(hll* best) { 
  int len;

  hll *bestPtr=best, *temp;
  int chl=0, i, bestscore=-1;
  gfc* tmpgf;
  for (temp = bestPtr; temp; temp = temp->bkptr) {
    chl++;
  }

  for (temp = bestPtr; temp; temp = temp->bkptr) {
    len = temp->seq1end - temp->seq1start + 1 ;
    if (!gapfreechunks || !temp->first) {
      printf("(%d %d)=",temp->seq2start, temp->seq2end);
      printf("(%d %d) %f\n",temp->seq1start, temp->seq1end, temp->score);
    }
    else {
      for (tmpgf = temp->first; tmpgf ; tmpgf = tmpgf->next) {
	printf("(%d %d)=(%d %d) %d\n", tmpgf->y, tmpgf->y + tmpgf->length-1, tmpgf->x, tmpgf->x + tmpgf->length-1, 
	       tmpgf->score);
	
      }
    }
  }
}

char* rolltonum(char* str) {
  char *got1=0, *got2=0;
  int in=0, i=0;
  while (1) {
    if (str[i] == 0) {
      break;
    }
    if (str[i] == ';' && got1 && got2){
      return got1;
    }
    if (isdigit(str[i])) {
      if (!in && (!i || isspace(str[i-1]))) { 
	if (got1) 
	  got2 = &str[i];
	else 
	  got1 = &str[i];
	in = 1;
      }
    }
    else if (in && (isspace(str[i]))) {
      if (got2) {
	got1 = got2; got2=0; in = 0;
      }
      in = 0;
    }

    else {
      in = 0;
      got1=got2=0;
    }
    i++;
  }
  return &str[i];
}

int getlineLagan(FILE* infile, hll* tt) {
  char temp[1024];
  char* help;
  int z, h;
  fgets(temp, 1024, infile);
   help = rolltonum(temp);
  z = sscanf(help, "%d %d;%n", &tt->seq2start, &tt->seq2end, &h);
  if (z<2)
    return 0;
  help = rolltonum(help+h);
  if (sscanf(help,"%d %d; score = %f (%*c)\n", &tt->seq1start,
	     &tt->seq1end,&tt->score)<3)
    return 0;
  return 1;
}


hll* parseCHAOS(FILE* infile, int* totnum) {
  hll *myres=0, *tt;
  gfc* temp;
  *totnum = 0;
  while(!feof(infile)) {
    tt = (hll*) malloc(sizeof(hll));
    while (!feof(infile) && !getlineLagan(infile, tt))
      ;
    if (feof(infile)) break;
    if (gapfreechunks) {
      tt->first = tt->last = temp = (gfc*) malloc(sizeof (gfc));
      temp->next = 0;
      while (fscanf(infile, "%d %d %d %d", &temp->y, &temp->x, &temp->length, &temp->score) == 4){
	tt->first = temp; 
	temp = (gfc*) malloc(sizeof (gfc));
	temp->next = tt->first;
      }
      free(temp);
      if (temp == tt->last) {
	tt->first = tt->last = 0;
      }
    }
    tt->next = myres;
    tt->bkptr = 0;
    tt->scoreSoFar = 0;
    (*totnum)++;
    myres = tt;
  }
  return myres;
}





