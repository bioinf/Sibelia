#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "fchaos.h"
#include "skiplist.h"
#include "thrtrie.h"
#include "global.h"
#include "translate.h"
#include "filebuffer.h"

#define VER_NUM "0.932"
#define BLOSUM_FILE "blosum62s.txt"
#define BLOSUM_FILE_SIZE 24
#define NUC_FILE "nucmatrix.txt"
#define NUC_FILE_SIZE 6

#define MAX2(x,y)   ( (x) >= (y) ? (x) : (y) )
#define MIN2(x,y)   ( (x) <= (y) ? (x) : (y) )
#define ABS(x)   ( ((x) >= (0)) ? (x) : (-x) )
#define WEQ2(x,y,a)  (((x)==(a))? 0: ((y)==(a))? 1:-1)
#define MIN(A,B) (A>B)?B:A
#define MAX(A,B) (A>B)?A:B

typedef struct SeqMatch {
  LList* myll;
  int offset;
} match;

extern int indeces[256];


void remElem(LList* tbf, int i);

int verbose = 0;
int wordlen = 10;
int ndegen = 1;
int cutoff = 25;
int lookback = 20;
int gapfreechunks = 0;
int mgaplen = 5;
int gappenc = -1;
int gappeno = 0 ;
int both = 0;
int translated = 0;
int s1start = 0;
int s1end = 0;
int s2start = 0;
int s2end = 0;

int extend = 0;
int reScoreCutoff = 0;

//int matchsco = 12;
//int mismatchsco = -8;

int gappenstart = -1500;
int gappenext = -50;
int dropcutoff = 1500;

int substmatrix[256][256];


hll* allhits = 0;
sklst* mylist;
int gapstart=20;
int gapcont=1;
char* alpha = "ATCGN";
char* triealpha = "ATCG";
char* protalpha = "PCMH[DE][KR][NQ][ST][ILV][FYW][AG]X*";
char* prottriealpha = "PCMH[DE][KR][NQ][ST][ILV][FYW][AG]";
char direction;

FILE* pairfile = 0;


char comp(char c) {
  switch(c) {
  case 'a': case 'A': return 'T'; 
  case 't': case 'T': return 'A'; 
  case 'c': case 'C': return 'G';
  case 'g': case 'G': return 'C'; 
  case 'n': case 'N': return 'N';
  default: printf("ERROR, Bad letter to RC: %c\n",c); return -1;
  }
}

void revComplement(char* a) {
  int length = strlen(a);
  char lft;
  int i;
  for (i=0; i < length/2; i++) {
    lft = a[i];
    a[i] = comp(a[length-i-1]);
    a[length-i-1] = comp(lft);
  }
  if (length % 2)
    a[length/2] = comp(a[length/2]);
}

void freeSeq (seq* tbf) {
  free(tbf->name);
  free(tbf->rptr);
  free(tbf);
}

void freeHLL (hll* tbf) {
  gfc *t = tbf->first;
  gfc *n;
  while (t) {
    n = t->next;
    free (t);
    t = n;
  }
  free (tbf);
}

void printHLL(hll* res,  seq* query, seq* dbase, int len) {
  hll* temp;
  align* myal;
  gfc* tmpgf;
  int currx, curry;
  char *qptr = query->lets, *dptr = dbase->lets;
  if (direction == '+') {
    while (res) {
      if (s1start > 0) {
	res->seq1start += (s1start-1);
	res->seq1end += (s1start-1);
	query->lets = query->rptr;
      }
      if (s2start > 0) {
	res->seq2start += (s2start-1);
	res->seq2end += (s2start-1);
	dbase->lets = dbase->rptr;
      }
      printf("%s %d %d; %s %d %d; score = %f (%c)\n", query->name, 
	     res->seq1start+1, res->seq1end+1, 
	     dbase->name, res->seq2start+1, res->seq2end+1, 
	     res->score,direction);
      if (verbose) {
	myal = global(query->lets, res->seq1start, 
		      res->seq1end, dbase->lets, res->seq2start, 
		      res->seq2end, gapstart, gapcont);
	printalign(query->lets, res->seq1start, 
		   res->seq1end, dbase->lets, res->seq2start, 
		   res->seq2end, myal);
      }
      if (gapfreechunks) {
	currx = res->seq1start+1;
	curry = res->seq2start+1;
	tmpgf = res->first;
	while (tmpgf) {
	  if (tmpgf->length) {
	    printf ("%d %d %d %d\n", currx, curry, tmpgf->length, tmpgf->score);
	    currx += tmpgf->length;
	    curry += tmpgf->length;
	  }
	  tmpgf = tmpgf->next;
	  if (!tmpgf)
	    break;
	  if (tmpgf->offset > 0) {
	    curry += tmpgf->offset;
	  }
	  else {
	    currx -= tmpgf->offset;
	  }
	}
      }
      temp = res;
      res = res->next;
      freeHLL(temp);
    }
  }
  else {
    while (res) {
      if (s1start > 0) {
	res->seq1start += (s1start-1);
	res->seq1end += (s1start-1);
	query->lets = query->rptr;
      }
      if (s2start > 0) {
	res->seq2start += (len-s2end);
	res->seq2end += (len-s2end);
      }

      printf("%s %d %d; %s %d %d; score = %f (%c)\n", query->name, 
	     res->seq1start+1, res->seq1end+1, 
	     dbase->name, len-(res->seq2start), len - (res->seq2end), 
	     res->score, direction);
      if (verbose) {
	myal = global(query->lets, res->seq1start, 
		      res->seq1end, dbase->lets, 
		      res->seq2start, res->seq2end, gapstart, gapcont);
	printalign(query->lets, res->seq1start, 
		   res->seq1end, dbase->lets, 
		   res->seq2start, res->seq2end, myal);
      }
      if (gapfreechunks) {
	currx = res->seq1start+1;
	curry = len - res->seq2start;
	tmpgf = res->first;
	while (tmpgf) {
	  if (tmpgf->length) {
	    printf ("%d %d %d %d \n", currx, curry, tmpgf->length, tmpgf->score);
	    currx += tmpgf->length;
	    curry -= tmpgf->length;
	  }
	  tmpgf = tmpgf->next;
	  if (!tmpgf)
	    break;
	  if (tmpgf->offset < 0) {
	    currx -= tmpgf->offset;
	  }
	  else {
	    curry -= tmpgf->offset;
	  }
	}
      }
      temp = res;
      res = res->next;
      freeHLL(temp);
    }
  }
  query->lets=qptr;
  dbase->lets = dptr;
}


void printList (hll *ptr){
  if (ptr){
    fprintf (stderr, "(%d %d)=(%d %d) %f\n", ptr->seq1start, ptr->seq1end, ptr->seq2start, ptr->seq2end, ptr->score);
    printList (ptr->next);
  }
}

int compare (hll *list1, hll *list2){
  return (list1->seq1start < list2->seq1start) ||
    (list1->seq1start == list2->seq1start && list1->seq1end > list2->seq1end);
}

hll* merge2(hll* list1, hll* list2) {
  hll *totallist = 0, *temp = 0;

  if (!list1) return list2;
  if (!list2) return list1;

  while (list1 || list2) {
    if (list1 && (!list2 || compare (list1, list2))){
      if (!totallist)
	totallist = temp = list1;
      else {
	temp->next = list1;
	temp = temp->next;
      }
      list1 = list1->next;
    }
    else {
      if (!totallist)
	totallist = temp = list2;
      else {
	temp->next = list2;
	temp = temp->next;
      }
      list2 = list2->next;
    }
  }
  temp->next = 0;
  return totallist;
}

hll* findmiddle(hll* mylist) {
  hll* other = mylist->next;
  while (other && other->next) {
    other = other->next->next;
    mylist = mylist->next;
  }
  return mylist;
}

hll* sortList(hll* mylist) {
  hll* premid; 
  hll* mid;

  if (!mylist || !mylist->next)
    return mylist;

  premid = findmiddle(mylist);
  mid = premid->next;
  premid->next = 0;
  mylist = sortList(mylist);
  mid = sortList(mid);
  return merge2(mylist,mid);
}

int duplicates(hll* f, hll* s) {
  return (s->seq2start >= f->seq2start) && (s->seq2end <= f->seq2end);
}

hll* removeDups(hll* allhits, seq* seq1, seq* seq2) {
  hll *i, *j, *jprev, *temp;
  for (i = allhits; i; i = i->next){
    jprev = i;
    for (j = i->next; j && (j->seq2start >= i->seq2end) ; j = j->next){
      if (duplicates (i, j) || mergeOverlap (i, j, seq1, seq2)){
	jprev->next = j->next;
	freeHLL (j);
	j = jprev;
      }
      else {
	jprev = j;
      }
    }
  }

  allhits = sortList (allhits);
  for (i = allhits; i; i = i->next){
    jprev = i;
    for (j = i->next; j && (j->seq1start <= i->seq1end) ; j = j->next){
      if (duplicates (i, j) || mergeOverlap (i, j, seq1, seq2)){
	jprev->next = j->next;
	freeHLL (j);
	j = jprev;
      }
      else {
	jprev = j;
      }
    }
  }

  return allhits;
}


seq* readfile(FILE* input, int seqnum) {
  char* res = (char*) malloc(sizeof(char));
  int ressize = 1, numread=0;
  char temp[256];
  seq* myseq = (seq*) malloc(sizeof(seq));
  char currchar;
  if (feof(input))
    return 0;
  fgets(temp, 255, input);
  if (temp[0] != '>') {
    fprintf(stderr, "File is not in FASTA format!!\n");
    exit(1);
  }
  myseq->name = (char*) malloc((strlen(temp))*sizeof(char));
  strcpy(myseq->name, temp+1);
  *(strchr(myseq->name, '\n')) = 0;
  currchar = fgetc(input);
  while ((currchar != '>') && (currchar != EOF)) {
    if (!isspace(currchar)) {
      currchar = toupper(currchar);
      if (!strchr(alpha, currchar)) {
	fprintf(stderr, "WARNING %c converted to N\n", currchar, alpha);
	currchar = 'N';
      }
      res[numread++] = currchar;
      if (numread >= ressize) {
	res=(char*)realloc(res, sizeof(char)*(ressize*=2)); 
      }
    }
    currchar = fgetc(input);
  }
  if (currchar == '>')
    ungetc(currchar, input);
  res[numread]=0;

  myseq->rptr = res;
  if (seqnum == 1) {
    if (s1start > 0) {
      res[s1end] = 0;
      res = &res[s1start-1];
      numread = s1end-s1start+1;
    }
  }
  else {
    if (s2start > 0) {
      res[s2end] = 0;
      res = &res[s2start-1];
      numread = s2end-s2start+1;

    }
  }
  myseq->lets = res;
  myseq->numlets = numread;
  return myseq;
}

int isin (char* arr, int size, int elem) {
  while (--size>=0) {
    if (arr[size] == elem)
      return 1;
  }
  return 0;
}

int chain(LList* second, int off2, LList* first, int off1, int diff1, int gap, float baseval) {
  int i, d1=0, d2=0;
  int diff2 = second->myloc->locs[off2] - first->myloc->locs[off1];
  int mindiff;
  int score=wordlen-second->degleft;

  gap = abs(gap)*gappenc + gappeno;

  if (diff2 <= 0  || diff2 >= lookback)
    return -1;

  if (diff1 >= wordlen && diff2 >= wordlen) {
    return score*baseval+gap;
  }
  mindiff = MIN(diff1, diff2);
  /* TODO
    for (i=second->degleft-1; i >=0; i--) {
    printf(" %d   %d %d \n", second->degloc[i], diff1, diff2);
    if (!d1 && second->degloc[i] - diff1 <= 0)
    d1 = 1;
    if (&d2 && second->degloc[i] - diff2 <= 0)
    d2 = 1;
    if (d1 || d2) {
    break;
    }
    }   
  */
  return mindiff*baseval+gap;
}

int tc =0;
int wc = 0;

inline void findPrev(LList* curr, int position, int offset, float baseval) {
  int j,k;
  LList* temp;
  sle* iterator;
  float bestscore = 0;
  LList* bestelem = 0;
  int bestoffset = -1;
  int doneset = 0;
  int tempscore, myscore = wordlen - curr->degleft;

  tc++;
  iterator = SLfind(mylist, position-curr->myloc->locs[offset]-mgaplen+1);
  if (iterator)  {
    curr->mysles[offset] = iterator;      
  }
  if (iterator && 
      iterator->index <= position-curr->myloc->locs[offset]-mgaplen) {
    iterator = iterator->next[0];
  }

  if (iterator && (iterator->index  < position-curr->myloc->locs[offset]))  {
    curr->mysles[offset] = iterator;      
  }

  while (iterator && 
	 (iterator->index < position-curr->myloc->locs[offset]+mgaplen)) {
    if (iterator->next[0] && (iterator->index  < position-curr->myloc->locs[offset]) && 
	(iterator->next[0]->index  >= position-curr->myloc->locs[offset]))  {
      curr->mysles[offset] = iterator;      
    }
    temp = ((match*)iterator->myelem)->myll;
    k = ((match*)iterator->myelem)->offset;
    j = position-temp->location;
    tempscore = chain(curr, offset, temp, k,j, iterator->index - position+curr->myloc->locs[offset], baseval);
    if (tempscore > 0) {
      if (temp->scores[k]+tempscore > bestscore) {
	bestscore = temp->scores[k]+tempscore;
	bestelem = temp;
	bestoffset=k;
      }
      else {
	temp->scores[k] = -1;
      }
    }
    /*    printf("it = %x next = %x\n", iterator, iterator->next[0]); */
    iterator = iterator->next[0];
    if (temp->toberemoved[k]) {
      remElem(temp, k);
      temp->mysles[k] = 0;
    }
  }
  if (bestelem) {
    wc++;
    curr->scores[offset] = bestscore;
    /*    printf("offs = %d, numlocs = %d\n",offset, curr->myloc->numlocs);*/
    curr->seq1startpnt[offset] = bestelem->seq1startpnt[bestoffset];
    curr->seq2startpnt[offset] = bestelem->seq2startpnt[bestoffset];
    curr->myhits[offset].inds1 = (int*) malloc (sizeof(int)*(bestelem->myhits[bestoffset].numind+1));
    curr->myhits[offset].inds2 = (int*) malloc (sizeof(int)*(bestelem->myhits[bestoffset].numind+1));
    curr->myhits[offset].numind = bestelem->myhits[bestoffset].numind+1;

    memcpy (curr->myhits[offset].inds2, bestelem->myhits[bestoffset].inds2,
	    bestelem->myhits[bestoffset].numind*sizeof(int));
    memcpy (curr->myhits[offset].inds1, bestelem->myhits[bestoffset].inds1,
	    bestelem->myhits[bestoffset].numind*sizeof(int));
    curr->myhits[offset].inds2[bestelem->myhits[bestoffset].numind] = position;
    curr->myhits[offset].inds1[bestelem->myhits[bestoffset].numind] = 
      (int) curr->myloc->locs[offset];

  }
  else { 
    curr->scores[offset] = myscore; 
    curr->seq2startpnt[offset] = position;
    curr->seq1startpnt[offset] = (int)curr->myloc->locs[offset];
    curr->myhits[offset].inds1 = (int*) malloc (sizeof(int));
    curr->myhits[offset].inds2 = (int*) malloc (sizeof(int));
    curr->myhits[offset].inds2[0] = position;
    curr->myhits[offset].inds1[0] = (int)curr->myloc->locs[offset];
    curr->myhits[offset].numind = 1;
  }
}

void connectToPrev(LList* curr, int index, float baseval) {
  int j;
  curr->scores = (float*) malloc(sizeof(float) * curr->myloc->numlocs);
  curr->myhits = (phits*) malloc(sizeof(phits) * curr->myloc->numlocs);
  curr->toberemoved = (char*) malloc(sizeof(char) * curr->myloc->numlocs);
  curr->seq1startpnt = (int*) malloc(sizeof(int) * curr->myloc->numlocs);
  curr->seq2startpnt = (int*) malloc(sizeof(int) * curr->myloc->numlocs);
  curr->seq1endpnt = (int*) malloc(sizeof(int) * curr->myloc->numlocs);
  curr->seq2endpnt = (int*) malloc(sizeof(int) * curr->myloc->numlocs);
  curr->mysles = (sle**) malloc(sizeof(sle*) * curr->myloc->numlocs);
  for (j = 0; j < curr->myloc->numlocs; j++) {
    curr->toberemoved[j] = 0;
    curr->myhits[j].numind = 0;
    curr->scores[j] = 0;
    curr->seq1startpnt[j] = 0;
    curr->seq2startpnt[j] = 0;
    curr->mysles[j] = 0;
    findPrev(curr,index,j,baseval);
  }
}

int doAlgo(TNode* root, seq* query, seq* dbase) {
  char* currword = dbase->lets;
  LList** LListArr = (LList**) malloc(sizeof(LList*) * dbase->numlets);
  LList* temp;
  match* mattemp;
  int i = 0, j;
  float bestscore=-1, baseval;
  int bestqueryloc=-1, bestdbaseloc=-1, numhits;
  while (*currword) {

    if (!(i%10000)) {
      //      fprintf(stderr,"WORKING %d\n",i); 
    }
    if (*currword == '.') {
      /*TODO */
    }
    LListArr[i] = temp = getNextWords(root, currword++, ndegen);

    /*****/
    numhits = 1;
    while (temp){
      numhits += temp->myloc->numlocs;
      temp = temp->next;
    }
    baseval = (float) log ((double) query->numsiglets / (double) numhits) / (float) wordlen;
    temp = LListArr[i];
    /*****/
    
    while (temp) {
      temp->location = i-wordlen+1;
      connectToPrev(temp, temp->location, baseval);
      for (j = 0; j < temp->myloc->numlocs; j++) {
	mattemp = (match*) malloc (sizeof(match));
	mattemp->myll = temp;
	mattemp->offset = j;
	if (temp->mysles[j])
	  temp->mysles[j] = SLinsertAfter(mylist, temp->mysles[j], temp->location-(int)temp->myloc->locs[j], mattemp);
	else
	  temp->mysles[j] = SLinsert(mylist, temp->location-(int)temp->myloc->locs[j], mattemp);
      }
      temp = temp->next;
    }
    if (i-lookback >= 0) {
      LListArr[i-lookback] = savenfreeLList(LListArr[i-lookback], query, dbase);
    }
    i++;
  }
  j = (i-lookback>=0)?i-lookback:0;
  for ( ; j < i; j++) {
    LListArr[j] = savenfreeLList(LListArr[j], query,dbase);
  }
  cleanJobQueue();
  free(LListArr);
  //  fprintf(stderr, "%d chained of %d\n", wc , tc);
  return 0;
}

char getLetter (FILE *file){
  char ch;

  while (!feof (file)){
    ch = fgetc (file);
    if (!isspace (ch)) return ch;
  }
  return 0;
}

void readSubstMatrix (char *filename, int size){
  FILE *file;
  char line[1024], *symbs;
  int i, j;

  sprintf (line, "%s/%s", getenv ("LAGAN_DIR"), filename);
  file = fopen (line, "r"); assert (file);
  
  for (i = 0; i < 256; i++){
    for (j = 0; j < 256; j++){
      substmatrix[i][j] = 0;
    }
  }
  
  symbs = (char *) malloc (sizeof (char) * size); assert (symbs);
  for (i = 0; i < size; i++) symbs[i] = getLetter (file);
  for (i = 0; i < size; i++){
    getLetter (file);
    for (j = 0; j < size; j++){
      fscanf (file, "%d", &(substmatrix[(unsigned char) symbs[i]][(unsigned char) symbs[j]]));
    }
  }

  
  fscanf (file, "%d", &gappenstart);
  fscanf (file, "%d", &gappenext);
  gappenstart = (gappenext *= 2);

  fclose (file);
}

void paramParse(int argc, char** argv) {
  int i = 3;

  for ( ; i < argc; i++) {
    if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "-P")) {
      alpha = protalpha;
      triealpha = prottriealpha;
      wordlen = 4;
      lookback = 8;
      //      dropcutoff = 50;
      readSubstMatrix (BLOSUM_FILE, BLOSUM_FILE_SIZE);
    }
    else if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "-V")) {
      verbose = 1;
    }
    else if (!strcmp(argv[i], "-b") || !strcmp(argv[i], "-B")) {
      both = 1;
    }
    else if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "-T")) {
      translated = 1;
      triealpha = prottriealpha;
      wordlen = 4;
      mgaplen = 3;
      lookback = 8;
      //      dropcutoff = 50;
      readSubstMatrix (BLOSUM_FILE, BLOSUM_FILE_SIZE);
    }
    else if (!strcmp(argv[i], "-rsc") || !strcmp(argv[i], "-RSC")) {
      reScoreCutoff  = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-gfc") || !strcmp(argv[i], "-GFC")) {
      gapfreechunks  = 1;
    }
    else if (!strcmp(argv[i], "-ext") || !strcmp(argv[i], "-EXT")) {
      extend = 1;
    }
    else if (!strcmp(argv[i], "-wl") || !strcmp(argv[i], "-WL")) {
      wordlen = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-nd") || !strcmp(argv[i], "-ND")) {
      ndegen = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-co") || !strcmp(argv[i], "-CO")) {
      cutoff = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-lb") || !strcmp(argv[i], "-LB")) {
      lookback = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-gl") || !strcmp(argv[i], "-GL")) {
      mgaplen = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-gs") || !strcmp(argv[i], "-GS")) {
      gappeno = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-gc") || !strcmp(argv[i], "-GC")) {
      gappenc = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-s1") || !strcmp(argv[i], "-S1")) {
      s1start = atoi(argv[++i]);
      s1end = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-s2") || !strcmp(argv[i], "-S2")) {
      s2start = atoi(argv[++i]);
      s2end = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-pairs") || !strcmp(argv[i], "-PAIRS")) {
      if (!(pairfile = fopen(argv[++i],"r"))) {
	printf("couldnt open pairs file %s\n",argv[i]);
	exit (2);
      }
    }
  }

  if (!translated) readSubstMatrix (NUC_FILE, NUC_FILE_SIZE);

}

void usage() {
  printf("usage: \nchaos queryfile dbasefile [options]\n\n");
  printf("Options:\n");
  printf("-p     = Peptide sequence [default genomic]\n");
  printf("-v     = Verbose mode [default brief]\n");
  printf("-b     = Both strands [default forward-only]\n");
  printf("-t     = Translated [default off]\n");
  printf("-ext   = do BLAST-like extention with given cutoff [default off]\n");
  printf("-wl #  = Word Length [default 10 for genomic, 4 for peptide]\n");
  printf("-nd #  = Number of Degeneracy [default 1 for genomic, 0 for peptide]\n");
  printf("-co #  = score CutOff [default 25]\n");
  printf("-rsc # = Rescoring cutoff [default 0]\n");
  printf("-lb #  = LookBack distance [default 20 for genomic, 8 for peptide]\n");
  printf("-gl #  = maximum Gap Length [default 5 for genomic, 3 for peptide]\n");
  printf("-gs #  = Gap Start penalty [default 0]\n");
  printf("-gc #  = Gap Continue penalty [default -1]\n");
  printf("-s1 # # = use the given substring of the query [default whole]\n");
  printf("-s2 # # = use the givensubstring of the dbase [default whole]\n");
  printf("-pairs pairfile = read \"-s1 # # -s2 # #\" from pairfile [default off]\n\t[This is not fully functional!!!]\n");
  printf("-version = prints the version of this CHAOS\n");
}

void rc(seq* dbase) {
  revComplement(dbase->lets);
}


int paircnt = 0;

char savs[2];
int savlocs[2] = {-1,-1};

void procPairs(seq* currquery, seq* currdbase) {
  //  int s1start, s1end, s2start, s2end;
  if (savlocs[0]>=0)
    currquery->rptr[savlocs[0]] = savs[0];
  if (savlocs[1]>=0)
    currdbase->rptr[savlocs[1]] = savs[1];
    
  do {
    //fprintf(stderr,"here\n");
    if (fscanf(pairfile, "-s1 %d %d -s2 %d %d\n", &s1start, &s1end, &s2start, &s2end) < 4) {
      pairfile = 0;
      return;
    }
    currquery->numlets = s1end-s1start+1;
    currdbase->numlets = s2end-s2start+1;
//         fprintf (stderr, "%d %d; %d\n",currquery->numlets,
//            currdbase->numlets, wordlen+1);
  }  while (currquery->numlets < wordlen+1 && currdbase->numlets < wordlen+1)
       ;

  savlocs[0] = s1end;
  savs[0] = currquery->rptr[s1end];
  currquery->rptr[s1end] = 0;
  currquery->lets = &(currquery->rptr[s1start-1]);
  currquery->numlets = s1end-s1start+1;
  savlocs[1] = s2end;
  savs[1] = currdbase->rptr[s2end];
  currdbase->rptr[s2end] = 0;
  currdbase->lets = &(currdbase->rptr[s2start-1]);
  currdbase->numlets = s2end-s2start+1;
  paircnt++;
  if (paircnt%20 ==19)
    fprintf(stderr, "done with %d\n", paircnt);
}

void transloc(hll* myhits, int frseq1, int frseq2, int seq1len, int seq2len) {
  int temp;
  while (myhits) {
    if (frseq1<=2) {
      myhits->seq1start = myhits->seq1start*3 + frseq1;
      myhits->seq1end = myhits->seq1end*3 + frseq1;
    }
    else {
      temp = (seq1len - myhits->seq1start)*3 + frseq1%3;
      myhits->seq1start = (seq1len - myhits->seq1end)*3 + frseq1%3;
      myhits->seq1end = temp;
    }

    if (frseq2<=2) {
      myhits->seq2start = myhits->seq2start*3 + frseq2;
      myhits->seq2end = myhits->seq2end*3 + frseq2;
    }
    else {
      temp = (seq2len - myhits->seq2start)*3 + frseq2%3;
      myhits->seq2start = (seq2len - myhits->seq2end)*3 + frseq2%3;
      myhits->seq2end = temp;
    }
    myhits = myhits->next;
  }
}

void doTranslated(FileBuffer query, FileBuffer dbase) {
  seq *currquery, *currdbase, *temp;
  seq *queryframes[6], *dbaseframes[6];
  char* currword; 
  TNode *roots[6];
  int i, j;
  currquery = FileRead(query, s1start, s1end, VER_FCHAOS);
  currdbase = FileRead(dbase, s2start, s2end, VER_FCHAOS);

  if (pairfile) {
    procPairs(currquery, currdbase);
    if (!pairfile) {
      FileClose (query);
      FileClose (dbase);
      return;
    }
  }
  do {
    for (i = 0; i < 6; i++) {
      queryframes[i] = transSeq(currquery,i);
      roots[i] = makeTrie(wordlen, triealpha);
      currword = queryframes[i]->lets;
      insertString(roots[i],currword);
    }
    mylist = makeSkLst();
    
    while (currdbase) {
      for (i = 0; i < 6; i++) {
	dbaseframes[i] = transSeq(currdbase,i);
      }
      direction = '+';
      for (i=0; i < 6; i++) 
	for (j=(i/3)*3; j < (i/3+1)*3; j++) {
	  //	  fprintf(stderr, "1DOING FRAME %d AGAINST %d\n",i,j);
	  doAlgo(roots[i], queryframes[i], dbaseframes[j]);
	  /****/
	  allhits = removeDups(allhits, queryframes[i], dbaseframes[j]);
	  transloc(allhits, i, j, queryframes[i]->numlets, dbaseframes[j]->numlets);
	  printHLL(allhits, queryframes[i], dbaseframes[j], currdbase->numlets);
	  allhits = 0;
	}
      if (both) {
	direction = '-';
	for (i=0; i < 6; i++) 
	  for (j=(i>2)?0:3; j < ((i>2)?3:6); j++) {
	    //	    fprintf(stderr, "2DOING FRAME %d AGAINST %d\n",i,j);
	    doAlgo(roots[i], queryframes[i], dbaseframes[j]);
	    /****/
	    allhits = removeDups(allhits, queryframes[i], dbaseframes[j]);
	    transloc(allhits, i, j, queryframes[i]->numlets, dbaseframes[j]->numlets);
	    printHLL(allhits, queryframes[i], dbaseframes[j], currdbase->numlets);
	    allhits = 0;
	  }
      }
      temp = currdbase;
      if (!pairfile)
	freeSeq(currdbase);
      currdbase = FileRead(dbase, s2start, s2end, VER_FCHAOS);
    }
    currdbase = temp;
    if (pairfile) {
      procPairs(currquery, currdbase);
      for (i=0; i < 6; i++) {
	freeSeq(queryframes[i]);
	freeTrie(roots[i]);
      }
    }
  } while (pairfile)
      ;
  
  FileClose (query);
  FileClose (dbase);
}

int main(int argc, char** argv) {
  FileBuffer query;
  FileBuffer dbase;

  seq *currquery, *currdbase, *temp; 
  char* currword; 
  TNode* root;
  int i;

  if (argc < 3) {
    if (argc == 2)
      if (!strcmp(argv[1], "-version") || !strcmp(argv[1], "-Version")) {
	printf("CHAOS version %s\n", VER_NUM);
	exit(0);
      }
    usage();
    return 1;
  }
  if (!(query = FileOpen(argv[1]))) {
    printf("couldnt open query file %s\n",argv[1]);
    usage();
    return 2;
  }
  if (!(dbase = FileOpen(argv[2]))) {
    printf("couldnt open dbase file %s\n",argv[2]);
    usage();
    return 2;
  }
  paramParse(argc, argv);
  initLib();

  if (translated) {
    doTranslated(query, dbase);
    return 0;
  }  

  currquery = FileRead(query, s1start, s1end, VER_FCHAOS);
  currdbase = FileRead(dbase, s2start, s2end, VER_FCHAOS);
  if (pairfile) {
    procPairs(currquery, currdbase);
    if (!pairfile) {
      FileClose (query);
      FileClose (dbase);
      return 0;
    }
  }

  do {
    root = makeTrie(wordlen, triealpha);
    mylist = makeSkLst();
    currword = currquery->lets;
    insertString(root,currword);

    while (currdbase) {
      direction = '+';
      doAlgo(root, currquery, currdbase);
      /***/
      allhits = removeDups(allhits, currquery, currdbase);
      printHLL(allhits, currquery, currdbase, currdbase->numlets);
      allhits = 0;
      if (both) {
	direction = '-';
	rc(currdbase);
	doAlgo(root, currquery, currdbase);
	/****/
	allhits = removeDups(allhits, currquery, currdbase);
	printHLL(allhits, currquery, currdbase, currdbase->numlets);
	allhits = 0;
      }
      temp = currdbase;
      if (!pairfile) {
	freeSeq(currdbase);
      }
      currdbase = FileRead(dbase, s2start, s2end, VER_FCHAOS);
    }
    currdbase = temp;
    if (pairfile) {
      procPairs(currquery, currdbase);
      freeTrie(root);
    }
  } while (pairfile)
      ;

  FileClose (query);
  FileClose (dbase);
  return 0;

}

void saveScore(LList* final, int index, gfc* first, gfc* last) {
  
  hll* myhit = (hll*) malloc(sizeof(hll));
  int temp;

  myhit->score = final->scores[index];
  myhit->seq1end = final->seq1endpnt[index]; 
  myhit->seq2end = final->seq2endpnt[index]; 
  myhit->seq1start = final->seq1startpnt[index];
  myhit->seq2start = final->seq2startpnt[index];
  myhit->last = last;
  myhit->first = first;
  myhit->next = allhits;
  allhits = myhit;
} 

void remElem(LList* tbf, int i) {
  free(tbf->mysles[i]->myelem);
  SLremove(mylist, tbf->mysles[i]);
}

inline int CHmatchscore(unsigned char a, unsigned char b) {
  return substmatrix[a][b];
  /*
  if (translated)
    return substmatrix[a][b];    
  if (a == 'N' || b == 'N' || a == 'X' || b == 'X')
    return 0;
  if ((a == '*' || b == '*') && a != b)
    return -50;
  if (indeces[a] == indeces[b])
      return matchsco;
  return mismatchsco;  
  */
}

int extendBLAST(int s1i, int s2i, char* s1, char* s2, int s1l, int s2l, int dir) {
  int peak=0, peakloc = 0, currscore=0, i = 1;
  while (peak - currscore < dropcutoff) {
    if (s1i+dir*i < 0 || s2i+dir*i < 0 || !s1[s1i+dir*i] || !s2[s2i+dir*i] || s1i+dir*i >= s1l || s2i+dir*i >= s2l)
      break;
    currscore += CHmatchscore (s1[s1i+dir*i], s2[s2i+dir*i]);
    //    fprintf(stderr, "%d(%c %c) ", currscore, s1[s1i+dir*i], s2[s2i+dir*i]);
    if (currscore > peak) {
      peak = currscore;
      peakloc = i;
    }
    i++;
  }
  //  fprintf(stderr, "got to %d, score %d(%d)\n", i, currscore, peak);
  return peakloc;
}

int extendMerge(int s1l, int s2l, int s1r, int s2r, char* s1, char* s2, int* dir) {

  int length, i;
  int *s1arr, *s2arr, bestscore=-9999999, bestloc=0;

  // HACK
  if (s1l < 0){ int err = -s1l; s1l += err; s2l += err; }
  if (s2l < 0){ int err = -s2l; s1l += err; s2l += err; }

  length = MIN2(s1r-s1l, s2r-s2l);

  //  fprintf(stderr,"extmerge (%d %d) (%d %d)\n", s1l, s2l, s1r, s2r);
  *dir = WEQ2(s1r-s1l, s2r-s2l, length);  //0 vertical, 1 horizontal
  if (length <= 0)
    return 0;
  s1arr = (int*) malloc (sizeof(int) * (length+1));
  s2arr = (int*) malloc (sizeof(int) * (length+1));
  s1arr[0] = s2arr[length] = 0;
  for (i = 1; i <= length; i++) {
    s1arr[i] = s1arr[i-1] + CHmatchscore(s1[s1l+i], s2[s2l+i]);
    s2arr[length-i] = s2arr[length-i+1] + CHmatchscore(s1[s1r-i], s2[s2r-i]);
  }
  for (i = 0; i < length; i++) {
    if (s1arr[i]+s2arr[i+1] > bestscore) {
      bestscore = s1arr[i]+s2arr[i+1];
      bestloc = i;
    }
  }
  //  fprintf(stderr, "extMer score = %d\n", bestscore);
  free (s1arr);
  free (s2arr);
  return bestloc;
}

int reScore(int s1l, int s2l, int len, char* s1, char* s2) {
  int i;
  int totscore = 0;

  // HACK
  if (s1l < 0){ int err = -s1l; s1l += err; s2l += err; len -= err; }
  if (s2l < 0){ int err = -s2l; s1l += err; s2l += err; len -= err; }

  for (i=0; i < len; i++) {
    totscore += CHmatchscore(s1[s1l+i], s2[s2l+i]);
  }
  return totscore;
}


void reScoreHit(LList* tbf, int index, char* s1, char* s2, int s1l, int s2l, gfc **frstgf, gfc **mygf) {
  int totscore = 0, myscore;  
  int ts1, ts2, te1, te2;
  int i=0, temp=0, offset, dir;


  if (extend) {
    temp = extendBLAST(tbf->myhits[index].inds1[i], tbf->myhits[index].inds2[i],
		       s1, s2, s1l, s2l, -1);
  }

  tbf->seq1startpnt[index] = ts1 = tbf->myhits[index].inds1[i] - temp;
  tbf->seq2startpnt[index] = ts2 = tbf->myhits[index].inds2[i] - temp;
  *frstgf = *mygf = (gfc*) malloc (sizeof (gfc));
  (*frstgf)->offset = 0;
  
  for (i = 0; i < tbf->myhits[index].numind-1; i++) {
    if (!(offset = ((tbf->myhits[index].inds1[i]-tbf->myhits[index].inds2[i]) -
		    (tbf->myhits[index].inds1[i+1]-tbf->myhits[index].inds2[i+1])))) {

      continue;
    }
    else {
      
      
      temp = extendMerge(tbf->myhits[index].inds1[i]+wordlen-1, 
			 tbf->myhits[index].inds2[i]+wordlen-1,
			 tbf->myhits[index].inds1[i+1], 
			 tbf->myhits[index].inds2[i+1], s1, s2, &dir);
      te1 = tbf->myhits[index].inds1[i] + wordlen - 1 + temp; 
      te2 = tbf->myhits[index].inds2[i] + wordlen - 1 + temp; 

      myscore = reScore(ts1, ts2, te1-ts1+1, s1, s2);
      totscore += myscore;
      totscore += (gappenstart + gappenext * ABS(offset));
      (*mygf)->length = te1-ts1+1;
      (*mygf)->score = myscore;
      (*mygf)->next = (gfc*) malloc (sizeof (gfc));
      (*mygf) = (*mygf)->next;
      (*mygf)->offset = offset;

      if (dir) {
	ts1 = te1+ABS(offset)+1;
	ts2 = te2+1;
      }
      else {
	ts2 = te2+ABS(offset)+1;
	ts1 = te1+1;
      }
    }
  }
  temp = 0;
  if (extend) {
    temp = extendBLAST(tbf->myhits[index].inds1[i]+wordlen-1, 
		       tbf->myhits[index].inds2[i]+wordlen-1, s1, s2, s1l, s2l, 1);
  }
  myscore = reScore(ts1, ts2, tbf->myhits[index].inds1[i]+wordlen-ts1+temp, s1, s2);
  (*mygf)->length = tbf->myhits[index].inds1[i]+wordlen-ts1+temp;
  (*mygf)->score = myscore;
  (*mygf)->next = 0;
  totscore += myscore;
  tbf->scores[index] = totscore;
  tbf->seq1endpnt[index] = tbf->myhits[index].inds1[i]+wordlen-1 + temp;
  tbf->seq2endpnt[index] = tbf->myhits[index].inds2[i]+wordlen-1 + temp;
}


LList* savenfreeLList(LList* tbf, seq* seq1, seq* seq2) {
  int i,j;
  LList* next;
  gfc *first, *last;
  if (!tbf)
    return 0;
  for (i=0; i < tbf->myloc->numlocs; i++) {
    if (tbf->scores[i] > cutoff) {
      tbf->seq1endpnt[i] = (int) tbf->myloc->locs[i] + wordlen - 1;
      tbf->seq2endpnt[i] = tbf->location +wordlen - 1;
      reScoreHit(tbf, i, seq1->lets, seq2->lets, seq1->numlets, seq2->numlets, &first, &last);
      j = tbf->scores[i];
      if (tbf->scores[i] > reScoreCutoff){
	saveScore(tbf,i, first, last);
      }
    }
  }
  for (i=0; i < tbf->myloc->numlocs; i++) {
    if (tbf->mysles[i]) {
      remElem(tbf,i);
    }
    free (tbf->myhits[i].inds1);
    free (tbf->myhits[i].inds2);
  }

  next = tbf->next;

  free (tbf->myhits);
  free (tbf->scores);
  free (tbf->mysles);
  free (tbf->seq1startpnt);
  free (tbf->seq2startpnt);
  free (tbf->seq1endpnt);
  free (tbf->seq2endpnt);
  free (tbf->toberemoved);
  free (tbf);
  return savenfreeLList(next, seq1, seq2);
}

int mergeOverlap(hll* h1, hll* h2, seq* seq1, seq* seq2) {
  int offset, myscore, nextscore, newscore, bestloc, dir, gappen;
  int s1l, s2l, s1r, s2r, s1n, s2n;

  //  return 0;
  //  fprintf (stderr, "(%d %d) (%d %d)", h1->seq1end, h1->seq2end, h2->seq1start, h2->seq2start);

  if ((h1->seq2end < h2->seq2start) && (h1->seq1end < h2->seq1start)) {
    //    fprintf (stderr, " no\n");
    return 0;
  }
  
  offset = (h1->seq1end-h1->seq2end) - (h2->seq1start-h2->seq2start);
  if (ABS(offset) > mgaplen)
     return 0;
  gappen = gappenstart + gappenext * ABS(offset);

  if ((-gappen) > h1-> score || (-gappen) > h2->score) {
    //    fprintf (stderr, " gap\n");
    return 0;
  }
  s1l = h1->seq1end - h1->last->length;
  s2l = h1->seq2end - h1->last->length;
  s1r = h2->seq1start + h2->first->length;
  s2r = h2->seq2start + h2->first->length;

  if (s1r <= s1l || s2r <= s2l) {
    //    fprintf (stderr, " swap\n");
    return 0;
  }
  if (offset) {
    bestloc =  extendMerge(s1l, s2l, s1r, s2r, seq1->lets, seq2->lets, &dir);
    myscore = reScore(s1l, s2l, bestloc,  seq1->lets, seq2->lets);
    if (dir) {
	s1n = s1l + bestloc + ABS(offset)+1;
	s2n = s2l + bestloc + 1;
      }
      else {
	s2n = s2l + bestloc + ABS(offset)+1;
	s1n = s1l + bestloc + 1;
      }
    nextscore = reScore(s1n, s2n, s2r - s2n,  seq1->lets, seq2->lets);
    //    fprintf (stderr, " %d %d %d\n", bestloc, myscore, nextscore);    
    //    fprintf (stderr, "a %d %d %d\n", s1l, s1n, s1r);
    newscore = h1->score + h2->score - (h2->first->score -  nextscore) - (h1->last->score - myscore) + gappen;
    if (newscore < h1-> score || newscore < h2->score) {
      //      fprintf (stderr, " score1\n");
      return 0;
    }
    h1->score = newscore;
    h1->last->length = bestloc;

    h2->first->score = nextscore;
    h2->first->offset = offset;
    h2->first->length = s2r - s2n;
    h1->last->score = myscore;
    h1->last->next = h2->first;
    if (h1->last->next)
      h1->last = h2->last;
    h2->first = 0;
  }
  else {
    myscore = reScore(s1l, s2l, s1r-s1l,  seq1->lets, seq2->lets);
    newscore = h1->score + h2->score - (h1->last->score - myscore) + gappen;
    if (newscore < h1-> score || newscore < h2->score) {
      //      fprintf (stderr, " score2\n");
      return 0;
    }
    h1->score = newscore;
    h1->last->score = myscore;
    h1->last->next = h2->first->next;
    h1->last->length = s1r - s1l;
    if (h1->last->next)
      h1->last = h2->last;
    h2->first->next = 0;
  }
  h1->seq2end = h2->seq2end;
  h1->seq1end = h2->seq1end;
  return 1;
}
