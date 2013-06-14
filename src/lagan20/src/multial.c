#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include "diagmatrix.h"
#include "multial.h"

#define INSERTION 1
#define DELETION 2
#define BOTH 3

#define MISMATCH_CUTOFF 8
#define ANCHOR_LENGTH_CUTOFF 10
#define ANCHOR_SCORE_CUTOFF 1500

#define MAX_SQ_SIZE (100 * (1 << 20))
#define BIG_SQ_WIDTH 20

#define CONS_FRAC 0.6

#define MIN2(x,y)   ( (x) >= (y) ? (y) : (x) )
#define MAX2(x,y)   ( (x) >= (y) ? (x) : (y) )
#define MAX3(x,y,z)  MAX2(MAX2(x,y),z)
#define MIN3(x,y,z)  MIN2(MIN2(x,y),z)
#define PROD(x,y)   ( (x) * (y) )

#define WEQ2(x,y,a)  (((x)==(a))? 0: ((y)==(a))? 1:-1)
#define WEQ3(x,y,z,a)  (((x)==(a))? 0: ((y)==(a))? 1: ((z)==(a))? 2:-1)

char* alpha = "ATCG.N";
char* nucmatrixfile = 0;

int s1start = 0;
int s1end = 0;
int s2start = 0;
int s2end = 0;
//int match = 18;
//int mismatch = -8;
int gapstart = -50;
int gapend = -50;
int gapcont = -5;
int gapperseq = -1;
int overlap = 0;
int glwidth= 15;
char dobin = 0;

float factor, offset;
int logs[MAX_SEQ*MAX_SEQ];

FILE* outfile;

static int substmatrix[256][256];
static int matchcache[1 << 24], gapcache[1 << 24];
int *freed = 0, freedsize, freedcap;
align **freedptr;

int normf;
int normprev;

inline int ismatch(char a, char b) {
  return (a == b);
}

inline int isGap(align* ali, int seqn, int loc) {
  int i = !((ali->algn[loc] >> seqn) & 1);
  return i;
}

inline int scoreLocal(int which, align* ali, int loc) {
  int i, lets = 0;
  for (i=0; i < 4; i++)
    lets += ali->cnts[i][loc];
  //  printf ("which is %d lets is %d, cnts[w] is %d \n",which, lets, ali->cnts[which][loc]);

  if (which <4)
    return (ali->cnts[which][loc]-1) * 100 + (lets - ali->cnts[which][loc]) * -70 +
      ali->cnts[CNTS_GS][loc] * gapstart + ali->cnts[CNTS_GC][loc] * gapcont;
  if (which == CNTS_GS)
    return lets * gapstart;
  if (which == CNTS_GC)
    return lets+ali->cnts[CNTS_GS][loc] * gapcont;
}

inline hll* reverseHLL(hll* tbr) {
  hll *nn, *prev=0;
  while (tbr) {
    nn = tbr->next;
    tbr->next = prev;
    prev = tbr;
    tbr = nn;
  }
  return prev;
}

hll* getAnchsFromAlign(align* current, int seqnum, int cutoff) {
  int i=0, j, newj=0;
  int currscore=0, oldscore, peakscore;
  hll *res = 0, *temp = (hll*) malloc (sizeof(hll));
  int which;
  long long int mask = ~(1<<seqnum);
  char ingap = 0, isfrst = 1;
  float peakfrac;

  assert (temp);

  for (j = 0; j < current->algnlen; j++) {
    if (!isGap(current, seqnum, j)) {
      ingap = 0;
      which = strchr(alpha, current->seqs[seqnum]->lets[i]) - alpha;
      which = (which>3)?CNTS_LEN:which;
      i++;
    }
    else {
      if (ingap)
	which = CNTS_GC;
      else {
	ingap = 1;
	which = CNTS_GS;
      }
    }


    currscore += scoreLocal(which, current, j);

    if (currscore > cutoff) {
      temp->score = currscore;
      temp->seq1end = newj;  temp->seq2start = i;
      temp->seq2end = i; temp->seq1start = newj; 
      currscore = 0;
      temp->next = res; res = temp;temp = (hll*) malloc (sizeof(hll));
      assert (temp);
    }
    if (currscore < 0)
      currscore = 0;
    if (current->algn[j]&mask)
      newj++;
  }

  if (currscore > cutoff) {
    temp->score = currscore;
    temp->seq1end = newj;  temp->seq2start = i;
    temp->seq2end = i; temp->seq1start = newj;
    temp->next = res; res = temp;
  }
  else free(temp);
  return reverseHLL(res);
}

int cons_cnt = 0;


seq* mkConsensus(align* ali) {
  int i, j;
  seq* res = (seq*) malloc (sizeof(seq));
  assert (res);
  res->name = (char*) malloc(sizeof(char)*64);
  assert (res->name);
  sprintf(res->name, "Consensus_%d", ++cons_cnt);
  res->numlets = ali->algnlen;
  res->rptr = res->lets = (char*) malloc (sizeof(char) * res->numlets);
  assert (res->lets);
  for (i=0; i< res->numlets; i++) {
    res->lets[i] = 'N';
    for (j=0; j< 4; j++) {
      if (ali->cnts[j][i] >= ((float)ali->numseq) * CONS_FRAC)
	res->lets[i] = alpha[j];
    }
  }
  return res;
}

inline void reverse (long long int* a, int length) {
  long long int lft;
  int i;
  for (i=0; i < length/2; i++) {
    lft = a[i];
    a[i] = a[length-i-1];
    a[length-i-1] = lft;
  }
}


align* unifyAlign(align* ali1, align* ali2, align* uni){
  char *mat[MAX_SEQ];
  int i,j,k, cbc, brcount;
  int s1 = 0, s2 = 0, tgs, tgc;
  align *res = (align*) malloc(sizeof(align));
  
  assert (res);
  res->score = uni->score;
  res->numseq = ali1->numseq + ali2->numseq;
  res->algnlen = uni->algnlen;
  res->nextalign = 0;
  res->dirty = 0;

  // memory allocation and alignment creation
  res->algn = (long long int*) malloc ((res->algnlen+1) * sizeof (long long int)); assert (res->algn);
  res->algn[0] = 0;
  for (j = 0; j < CNTS_LEN; j++){
    res->cnts[j] = (char*) malloc((res->algnlen+1) * sizeof(char));
    assert (res->cnts[j]);
  }
  for (i=0; i<= res->algnlen; i++){
    res->algn[i] = 0;
    for (j=0; j<CNTS_LEN; j++)
      res->cnts[j][i] = 0; 
    if (!isGap(uni, 0, i)) res->algn[i] |= ali1->algn[s1++];
    if (!isGap(uni, 1, i)) res->algn[i] |= (ali2->algn[s2++] << ali1->numseq);
  }

  for (i = 0; i < res->numseq; i++){
    res->seqs[i] = (i < ali1->numseq) ? ali1->seqs[i] : ali2->seqs[i - ali1->numseq];
    mat[i] = (char *) malloc (sizeof (char) * (res->algnlen + 1)); assert (mat[i]);
    mat[i][0] = 0;
    for (j = 0, k = 0; j <= res->algnlen; j++)
      mat[i][j] = isGap (res, i, j) ? '-' : res->seqs[i]->lets[k++];
  }  

  s1 = s2 = 1;
  
  for (i=0; i<=res->algnlen; i++){
    for (j = 0; j < res->numseq; j++){
      switch (mat[j][i]){
      case 'A': res->cnts[CNTS_A][i]++; if (i > 1 && mat[j][i-1] == '-') res->cnts[CNTS_GE][i]++; break;
      case 'T': res->cnts[CNTS_T][i]++; if (i > 1 && mat[j][i-1] == '-') res->cnts[CNTS_GE][i]++; break;
      case 'C': res->cnts[CNTS_C][i]++; if (i > 1 && mat[j][i-1] == '-') res->cnts[CNTS_GE][i]++; break;
      case 'G': res->cnts[CNTS_G][i]++; if (i > 1 && mat[j][i-1] == '-') res->cnts[CNTS_GE][i]++; break;
      case '-':
	if (i > 0 && mat[j][i-1] == '-')
	  res->cnts[CNTS_GC][i]++;
	else
	  res->cnts[CNTS_GS][i]++;
	break;
      }
    }
  }
  
  for (i = 0; i < res->numseq; i++) free (mat[i]);

  return res;
}


align* getChain(dmat* mydm, int x, int y, int j) {
  int temp;
  align *res = (align*) malloc (sizeof(align)), *help; 
  long long int* almt = (long long int*) malloc ( sizeof(long long int));
  int i=0, almtsize = 1, which, inrun = j;
  char zz = DMgetPtr(mydm, x, y); 
  assert (res);
  assert (almt);
  
  for (i=0; i<CNTS_LEN; i++)
    res->cnts[i] = 0;
  i = 0;

  ///////////////
  res->dirty = 0;
  res->nextalign = 0;
  res->algn = 0;
  res->algnlen = 0;

  res->num = freedsize;
  freed[freedsize] = 0;
  freedptr[freedsize] = res;
  if (++freedsize == freedcap){
    freedcap *= 2;
    freed = (int *) realloc (freed, sizeof (int) * freedcap);
    freedptr = (align **) realloc (freedptr, sizeof (align *) * freedcap);
  }

  do { 
    //    printf("I am at %d,%d  %x\n", x,y, zz);
    which = zz & Mmask;

    if (which == 0x3) {
      help = DMgetNeck(mydm, x, y, inrun);
      if (!help) {
	if (i > 2)
	  fprintf (stderr, "PROBLEM %d %d after %d (norm %d, %d)\n", x, y,i, normf, normprev);
	free(almt);
	res->algn = 0;
	res->algnlen = i;
	return res;
      }
      /*      if (! help->nextalign)
	fprintf (stderr, "check %d %d after %d\n", x, y,i);
      */
      help->dirty++;
      res->nextalign = help;
      break;
    }

    
    if (inrun == 1 && (zz & Nmask))
      which = 1;
    else if (inrun == 2 && (zz & Omask))
      which = 2;
    else
      which = 0;
    
    
    /*
    if (inrun == 1) {
      if (zz & Nmask) {
	which = 1;
      }
    }
    else if (inrun == 2) {
      if (zz & Omask) {
	which = 2;
      }
    }
    */

    if (which == 0) {
      inrun = zz & Mmask;
      almt[i++] = BOTH;
      zz = DMgetPtr(mydm,--x,--y);
    }

    else if (which == 1) {  /*N*/
      inrun = 1;
      almt[i++] = INSERTION;
      zz = DMgetPtr(mydm, --x, y);
    }
    
    else if (which == 2) {
      inrun = 2;
      almt[i++] = DELETION;
      zz = DMgetPtr(mydm, x, --y);
    }
    else 
      printf("a really dumb error %d\n", i);
 
    if (i >= almtsize) {
      almt = realloc (almt, sizeof(long long int)* (almtsize *= 2));
    }
    //   printf ("retrace %d %d after %d\n", x, y,i);

  } while (x > 0 && y > 0);
    reverse(almt, i);

  //  fprintf(stderr, "getChain done at %d %d after %d\n", x , y , i);
  //  printf("gotChain\n");
  res->algn = almt;
  res->algnlen = i;
  //  printf("done w it\n");
  return res;
}


void saveNeck(dmat* mydm, int neckdiag) {
  int size1, size2, x1, x2, y1, y2;
  alel *first = DMgetDiagStart(mydm, neckdiag-1, &size1, &x1, &y1),
    *second = DMgetDiagStart(mydm, neckdiag, &size2, &x2, &y2);
  int i, j;
  align* a;

  //  printf("saving neck %d\n", neckdiag);
  normprev = normf;
  normf = DMnextNecks(mydm, neckdiag);

  for (i=0; i<size2; i++,x2++,y2--) {
    for (j=0; j<3; j++) {
      a = getChain(mydm, x2, y2, j);
      DMsetNeck(mydm, a, x2, y2, j);
    }
  }
  for (i=0; i<size1; i++,x1++,y1--) {
    for (j=0; j<3; j++) {
      a = getChain(mydm, x1, y1, j);
      DMsetNeck(mydm, a, x1, y1, j);
    }
  }
}

void joinAligns (align* a) {
  align *n = a->nextalign, *t;
  long long int* temp,  *temp2;
  int totsize=0;
  int i =0;
  for (t = a; t; t = t->nextalign) {
    totsize += t->algnlen;
    i++;
  }

  temp = malloc ((totsize+1)*sizeof(long long int));
  assert (temp);
  temp[totsize] = 0;
  temp2 = temp + totsize;
  totsize = 0;
  for (t=a; t; t = t->nextalign) {
    totsize += t->algnlen;
    memcpy(temp2-totsize, t->algn, t->algnlen*sizeof(long long int));
  }
  free (a->algn);
  a->algn = temp;
  a->algnlen = totsize;
  a->nextalign = 0;
  /*
  for (a = a->nextalign; a;) {
    t = a;
    a = a->nextalign;
    freeAlign(t);
  }
  */
}

inline int scoreGap(int numgs, int numgc, int numge, int numseq) {
  return (MIN2(numgc, numseq-numgc) * gapcont) +
    (MIN2(numgs, numseq-numgs) * gapstart) +
    (MIN2(numge, numseq-numge) * gapend);
}

void printcache(){
  int a, b, c, d;
  for (a = 0; a < 3; a++){
    for (b = 0; b < 3; b++){
      for (c = 0; c < 3; c++){
	for (d = 0; d < 3; d++){
	  fprintf (stderr, "%d %d %d %d -- %d\n", a, b, c, d, matchcache[a | (b << 6) | (c << 12) | (d << 18)]);
	}
      }
    }
  }
}

char getLetter (FILE *file){
  char ch;

  while (!feof (file)){
    ch = fgetc (file);
    if (!isspace (ch)){
      //      fprintf (stderr, "LETTER READ: \"%c\"\n", ch);
      return ch;
    }    
  }

  assert (0);
  return 0;
}

int readit = 0;

void readSubstMatrix (char *filename, int size, int substmatrix[256][256]){
  FILE *file;
  char line[1024];
  unsigned char *symbs, ch;
  int i, j, k;

  if (readit) return;
  readit = 1;

  if (!nucmatrixfile) {
    sprintf (line, "%s/%s", getenv ("LAGAN_DIR"), filename);
    file = fopen (line, "r"); assert (file);
  }
  else {
    file = fopen (nucmatrixfile, "r"); assert (file);
    
  }

  for (i = 0; i < 256; i++){
    for (j = 0; j < 256; j++){
      substmatrix[i][j] = 0;
    }
  }
  
  symbs = (unsigned char *) malloc (sizeof (unsigned char) * size); assert (symbs);
  for (i = 0; i < size; i++) symbs[i] = (unsigned char) getLetter (file);
  for (i = 0; i < size; i++){
    ch = getLetter (file);
    assert (ch == symbs[i]);
    for (j = 0; j < size; j++){
      fscanf (file, "%d", &k);
      //      fprintf (stderr, "NUMBER READ: %d\n", k);
      substmatrix[(int) symbs[i]][(int) symbs[j]] = k;
      assert ((int) symbs[i] > 0);
      assert ((int) symbs[j] > 0);
    }
  }

  fscanf (file, "%d", &gapstart);
  fscanf (file, "%d", &gapcont);
  //  fprintf (stderr, "GAP SCORES: %d %d\n", gapstart, gapcont);
  gapend = gapstart / 2;
  gapstart -= gapend;
  
  free (symbs);
  fclose (file);
}

inline int chmatchscore (unsigned char a, unsigned char b, int substmatrix[256][256]) {
  return substmatrix[a][b];
}

void buildcache (){
  int score, i, j;
  int gs, gc, ge, ns;
  char *lets = "ATCG";
  int num[4];
  int numseqs = MAX_SEQ;

  readSubstMatrix (NUC_FILE, NUC_FILE_SIZE, substmatrix);

  for (num[0] = 0; num[0] <= numseqs; num[0]++){ // A
    for (num[1] = 0; num[1] <= numseqs; num[1]++){ // T
      for (num[2] = 0; num[2] <= numseqs; num[2]++){ // C
	for (num[3] = 0; num[3] <= numseqs; num[3]++){ // G

	  score = 0;
	  for (i = 0; i < 4; i++){
	    score += num[i] * (num[i] - 1) / 2 * chmatchscore ((unsigned char)lets[i], (unsigned char)lets[i], substmatrix);
	    for (j = i + 1; j < 4; j++){
	      score += num[i] * num[j] * chmatchscore ((unsigned char) lets[i], (unsigned char) lets[j], substmatrix);
	    }
	  }
	  matchcache[num[0] | (num[1] << 6) | (num[2] << 12) | (num[3] << 18)] = score;
	}
      }
    }
  }

  for (gs = 0; gs <= numseqs; gs++){
    for (gc = 0; gc <= numseqs; gc++){
      for (ge = 0; ge <= numseqs; ge++){
	for (ns = 0; ns <= numseqs; ns++){
	  gapcache[gs | (gc << 6) | (ge << 12) | (ns << 18)] = scoreGap (gs, gc, ge, ns);
	}
      }
    }
  }

  //  builtcache = 1;

  // printcache();
}

inline int v (int y){
  if (y >= 0 && y <= MAX_SEQ) return y;
  fprintf(stderr, "Got %d in v\n", y);
  assert (0);
  return 0;
}

inline int matchscore (align*a, int ai, align *b, int bi){
  
  return
    matchcache[v(a->cnts[0][ai] + b->cnts[0][bi]) | 
	      (v(a->cnts[1][ai] + b->cnts[1][bi]) << 6) |
	      (v(a->cnts[2][ai] + b->cnts[2][bi]) << 12) |
	      (v(a->cnts[3][ai] + b->cnts[3][bi]) << 18)] +
    gapcache[v(a->cnts[CNTS_GS][ai] + b->cnts[CNTS_GS][bi]) |
	    (v(a->cnts[CNTS_GC][ai] + b->cnts[CNTS_GC][bi]) << 6) |
	    (v(a->cnts[CNTS_GE][ai] + b->cnts[CNTS_GE][bi]) << 12) |
	    (v(a->numseq + b->numseq - (a->cnts[CNTS_CB][ai] + b->cnts[CNTS_CB][bi])) << 18)];
}

inline int scoreOpp (align *other, int ow, int oppnum){
  return matchcache[v(other->cnts[0][ow]) | 
		   (v(other->cnts[1][ow]) << 6) |
		   (v(other->cnts[2][ow]) << 12) |
		   (v(other->cnts[3][ow]) << 18)];
}

inline int endGap0 (align* a, int ai, align* b, int bi){
  return gapcache[(v(a->cnts[CNTS_GE][ai]+b->cnts[CNTS_GE][bi])<<12) | 
		  (v(a->numseq + b->numseq-(b->cnts[CNTS_CB][bi]+a->cnts[CNTS_CB][ai])) << 18)];
}

inline int endGap1 (align* a, int ai, align* b, int bi){

  return gapcache[(v((b->numseq - b->cnts[CNTS_GS][bi] - b->cnts[CNTS_GC][bi]) + a->cnts[CNTS_GE][ai]) << 12) | 
		  (v(a->numseq + b->numseq - (b->cnts[CNTS_CB][bi]+a->cnts[CNTS_CB][ai])) << 18)];
}

inline int endGap2 (align* a, int ai, align* b, int bi){
  return gapcache[(v((a->numseq - a->cnts[CNTS_GS][ai] - a->cnts[CNTS_GC][ai]) + b->cnts[CNTS_GE][bi])<<12) | 
		  (v(a->numseq + b->numseq - (b->cnts[CNTS_CB][bi]+a->cnts[CNTS_CB][ai])) << 18)];
}

inline int contGap(align* ali, int myw, align* other, int ow, int *sopp) {
  return gapcache[(v(other->cnts[CNTS_GS][ow])) |
		  (v(ali->numseq + other->cnts[CNTS_GC][ow]) << 6) |
		  (v(other->cnts[CNTS_GE][ow]) << 12) |
		  (v(ali->numseq + other->numseq - (ali->cnts[CNTS_CB][myw] + other->cnts[CNTS_CB][ow])) << 18)] +
    sopp[ow];
}

inline int openGap(align* ali, int w, align* other, int ow, int *sopp, char *desc) {
  int alopen, pen, sav, i;

  alopen = ali->cnts[CNTS_GC][w] + ali->cnts[CNTS_GE][w];
  /**
   * Watch out for running off end of array.
   */
  //  if (w < ali->algnlen) alopen += ali->cnts[CNTS_GS][w+1];

  
  sav = gapcache[(v(ali->numseq - (alopen + ali->cnts[CNTS_CB][w]) + other->cnts[CNTS_GS][ow])) |
		 (v(alopen + other->cnts[CNTS_GC][ow]) << 6) |
		 (v(other->cnts[CNTS_GE][ow]) << 12) |
		 (v(ali->numseq+other->numseq - (ali->cnts[CNTS_CB][w]+other->cnts[CNTS_CB][ow])) << 18)];

  return sav;
}


void mkBarrel(int s1, int s2, int e1, int e2, int width, int *dn, int dt, int* starts, int *ends, dmat* mydm) {
  int sd = s1+s2-1, dlen;
  int elem = (sd < mydm->d2)? s1: mydm->d2-s2;
  int incr;
  double fl = 0;
  double slope = (double)(e2-s2)/(double)(e1-s1);
  double cloc = elem;

  if ((e2-s2 == 0) && (e1-s1 == 0))
    slope = 1;
  else if (e1-s1 == 0)
    slope = 100000;
  //  // printf("dt = %d\n", dt);
  //  printf("BA: %d, %d to %d, %d %f\n", s1,s2,e1,e2,slope);
  for ( ; sd <(*dn); sd++) {
    if (fl>=slope || (int)(cloc) == (int)(cloc+slope)) {
      cloc+=slope;
      fl -= slope;
    }
    else {
      elem--;
      fl++;
    }
    if (sd <= mydm->d2)     
      elem++;
  }
  fl = 0;
  for ( ; *dn < dt; (*dn)++) {
    //    // printf("dn =%d  ", *dn);
    if (fl>=slope || (int)(cloc) == (int)(cloc+slope)) {
      cloc+=slope;
      fl -= slope;
    }
    else {
      elem -=1;
      fl++;
    }
    if (*dn <= mydm->d2) 
      elem++;

    if (*dn < MIN2(mydm->d2, mydm->d1))
      dlen = *dn;
    else if (*dn < MAX2(mydm->d2, mydm->d1))
      dlen = MIN2(mydm->d2, mydm->d1);
    else 
      dlen = mydm->d2 + mydm->d1 - *dn;
    starts[*dn] = MAX2(elem - width, 0);
    ends[*dn] = MIN2(elem+width, dlen-1);
  }
}



void mkSquare(int s1, int s2, int e1, int e2, int *dn, int dt, int* starts, int *ends, dmat* mydm) {
  int dists[2], dlen;
  long long int size = ((long long int)e1-(long long int)s1)
    * ((long long int)e2-(long long int)s2);
  int dn2;
  int eval, sval;
  
  if (size > MAX_SQ_SIZE) {
    fprintf (stderr, "SQUARE TOO BIG: %d,%d to %d,%d\n", s1, e1,s2,e2);
    mkSquare(s1, s2, (s1+e1)/2+glwidth, (s2+e2)/2+glwidth, dn, (*dn+dt)/2, starts, ends, mydm);
    mkSquare((s1+e1)/2-glwidth, (s2+e2)/2-glwidth, e1, e2, dn, dt, starts, ends, mydm);
    return;
  }
  //  // printf("dt = %d\n", dt);
  //  // printf("SQ: %d, %d to %d, %d\n", s1,s2,e1,e2);

  // fill in part before square
  dn2 = *dn - 1;
  while (1){
    if (dn2 < mydm->d2) {
      dists[0] = s1-1;
      dists[1] = dn2 - e2;
    }
    else {
      dists[0] = mydm->d2 - e2;
      dists[1] = s1 - (dn2 - mydm->d2)-1;
    }
    starts[dn2] = MIN2(starts[dn2], sval = MAX3(dists[0], dists[1],0));

    if (dn2 < mydm->d2) {
      dists[0] = e1-1;
      dists[1] = dn2 - s2;
    }
    else {
      dists[0] = mydm->d2 - s2;
      dists[1] = e1 - (dn2-mydm->d2)-1;
    }
    if (dn2 < MIN2(mydm->d2, mydm->d1))
      dlen = dn2;
    else if (dn2 < MAX2(mydm->d2, mydm->d1))
      dlen = MIN2(mydm->d2, mydm->d1);
    else 
      dlen = mydm->d2 + mydm->d1 - dn2;
    ends[dn2] = MAX2(ends[dn2], eval = MIN3(dists[0], dists[1],dlen-1));
    if (eval - sval <= 5) break; // break after fill in
    dn2--;
  }

  for ( ; *dn < dt; (*dn)++) {
    //    // printf("square dn = %d\n", *dn);
    if (*dn < mydm->d2) {
      dists[0] = s1-1;
      dists[1] = *dn - e2;
    }
    else {
      dists[0] = mydm->d2 - e2;
      dists[1] = s1 - (*dn - mydm->d2)-1;
    }
    starts[*dn] = MAX3(dists[0], dists[1],0);

    if (*dn < mydm->d2) {
      dists[0] = e1-1;
      dists[1] = *dn - s2;
    }
    else {
      dists[0] = mydm->d2 - s2;
      dists[1] = e1 - (*dn-mydm->d2)-1;
    }
    if (*dn < MIN2(mydm->d2, mydm->d1))
      dlen = *dn;
    else if (*dn < MAX2(mydm->d2, mydm->d1))
      dlen = MIN2(mydm->d2, mydm->d1);
    else 
      dlen = mydm->d2 + mydm->d1 - *dn;
    ends[*dn] = MIN3(dists[0], dists[1],dlen-1);
  }
}

void doShapes(hll* myres, dmat* mydm, int* starts, int *ends) {
  int p1=MAX2(overlap,glwidth)+1, p2=MAX2(overlap,glwidth)+1; 
  int t1, t2;
  int dn = 1, dt;
  int width = glwidth;
  while (myres) {

    while (1){
      if (!myres || (myres->seq1start >= 1 && myres->seq2start >= 1 &&
		     myres->seq1end >= 1 && myres->seq2end >= 1 &&
		     myres->seq1start < mydm->d1 && myres->seq2start < mydm->d2 &&
		     myres->seq1start < myres->seq1end && myres->seq2start < myres->seq2end &&
		     myres->seq1end < mydm->d1 && myres->seq2end < mydm->d2 &&
		     abs((myres->seq1end-myres->seq1start) -
			 (myres->seq2end-myres->seq2start)) <= MISMATCH_CUTOFF))
	break;
      myres = myres->next;
    }
    if (!myres) break;

    /*
    printf("--> (%d %d)=(%d %d)\n", 
	   myres->seq1start, myres->seq1end,
	   myres->seq2start, myres->seq2end);
    */
    t1 = myres->seq1start;   /* between hits */
    t2 = myres->seq2start;
    dt = t1 + t2 - 1 + overlap;    
    mkSquare(p1-MAX2(overlap, width), p2-MAX2(overlap, width), 
	     t1+MAX2(overlap, width), t2+MAX2(overlap, width), 
	     &dn, dt, starts, ends, mydm);
    p1 = myres->seq1end;   /* within a hit */
    p2 = myres->seq2end;
    dt = p1 + p2 - 1 - overlap; 
    mkBarrel(t1, t2, p1, p2, width, &dn, dt, starts, ends, mydm);
    myres = myres->next;
  }
  t1 = mydm->d1; 
  t2 = mydm->d2; 
  dt = t1 + t2;     
  mkSquare(p1-MAX2(overlap,width), p2-MAX2(overlap,width), t1, t2, &dn, dt, starts, ends, mydm);
}


void doAncs(dmat* mydm, align* ali1, align* ali2, hll* ancs) {
  int *starts, *ends;

  starts = (int*) malloc(sizeof(int)*(ali1->algnlen + ali2->algnlen+2)); assert (starts);
  ends = (int*) malloc(sizeof(int)*(ali1->algnlen + ali2->algnlen+2)); assert (ends);
  doShapes(ancs, mydm, starts, ends);
  DMinitDiag(mydm, starts,ends);
  free(starts);
  free(ends);
}


align* doNW(dmat* mydm, align* ali1, align* ali2) {
  int i, j;
  int x, y, size;
  int gapstartN = 0, gapstartO = 0;
  int gapcontN, gapcontO; 
  int gapend[3];
  int tt, prevgap;
  alel *curr, *pasts0, *pasts1, *pasts2; 
  align* a, *b;
  char rh, ptr=0, isneck;
  int ndiags = mydm->d1 + mydm->d2 -1;
  int *sopp1, *sopp2;
  int numNecks =0, oldneck =0;
  register int s1, s2, s3, z1, z2,z3;

  //  int M[20][20][6];

  
  isneck = DMnextDiag(mydm);
  curr = DMgetDiagStart(mydm, 1, &size, &x, &y);
  curr->N = curr->O = 0;
  curr->M = 0;
  DMsetPtr(mydm, 0, 1, 1);

  buildcache();

  sopp1 = (int*) malloc (sizeof (int) * (ali1->algnlen+1));
  sopp2 = (int*) malloc (sizeof (int) * (ali2->algnlen+1));
  assert (sopp1); assert (sopp2);

  for (i = 0; i < ali1->algnlen; i++) sopp1[i] = scoreOpp (ali1, i, 0);
  for (i = 0; i < ali2->algnlen; i++) sopp2[i] = scoreOpp (ali2, i, 0);

  /*fprintf (stderr, "Checking diagonals...\n");
  for (i = ndiags - 50; i <= ndiags; i++){
  DMgetDiagStart (mydm, i, &size, &x, &y); */

  //  fprintf (stderr, "ndiag = %d (%d %d)\n", ndiags, ali1->algnlen, ali2->algnlen);
 
  for (i = 2; i <= ndiags; i++) {
    isneck = DMnextDiag(mydm);
    if (!(i%10000))
      fprintf(stderr, "WORKING %d/%d\n", i/10000,ndiags/10000 );
    
    curr = DMgetDiagStart(mydm, i, &size, &x, &y);
    pasts2 = DMgetElem(mydm, x-1, y);
    pasts1 = DMgetElem(mydm, x-1, y-1);

    for (j = 0; j < size; j++) {
      gapstartN = openGap(ali2, y, ali1, x, sopp1, "gapstartN");
      gapstartO = openGap(ali1, x, ali2, y, sopp2, "gapstartO");

      gapcontN = contGap(ali2, y, ali1, x-1, sopp1);
      gapcontO = contGap(ali1, x, ali2, y-1, sopp2);

      pasts0 = pasts2;
      pasts2 = DMgetElem2(mydm, x, y-1, pasts2);

      curr->M = matchscore (ali1, x - 1, ali2, y - 1);

      z1 = pasts1->M + endGap0 (ali1, x - 1, ali2, y - 1);
      z2 = pasts1->N + endGap1 (ali1, x - 1, ali2, y - 1);
      z3 = pasts1->O + endGap2 (ali1, x - 1, ali2, y - 1);

      if (z1 >= z2){
	if (z1 >= z3){ curr->M += z1; ptr = 0; }// + endGap0 (ali1, x - 0, ali2, y - 0); }
	else         { curr->M += z3; ptr = 2; }// + endGap2 (ali1, x - 0, ali2, y - 0); }
      }
      else {
	if (z2 >= z3){ curr->M += z2; ptr = 1; } // + endGap1 (ali1, x - 0, ali2, y - 0); }
	else         { curr->M += z3; ptr = 2; } // + endGap2 (ali1, x - 0, ali2, y - 0); }
      }

      s2 = pasts0->N + gapcontN;
      s3 = pasts2->O + gapcontO;

      s1 = curr->M + gapstartN;
      if (s1 >= s2){ curr->N = s1; }
      else         { curr->N = s2; ptr |= 4; }
      s1 = curr->M + gapstartO;
      if (s1 >= s3){ curr->O = s1; }
      else         { curr->O = s3; ptr |= 8; }

      DMsetPtr(mydm, ptr, x, y);

      curr++; x++; y--;

      pasts1 = DMgetElem2(mydm, x-1, y-1, pasts1);
    }
    if (isneck) {
      numNecks++;
      saveNeck(mydm, i);
      oldneck = i;
    }
  }
  
  free (sopp1);
  free (sopp2);

  mydm->currneck++;
  a = getChain(mydm, mydm->d1, mydm->d2, 0);
  curr--;
  a->score = MAX3(curr->M, curr->N, curr->O);
  freed[a->num] = 1;  
  joinAligns(a);



  //  fprintf(stderr, "done NW\n");
  return a;
}

align* makeAlign(align* ali1, align* ali2, hll* anchors, align **uni) {
  align *res;
  dmat* mydm;
  int numseq = ali1->numseq + ali2->numseq, i;
  int oldgapstart = gapstart, oldgapcont = gapcont, oldgapend = gapend;

  mydm = makeDM(ali1->algnlen, ali2->algnlen);

  gapstart *= (numseq-1); gapend *= (numseq-1); 
  gapcont *= (numseq-1);
  fprintf (stderr, "gs ge gc %d %d %d\n", gapstart, gapend, gapcont);
  //  initEntropy(ali1, ali2);

  doAncs(mydm, ali1, ali2, anchors);

  freedsize = 0; freedcap = 1;
  freed = (int *) malloc (sizeof (int) * freedcap);
  freedptr = (align **) malloc (sizeof (align *) * freedcap);
  assert (freed);
  assert (freedptr);

  *uni = doNW(mydm, ali1, ali2);
  res = unifyAlign(ali1, ali2, *uni);
  //  printf("firstlen = %d, seclen = %d, relen = %d\n", ali1->algnlen, ali2->algnlen, res->algnlen);
  freeDM(mydm);

  //  fprintf(stderr, "Final freeing\n");
  for (i = freedsize-1; i >= 0; i--){
    if (!freed[i]){
      freeAlign (freedptr[i]);
      freedptr[i] = 0;
    }
  }
  //  fprintf(stderr, "Final freeing done\n");
  free (freed); free (freedptr);
  freed = 0;
  gapstart = oldgapstart; gapend = oldgapend; gapcont = oldgapcont;
  
  return res;
}

align* mkSimAlign(seq* seq1) {
  int i,j,k,oldk=-1;
  align* res = (align*) malloc( sizeof(align));
  assert (res);

  res->score = 0;
  res->nextalign = 0;
  res->dirty = 0;
  res->numseq = 1;
  res->algnlen = seq1->numlets;
  res->seqs[0] = seq1;

  /**
   * Evidence that you need one more character.
   */
  res->algn = (long long int*) malloc((res->algnlen+1) * sizeof(long long int));
  assert (res->algn);
  for (j=0; j<CNTS_LEN; j++){
    res->cnts[j] = (char*) malloc((res->algnlen+1) * sizeof(char));    
    assert (res->cnts[j]);
  }
  for (i=0; i< res->algnlen;i++) {
    for (j=0; j<CNTS_LEN; j++)
      res->cnts[j][i] = 0; 
    res->algn[i] = 1;
    k=strchr(alpha,seq1->lets[i])-alpha;
    if (k<5)
      res->cnts[k][i]++;
    if (oldk == 4)
      res->cnts[4][i]++;
    oldk = k;
  }
  for (j=0; j<CNTS_LEN; j++)
    res->cnts[j][i] = 0; 
  res->algn[i] = 0;
  return res;
}

 
align* removeSeq(align* ali, int seqnum) {
  int i,j, k, n, p, bit = (1 << seqnum);
  int mask = bit - 1, resint, flag = 0;
  align* res = (align*) malloc(sizeof(align));
  res->score = 0;
  res->numseq = ali->numseq-1;
  for (i=0; i< seqnum; i++)
    res->seqs[i] = ali->seqs[i];
  for (i++; i< ali->numseq; i++)
    res->seqs[i-1] = ali->seqs[i];

     res->algn = (long long int*) malloc(ali->algnlen * sizeof(long long int));  
  for (j=0; j<CNTS_LEN; j++)
    res->cnts[j] = (char*) malloc(ali->algnlen * sizeof(char));    

  for (i=0, j=0, n=0; i < ali->algnlen; i++) {
    resint = (ali->algn[i] & mask) | ((ali->algn[i] & ~(mask|bit)) >> 1);
    if (resint) {
      for (k=0; k<CNTS_LEN; k++)
	res->cnts[k][j] = ali->cnts[k][i]; 
      res->algn[j] = resint;
      if (!isGap(ali, seqnum, i)) {
	k=strchr(alpha,ali->seqs[seqnum]->lets[n])-alpha;
	if (k<5)
	  res->cnts[k][j]--;
	if (i && isGap(ali, seqnum, i-1))
	  res->cnts[CNTS_GE][j]--;
	n++;
      }
      else {
	if (i && isGap(ali, seqnum, i-1))
	  res->cnts[CNTS_GC][j]--;
	else
	  res->cnts[CNTS_GS][j]--;
      }
      if (flag) {
	
	res->cnts[CNTS_GS][j] = 0;
	res->cnts[CNTS_GC][j] = 0;
	res->cnts[CNTS_GE][j] = 0;
	for (p = 0; p < res->numseq; p++) {
	  if (j<=1 || isGap(res, p, j-1)) {
	    if (!isGap(res, p, j))
	      res->cnts[CNTS_GE][j]++;
	    else
	      res->cnts[CNTS_GC][j]++;
	  }
	  else {
	    if (j && isGap(res, p, j))
	      res->cnts[CNTS_GS][j]++;
	  }
	}
      }
      j++;
    }
    else { n++; flag = 1;}
  }

  res->algnlen = j;

  for (i=0; i<CNTS_LEN; i++)
    res->cnts[i][j] = 0;

  //  printf("%d squished to %d\n", ali->algnlen, res->algnlen);
  return res;
}


align* removeSeqByName(align* ali, char *name) {
  int i=0;

  seq *removed;

  while (strcmp(ali->seqs[i]->name, name)) { i++; }
  removed = ali->seqs[i];

  removeSeq(ali, i);
}

int getSeqNum(align* ali, seq* trgt) {
  int i=0;

  seq *removed;

  while (ali->seqs[i] != trgt) { i++; }
  return i;
}


void swapHLL(hll* h1) {
  int i, j;
  
  while(h1) {
    i=h1->seq1start;
    j=h1->seq1end;
    h1->seq1start=h1->seq2start;
    h1->seq1end=h1->seq2end;
    h1->seq2start=i;
    h1->seq2end=j;
    h1=h1->next;
  }
}


int countpos (align* aln, int seqnum){
  int i, j = 0;
  for (i = 0; i < aln->algnlen; i++){
    if (!isGap (aln, seqnum, i)) j++;
  }
  return j;
}

hll* remapHLLs(hll* anchs, int which, align* aln, int seqnum) {
  int mybp, i, *searchint, stmybp, mylen, olen, osize;
  hll *wlist = anchs, *temp, *prev;
  float scale;
  char isfrst=1;

  // fprintf (stderr, "which=%d\n", which);
  //
  //    fprintf (stderr, "This is a list of the entries before going into remapHLLs:\n");
  //    printHLL (anchs);

  if (!anchs)
    return anchs;

  mylen = countpos (aln, seqnum);
  //    olen = countpos (aln, !seqnum);

  //   fprintf (stderr, "Here is some information about the alignment:\n");
  //   fprintf (stderr, "   alignment length = %d\n", aln->algnlen);
  //   fprintf (stderr, "   number of positions in sequence to remap = %d\n", mylen);
  //   fprintf (stderr, "   number of positions in other sequence = %d\n", olen);

  prev = NULL;
  for (temp = wlist; temp; temp = temp->next){
    if (temp->seq1start < 1) temp->seq1start = 1;
    if (temp->seq2start < 1) temp->seq2start = 1;
    if (!which && temp->seq1end > mylen) temp->seq1end = mylen;
    else if (which && temp->seq2end > mylen) temp->seq2end = mylen;

    if (temp->seq1start > temp->seq1end) {
      fprintf(stderr, "1 (%d %d)(%d %d)", temp->seq1start, temp->seq1end, temp->seq2start, temp->seq2end);
      assert(0);
    }

    if  (temp->seq2start > temp->seq2end) {
      fprintf(stderr, "2 (%d %d)(%d %d)", temp->seq1start, temp->seq1end, temp->seq2start, temp->seq2end);
      assert(0);
    }
  }

  wlist = (hll*)malloc(sizeof(hll)); assert (wlist);
  wlist->next = anchs;
  prev = wlist;

  mybp = stmybp = 0;
  searchint = (!which)?&(anchs->seq1start):&(anchs->seq2start);
  
  for (i=1; i<=aln->algnlen; i++) {
    if (isGap(aln,seqnum,i)){
      if (isfrst) continue;

      scale = (!which) ? 
	((anchs->seq1end == stmybp) ? 0 : (float)(mybp - stmybp) / (float)(anchs->seq1end - stmybp)) :
	((anchs->seq2end == stmybp) ? 0 : (float)(mybp - stmybp) / (float)(anchs->seq2end - stmybp));
      osize = (!which) ?
	(int)((anchs->seq2end - anchs->seq2start) * scale) :
	(int)((anchs->seq1end - anchs->seq1start) * scale);
      assert (osize >= 0);
      
      if (//mybp - stmybp < ANCHOR_LENGTH_CUTOFF || osize < ANCHOR_LENGTH_CUTOFF ||
	  anchs->score * scale < ANCHOR_SCORE_CUTOFF){
	  
	//	fprintf (stderr, "1. The region from %d to %d was cropped.\n", stmybp, mybp);

	if (!which){	  
	  anchs->score -= anchs->score * scale;
	  anchs->seq1start = mybp+1;
	  anchs->seq2start = anchs->seq2start + osize + 1;
	  isfrst = 1;
	  searchint = &(anchs->seq1start);
	}
	else {
	  anchs->score -= anchs->score * scale;
	  anchs->seq1start = anchs->seq1start + osize + 1;
	  anchs->seq2start = mybp+1;
	  isfrst = 1;
	  searchint = &(anchs->seq2start);
	}

	if (anchs->seq1start >= anchs->seq1end || anchs->seq2start >= anchs->seq2end){
	  //	  fprintf (stderr, "6. The region from %d to %d was thrown away.\n", stmybp, mybp);
	  temp = anchs;
	  prev->next = anchs->next;
	  anchs = anchs->next;
	  free (temp);
	  if (!anchs) break;
	  searchint = (!which)?&(anchs->seq1start):&(anchs->seq2start);
	}	
	continue;
      }

      temp = (hll*) malloc(sizeof(hll)); assert (temp);
      temp->next = anchs->next;
      anchs->next = temp;
      temp->seq1end = anchs->seq1end;
      temp->seq2end = anchs->seq2end;


      //      fprintf (stderr, "2. A new region from %d to %d was created.\n", stmybp, mybp);
      //fprintf (stderr, "Currently looking at (%d %d)=(%d %d)\n", anchs->seq1start, anchs->seq1end, anchs->seq2start, anchs->seq2end);


      if (!which){
	temp->score = anchs->score * scale;
	anchs->score -= temp->score;
	anchs->seq1end = i;
	anchs->seq2end = anchs->seq2start + osize;
	temp->seq1start = mybp+1;
	temp->seq2start = anchs->seq2end + 1;
	isfrst = 1;
	searchint=&(temp->seq1start);
      }
      else {
	temp->score = anchs->score * scale;
	anchs->score -= temp->score;
	anchs->seq1end = anchs->seq1start + osize;
	anchs->seq2end = i;
	temp->seq1start = anchs->seq1end + 1;
	temp->seq2start = mybp+1;
	isfrst = 1;
	searchint=&(temp->seq2start);
      }
      assert (anchs->seq1start <= anchs->seq1end);
      assert (anchs->seq2start <= anchs->seq2end);
      prev = anchs;
      anchs = temp;

      if (anchs->seq1start >= anchs->seq1end || anchs->seq2start >= anchs->seq2end){
	//	fprintf (stderr, "5. The region from %d to %d was thrown away.\n", stmybp, mybp);
	temp = anchs;
	prev->next = anchs->next;
	anchs = anchs->next;
	free (temp);
	if (!anchs) break;
	searchint = (!which)?&(anchs->seq1start):&(anchs->seq2start);
      }	

      //      fprintf (stderr, "Now, I am looking for %d, isfrst=%d (%d %d).\n", *searchint, isfrst, temp->seq1start, temp->seq1end);
      //      fprintf (stderr, "Currently, we are position %d in the sequence.\n", mybp);
      continue;
    }
    mybp++;
    if (mybp==*searchint){
      if (isfrst) {
	*searchint = i;
	searchint = (!which)?&(anchs->seq1end):&(anchs->seq2end);
	stmybp = mybp;
	isfrst = !isfrst;
	//	fprintf (stderr, "2) Now, I am looking for %d, isfrst=%d.\n", *searchint, isfrst);
	//	fprintf (stderr, "Currently, we are position %d in the sequence.\n", mybp);
      }
    }
    if (mybp==*searchint){
      if (!isfrst){
	*searchint = i;

	assert (anchs->seq1start <= anchs->seq1end);
	assert (anchs->seq2start <= anchs->seq2end);
	
	if (which == 0 && anchs->seq1end - anchs->seq1start < ANCHOR_LENGTH_CUTOFF ||
	    which == 1 && anchs->seq2end - anchs->seq2start < ANCHOR_LENGTH_CUTOFF){
	  //	  fprintf (stderr, "4. The region from %d to %d was thrown away.\n", stmybp, mybp);
	  temp = anchs;
	  prev->next = anchs->next;
	  anchs = anchs->next;
	  free (temp);
	}
	else {
	  //	  fprintf (stderr, "3. The region from %d to %d was saved.\n", stmybp, mybp);
	  prev = anchs;
	  anchs = anchs->next;
	}
	if (!anchs)
	  break;
	searchint = (!which)?&(anchs->seq1start):&(anchs->seq2start);

	isfrst = !isfrst;
	//	fprintf (stderr, "Now, I am looking for %d, isfrst=%d.\n", *searchint, isfrst);
	//	fprintf (stderr, "Currently, we are position %d in the sequence.\n", mybp);
      }
    }
  }

  //  fprintf (stderr, "By the end, I have reached mybp=%d, stmybp=%d.\n", mybp, stmybp);
  //  fprintf (stderr, "   number of positions in sequence to remap = %d\n", mylen);
  //  fprintf (stderr, "   number of positions in other sequence = %d\n", olen);
  
  temp = wlist;
  wlist = wlist->next;
  free (temp);

  for (temp = wlist; temp; temp = temp->next){
    // fprintf (stderr, "(%d %d)=(%d %d) %f\n", temp->seq1start, temp->seq1end, temp->seq2start, temp->seq2end, temp->score);
    assert (temp->seq1start <= temp->seq1end);
    assert (temp->seq2start <= temp->seq2end);
    assert (temp->seq1start >= 0);
    assert (temp->seq2start >= 0);
    assert (temp->seq1end >= 0);
    assert (temp->seq2end >= 0);
  }

  return wlist;
}


int hllIntersection(hll *h1, hll *h2) {
  int i, j;
  int r1, r2;

  if (!h1 || !h2) return 0;

  i=MAX2(h1->seq1start, h2->seq1start);
  j=MIN2(h1->seq1end, h2->seq1end);
    
  r1 = ((i<j) ? j-i : 0);

  i=MAX2(h1->seq2start, h2->seq2start);
  j=MIN2(h1->seq2end, h2->seq2end);
    
  r2 = ((i<j) ? j-i : 0);

  return (MIN2(r1, r2));
}

int hllUnion(hll *h1, hll *h2) {
  int i, j;
  int r1, r2;

  if (!h1 && !h2) return 0;
  if (!h1) return MAX2(h2->seq1end - h2->seq1start,
		       h2->seq2end - h2->seq2start);
  if (!h2) return MAX2(h1->seq1end - h1->seq1start,
		       h1->seq2end - h1->seq2start);

  i=MIN2(h1->seq1start, h2->seq1start);
  j=MAX2(h1->seq1end, h2->seq1end);
    
  r1 = ((i<j) ? j-i : 0);

  i=MIN2(h1->seq2start, h2->seq2start);
  j=MAX2(h1->seq2end, h2->seq2end);
    
  r2 = ((i<j) ? j-i : 0);
  
  return (MAX2(r1, r2));
}


hll* hllJoin(hll *h1, hll *h2, int score) {
  int i, j;
  hll *res = malloc (sizeof(hll));

  
  res->seq1start=MIN2(h1->seq1start, h2->seq1start);
  res->seq1end=MAX2(h1->seq1end, h2->seq1end);
    
  res->seq2start=MIN2(h1->seq2start, h2->seq2start);
  res->seq2end=MAX2(h1->seq2end, h2->seq2end);
  res->score = score;

  return res;
}


int minHLL(hll *h1, hll *h2){
  int i, j;

  i=MIN2(h1->seq1end, h2->seq1end);
  return (i==h2->seq1end);
}


float scoreMerge(hll* h1, hll *h2) {
  float i, u;
  i = hllIntersection(h1, h2);
  u = hllUnion(h1, h2);

  return (h1->score + h2->score)*(i/u);
}


void printSeqsNames(align *a) {
  int i;
  printf("( ");
  for (i=0; i<a->numseq; i++) {
    printf("%s ", a->seqs[i]->name);
  }
  printf(")\n");
}


void printMyHLL(hll *myres) {
  /* 
  while(myres) {

    printf("***: (%d %d)=(%d %d)\n", 
	   myres->seq1start, myres->seq1end,
	   myres->seq2start, myres->seq2end);    

    myres=myres->next;
  }
  */
}

hll* mergeHLLs(hll* anchs1, int wh1, hll* anchs2, int wh2) {
  int i, j, mscore;
  hll* res=0, *temp;
  if(wh1) swapHLL(anchs1);
  if(wh2) swapHLL(anchs2);
  /*
  printf("anchs1: \n");
  printMyHLL(anchs1);
  printf("anchs2: \n");
  printMyHLL(anchs2);
  */
  if (anchs1==anchs2) {
    //    fprintf(stderr, "mergeHLLs called on same hll!\n");
    return anchs1;
  }

  while((anchs1 && anchs2)) {
    //    printf("calling printMyHLL!\n");
    // printMyHLL(res);
    if (hllIntersection(anchs1, anchs2)) {
      mscore = scoreMerge(anchs1, anchs2);
      if (MAX3(anchs1->score, anchs2->score, mscore) == mscore) {
	temp = hllJoin(anchs1, anchs2, mscore);
	temp->next = res;
	res = temp;
      }
    }
    if (minHLL(anchs1, anchs2)) {
      temp = anchs2->next;
      anchs2->next = res;
      res = anchs2;
      anchs2 = temp;
    }
    else {
      temp = anchs1->next;
      anchs1->next = res;
      res = anchs1;
      anchs1 = temp;
    }
  }
  if (anchs1 && !anchs2)
    while (anchs1) {
      temp = anchs1->next;
      anchs1->next = res;
      res = anchs1;
      anchs1 = temp;
    }
  if (!anchs1 && anchs2)
    while (anchs2) {
      temp = anchs2->next;
      anchs2->next = res;
      res = anchs2;
      anchs2 = temp;
    }
  return res;
}

int printTextAlign(FILE* outfile, align* myalign) {
  int s1=0, s2=0, c, k, i;
  int nlets=0;
  int* inds = (int*) malloc (sizeof(int)* myalign->numseq);
  if (!outfile)
    outfile = stdout;

  for (i=0; i< myalign->numseq; i++) {
    inds[i] = 1;
  }

  //  fprintf(outfile, "ALIGNMENT LENGTH=%d\n\n", myalign->algnlen);

  for (c = 1; c < myalign->algnlen; c = c + 60) {

    for (i=0; i< myalign->numseq; i++) {

      for (k = c; (k < (c + 60)) && (k < myalign->algnlen); k++) {

	if (myalign->algn[k] & (1<<i))
	  fprintf(outfile, "%c", myalign->seqs[i]->lets[inds[i]++]);
	else 
	  fprintf(outfile,"-");
	
      }
      fprintf(outfile,"\n");

    }
    for (i=4; i < CNTS_LEN; i++) {
      for (k = c; (k < (c + 60)) && (k < myalign->algnlen); k++) {
	fprintf(outfile, "%d", myalign->cnts[i][k] % 10 );
      }
      fprintf(outfile,"\n");
    }

    /*
    fprintf(outfile,"\n"); 
    for (k=c;(k < (c + 60)) && (k < myalign->algnlen); k++) {
      fprintf(outfile, "%d", k/100);
    }
    fprintf(outfile,"\n"); 
    for (k=c;(k < (c + 60)) && (k < myalign->algnlen); k++) {
      fprintf(outfile, "%d", (k/10)%10);
    }
    fprintf(outfile,"\n"); 
    for (k=c;(k < (c + 60)) && (k < myalign->algnlen); k++) {
      fprintf(outfile, "%d", k%10);
    }
    fprintf(outfile,"\n"); 
    */

    fprintf(outfile,"\n\n");
  }


  fprintf(outfile,"\n");
  free(inds);
}

int printFASTAAlign(FILE* outfile, align* myalign) {
  int s1=0, s2=0, c, k, i;
  int nlets=0;
  int* inds = (int*) malloc (sizeof(int)* myalign->numseq);
  if (!outfile)
    outfile = stdout;

  for (i=0; i< myalign->numseq; i++) {
    inds[i] = 1;
  }

  for (i=0; i< myalign->numseq; i++) {
    fprintf(outfile, ">%s\n", myalign->seqs[i]->name);
    for (c = 1; c < myalign->algnlen; c = c + 60) {
      for (k = c; (k < (c + 60)) && (k < myalign->algnlen); k++) {
	if (myalign->algn[k] & (1<<i))
	  fprintf(outfile, "%c", myalign->seqs[i]->lets[inds[i]++]);
	else 
	  fprintf(outfile,"-");
      }
      fprintf(outfile,"\n");
    }
  }
  fprintf(outfile,"\n");

  free (inds);
}

int printXMFAAlign(FILE* outfile, align* myalign) {
  int s1=0, s2=0, c, k, i;
  int nlets=0;
  int* inds = (int*) malloc (sizeof(int)* myalign->numseq);
  if (!outfile)
    outfile = stdout;

  for (i=0; i< myalign->numseq; i++) {
    inds[i] = 1;
  }

  for (i=0; i< myalign->numseq; i++) {
    fprintf(outfile, ">%d:%d-%d + %s\n", myalign->seqs[i]->index, myalign->seqs[i]->leftbound,
	    myalign->seqs[i]->rightbound-1, myalign->seqs[i]->name);
    for (c = 1; c < myalign->algnlen; c = c + 60) {
      for (k = c; (k < (c + 60)) && (k < myalign->algnlen); k++) {
	if (myalign->algn[k] & (1<<i))
	  fprintf(outfile, "%c", myalign->seqs[i]->lets[inds[i]++]);
	else 
	  fprintf(outfile,"-");
      }
      fprintf(outfile,"\n");
    }
    fprintf(outfile,"\n");

  }

  free (inds);
}




void freeHLLs(hll *myHLL) {
  hll* a = myHLL;
  while (a) {
    myHLL = myHLL->next;
    free (a);
    a = myHLL;
  }
}


void freeSequence(seq *mySeq) {
  free(mySeq->rptr);
  free(mySeq->name);
  // rptr is a utility pointer, do not free
  // filename is not allocated, do not free
  free(mySeq);
}

void freeAlign(align *myAlign) {
  int i;
  //  if (freed[myAlign->num]) {
  //    fprintf (stderr, "Something very wrong... %d/%d", myAlign->num, freedsize);
  //  }
  assert (myAlign->dirty != 23);

  if (myAlign->nextalign) {
    myAlign->nextalign->dirty--;
    if (!myAlign->nextalign->dirty){
      freeAlign(myAlign->nextalign);
    }
  }
  myAlign->nextalign = 0;
  myAlign->dirty = 23;
  
  if (myAlign->algn){
    free(myAlign->algn);
    myAlign->algn = (long long int *) 0;
  }

  for (i=0; i<CNTS_LEN; i++) {
    if (myAlign->cnts[i]){
      free(myAlign->cnts[i]);
      myAlign->cnts[i] = (char *) 0;
    }
  }
  
  // sequences not freed
  // HLLs not freed
  if (freed)
    freed[myAlign->num] = 1;
  free(myAlign);
}

/*
void setScores(int gapstartV, int gapcontV, int gapendV, int gapperseqV, int overlapV, int glwidthV) {
  gapstart = gapstartV;
  gapcont = gapcontV;
  gapend = gapendV;
  gapperseq = gapperseqV;
  overlap = overlapV;
  glwidth = glwidthV;
  }*/



