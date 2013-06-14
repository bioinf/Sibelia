#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "diagmatrix.h"
#include "filebuffer.h"

#define NUC_FILE "nucmatrix.txt"
#define NUC_FILE_SIZE 6

#define MAX_SQ_SIZE (500 * (1 << 20))
#define BIG_SQ_WIDTH 20

#define VER_NUM "1.1"

#define INSERTION 2
#define DELETION 3

#define ISCB(c) ((c)=='.')

#define MIN2(x,y)   ( (x) >= (y) ? (y) : (x) )
#define MAX2(x,y)   ( (x) >= (y) ? (x) : (y) )
#define MAX3(x,y,z)  MAX2(MAX2(x,y),z)

#define WEQ2(x,y,a)  ((x==a)? 0: (y==a)? 1:-1)
#define WEQ3(x,y,z,a)  ((x==a)? 0: (y==a)? 1: (z==a)? 2:-1)

align* makeAlign(dmat* mydm, char* seq1, char* seq2);


char* alpha = "ATCGN.";

int s1start = 0;
int s1end = 0;
int s2start = 0;
int s2end = 0;
int gapstart = -1500;
int gapcont = -50;
//int match =12;
//int mismatch = -8;
int overlap = 0;
int glwidth= 15;
char dobin = 0;
char domfa = 0;
char doxmfa = 0;
FILE* ancfile = 0;
FILE* outfile;

int substmatrix[256][256];


seq* readfile(FILE* input, int seqnum) {
  char* res = (char*) malloc(sizeof(char)*2);
  int ressize = 2, numread=1;
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
  res[0] = 0;
  currchar = fgetc(input);
  while ((currchar != '>') && (currchar != EOF)) {
    if (!isspace(currchar)) {
      currchar = toupper(currchar);
      if (!strchr(alpha, currchar)) {
	fprintf(stderr, "WARNING %c converted to 'N'\n", currchar);
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
      res = &res[s1start-1];
      res[s1end-s1start+1] = 0;
      numread = s1end-s1start+1;
    }
    else {
      s1start = 1;
      s1end = numread;
    }
  }
  else {
    if (s2start > 0) {
      res = &res[s2start-1];
      res[s2end-s2start+1] = 0;
      numread = s2end-s2start+1;
    }
    else {
      s2start = 1;
      s2end = numread;
    }
  }
  myseq->lets = res;
  myseq->numlets = numread-1;
  //  printf("red %d lets\n",numread);
  return myseq;
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

  fscanf (file, "%d", &gapstart);
  fscanf (file, "%d", &gapcont);
  
  fclose (file);
}

void paramParse(int argc, char** argv) {
  int i = 3;
  for ( ; i < argc; i++) {
    if (!strcmp(argv[i], "-gs") || !strcmp(argv[i], "-GS")) {
      gapstart = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-gc") || !strcmp(argv[i], "-GC")) {
     gapcont = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-bin") || !strcmp(argv[i], "-BIN")) {
      dobin =1;
    }
    else if (!strcmp(argv[i], "-mfa") || !strcmp(argv[i], "-MFA")) {
      domfa =1;
    }
    else if (!strcmp(argv[i], "-xmfa") || !strcmp(argv[i], "-XMFA")) {
      doxmfa =1;
    }
    /*    else if (!strcmp(argv[i], "-mt") || !strcmp(argv[i], "-MT")) {
      match = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-ms") || !strcmp(argv[i], "-MS")) {
      mismatch = atoi(argv[++i]);
      }*/
    else if (!strcmp(argv[i], "-bw") || !strcmp(argv[i], "-BW")) {
      glwidth = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-s1") || !strcmp(argv[i], "-S1")) {
      s1start = atoi(argv[++i]);
      s1end = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-s2") || !strcmp(argv[i], "-S2")) {
      s2start = atoi(argv[++i]);
      s2end = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-anc") || !strcmp(argv[i], "-ANC")) {
      if (!(ancfile = fopen(argv[++i],"r"))) {
	printf("couldnt open anchors file %s\n",argv[i]);
	exit(2);
      }
    }
    else if (!strcmp(argv[i], "-out") || !strcmp(argv[i], "-OUT")) {
      if (!(outfile = fopen(argv[++i],"w"))) {
	printf("couldnt open output file %s\n",argv[i]);
	exit(2);
      }
    }
  }

  readSubstMatrix (NUC_FILE, NUC_FILE_SIZE);
}

void usage() {
  printf("usage: \norder seq1file seq2file [options]\n\n");
  printf("Options:\n");
  printf("-gs #  = Gap Start [default -100]\n");
  printf("-gc #  = Gap Continue [default -2]\n");
  /*  printf("-mt #  = MaTch [default 12]\n");
      printf("-ms #  = MiSmatch [default -8]\n");*/
  printf("-bw #  = Barrel Width around conserved regions [default 15]\n");
  printf("-anc anchorfile  = specify an anchorfile to use [default no file]\n");
  printf("-out outfile  = write output to outfile [default screen]\n");
  printf("-bin   = write output in BINary format [default text]\n");
  printf("-mfa   = write output in MultiFAsta format [default text]\n");
  printf("-s1 # # = use the given substring of the query [default whole]\n");
  printf("-s2 # # = use the givensubstring of the dbase [default whole]\n");
  printf("-version = prints the version of this ORDER\n");
}

hll* readAncFile(seq* seq1, seq* seq2) {
  hll *myres = 0, *tt;
  char buff[256];
  int i=0;
  
  while (!feof(ancfile)) {
    if (!fgets(buff, 256, ancfile)) {
      break;
    }
    tt = (hll*) malloc(sizeof(hll));
    sscanf(buff, "(%d %d)=(%d %d) %*f", &tt->seq1start, &tt->seq1end,
	   &tt->seq2start, &tt->seq2end);

    if ((tt->seq1start >= s1start && tt->seq1end <= s1end || s1start == 0 && s1end == 0) &&
	(tt->seq2start >= s2start && tt->seq2end <= s2end || s2start == 0 && s2end == 0)){
      
      if (tt->seq1start <= 0 && tt->seq1end <= 0) continue;
      if (tt->seq2start <= 0 && tt->seq2end <= 0) continue;
      if (tt->seq1start > s1start + seq1->numlets && tt->seq1end > s1start + seq1->numlets) continue;
      if (tt->seq2start > s2start + seq2->numlets && tt->seq2end > s2start + seq2->numlets) continue;

      if (s1start > 0){
	tt->seq1start = MAX2 (tt->seq1start - s1start + 1, 1);
	tt->seq1end = MIN2 (tt->seq1end - s1start + 1, s1end);
      }
      if (s2start > 0){
	tt->seq2start = MAX2 (tt->seq2start - s2start + 1, 1);
	tt->seq2end = MIN2 (tt->seq2end - s2start + 1, s2end);
      }
      
      tt->seq1start = MAX2 (tt->seq1start, 1);
      tt->seq2start = MAX2 (tt->seq2start, 1);
      tt->seq1end = MIN2 (tt->seq1end, seq1->numlets);
      tt->seq2end = MIN2 (tt->seq2end, seq2->numlets);

      tt->next = myres;
      i++;
      myres = tt;      



    }
  }
  fprintf(stderr,"read %d anchs\n", i);
  return myres;
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
  //  printf("dt = %d\n", dt);
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
    //    printf("dn =%d  ", *dn);
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
//    if (*dn < 0 || *dn >= 34939) fprintf (stderr, "%d %d\n", *dn, dt);
    starts[*dn] = MAX2(elem - width, 0);
    ends[*dn] = MIN2(elem+width, dlen-1);
    //    printf("BARREL %d  %d %d\n",*dn,starts[*dn],ends[*dn]);
  }
}

void mkSquare(int s1, int s2, int e1, int e2, int *dn, int dt, int* starts, int *ends, dmat* mydm) {
  int dists[2];
  long long int size = ((long long int)e1-(long long int)s1)
    * ((long long int)e2-(long long int)s2);
  //  printf("dt = %d\n", dt);
  //  printf("SQ: %d, %d to %d, %d\n", s1,s2,e1,e2);
  if (size > MAX_SQ_SIZE) {
    fprintf (stderr, "SQUARE TOO BIG: %d,%d to %d,%d\n", s1, e1,s2,e2);
    mkSquare(s1, s2, (s1+e1)/2+glwidth, (s2+e2)/2+glwidth, dn, (*dn+dt)/2, starts, ends, mydm);
    mkSquare((s1+e1)/2-glwidth, (s2+e2)/2-glwidth, e1, e2, dn, dt, starts, ends, mydm);
    return;
  }
  for ( ; *dn < dt; (*dn)++) {
    //    printf("square dn = %d\n", *dn);
    if (*dn < mydm->d2) {
      dists[0] = s1-1;
      dists[1] = *dn - e2;
    }
    else {
      dists[0] = mydm->d2 - e2;
      dists[1] = s1 - (*dn - mydm->d2)-1;
    }
//    if (*dn < 0 || *dn >= 34939) fprintf (stderr, "%d\n", *dn);
    starts[*dn] = MAX2(dists[0], dists[1]);

    if (*dn < mydm->d2) {
      dists[0] = e1-1;
      dists[1] = *dn - s2;
    }
    else {
      dists[0] = mydm->d2 - s2;
      dists[1] = e1 - (*dn-mydm->d2)-1;
    }
    ends[*dn] = MIN2(dists[0], dists[1]);
    //    printf("SQUARE %d  %d %d\n",*dn, starts[*dn],ends[*dn]);
  }
}

void doShapes(hll* myres, dmat* mydm, int* starts, int *ends) {
  int p1=MAX2(overlap,glwidth)+1, p2=MAX2(overlap,glwidth)+1; 
  int t1, t2;
  int dn = 1, dt;
  int width = glwidth;
  while (myres) {
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


void parseAncs(dmat* mydm, seq* seq1, seq* seq2) {
  int *starts = (int*) malloc(sizeof(int)*(seq1->numlets + seq2->numlets+2));
  int *ends = (int*) malloc(sizeof(int)*(seq1->numlets + seq2->numlets+2));
  hll* myres = 0;
  if (ancfile) {
    myres = readAncFile(seq1, seq2);
  }
  //  printf("khe0\n");
  doShapes(myres, mydm, starts, ends);
  //  printf("khe1\n");
  DMinitDiag(mydm, starts,ends);
  //  printf("khe2\n");
  free(starts);
  free(ends);
}

void doAlign(dmat* mydm, seq* seq1, seq* seq2) {
  align *a = (align*) makeAlign(mydm, seq1->lets, seq2->lets);
  //  printf("into printing\n");
  if (!dobin && !domfa && !doxmfa)
    printTextAlign(seq1->lets, seq2->lets, a);
  else if (!domfa && !doxmfa)
    printBinAlign(seq1->lets, seq2->lets, a);
  else if (!doxmfa)
    printMFAAlign(seq1->lets, seq2->lets, a, seq1->name, seq2->name);
  else 
    printXMFAAlign(seq1->lets, seq2->lets, a, seq1->name, seq2->name);
  //  printf("doneprinting\n");
}

int main(int argc, char** argv) {
  FileBuffer fseq1, fseq2;
  seq *seq1, *seq2;
  dmat* mydm;
  if (argc < 3) {
    if (argc == 2)
      if (!strcmp(argv[1], "-version") || !strcmp(argv[1], "-Version")) {
	printf("ORDER version %s\n", VER_NUM);
	exit(0);
      }
    usage();
    return 1;
  }
  if (!(fseq1 = FileOpen(argv[1]))) {
    printf("couldnt open query file %s\n",argv[1]);
    usage();
    return 2;
  }
  if (!(fseq2 = FileOpen(argv[2]))) {
    printf("couldnt open dbase file %s\n",argv[2]);
    usage();
    return 2;
  }
  outfile = stdout;
  paramParse(argc, argv);
  seq1 = FileRead(fseq1, s1start, s1end, VER_ORDER);
  seq2 = FileRead(fseq2, s2start, s2end, VER_ORDER);
  if (s1start == s1end && s1end == 0) {
    s1start = 1;
    s1end = seq1->numlets;
  }
  if (s2start == s2end && s2end == 0) {
    s2start = 1;
    s2end = seq2->numlets;
  }
  mydm = makeDM(seq1->numlets+1, seq2->numlets+1);
  parseAncs(mydm, seq1, seq2);
  doAlign(mydm, seq1, seq2);
  return 0;
}


inline int ismatch(char a, char b) {
  return a == b;
}

inline int matchscore (unsigned char a, unsigned char b) {
  return substmatrix[a][b];
  /*
    
  if (!a || !b)
    return 0;
  if (a == 'N' || b == 'N')
    return 0;
  if (a == b)
    return match;
  return mismatch;
  */
}

void reverse (char* a, int length) {
  char lft;
  int i;
  for (i=0; i < length/2; i++) {
    lft = a[i];
    a[i] = a[length-i-1];
    a[length-i-1] = lft;
  }
}

align* getChain(dmat* mydm, char* seq1, char* seq2, int x, int y, int inrun) {
  int temp;
  align *res = (align*) malloc (sizeof(align)), *help; 
  char* almt = (char*) malloc ( sizeof(char));
  int i=0, almtsize = 1, which;
  char zz;
  zz = DMgetPtr(mydm, x, y); 
  
  res->dirty = 0;
  res->nextalign = 0;
  res->algn = 0;
  res->algnlen = 0;

  do { 
    //    printf("I am at %d,%d %x\n", x,y, zz);
    which = zz & Mmask;

    if (which == 0x3) {
      help = DMgetNeck(mydm, x, y,inrun);
      if (!help) {
	return res;
      }
      help->dirty = 1;
      res->nextalign = help;
      break;
    }

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

    if (which == 0) {
      inrun = 0;
      almt[i++] = ismatch(seq1[x-1], seq2[y-1]);
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
     almt = (char *) realloc (almt, sizeof(char)* (almtsize *= 2));
   }

  }  while (x > 0 && y > 0);


  //  printf("gotChain\n");
  reverse(almt, i);
  res->algn = almt;
  res->algnlen = i;
  //  printf("done w it\n");
  return res;
}

void saveNeck(dmat* mydm, char* seq1, char* seq2, int neckdiag) {
  int size1, size2, x1, x2, y1, y2;
  alel *first = DMgetDiagStart(mydm, neckdiag-1, &size1, &x1, &y1),
    *second = DMgetDiagStart(mydm, neckdiag, &size2, &x2, &y2);
  int i, j;
  align* a;

  DMnextNecks(mydm, neckdiag);
  for (i=0; i<size2; i++,x2++,y2--) {
    for (j=0; j<3; j++) {
      a = getChain(mydm, seq1, seq2, x2, y2, j);
      DMsetNeck(mydm, a, x2, y2, j);
    }
  }
  for (i=0; i<size1; i++,x1++,y1--) {
    for (j=0; j<3; j++) {
      a = getChain(mydm, seq1, seq2, x1, y1, j);
      DMsetNeck(mydm, a, x1, y1, j);
    }
  }
}

void freeAlign(align* t) {
  free(t->algn);
  free(t);
}

void joinAligns (align* a) {
  align *n = a->nextalign, *t;
  char* temp, *temp2;
  int totsize=0;
  for (t = a; t; t = t->nextalign) {
    totsize += t->algnlen;
  }
  temp = (char *) malloc (totsize*sizeof(*temp));
  temp2 = temp + totsize;
  totsize = 0;
  for (t=a; t; t = t->nextalign) {
    totsize += t->algnlen;
    memcpy(temp2-totsize, t->algn, t->algnlen*sizeof(*temp));
  }
  free (a->algn);
  a->algn = temp;
  a->algnlen = totsize;
  for (a = a->nextalign; a;) {
    t = a;
    a = a->nextalign;
    freeAlign(t);
  }
}

align* makeAlign(dmat* mydm, char* seq1, char* seq2) {
  int i, j;
  int x, y, size;
  alel *curr, *pasts0, *pasts1, *pasts2;
  align* a;
  char isneck;
  int ndiags = mydm->d1 + mydm->d2 -1;
  register int s1, s2, s3;
  register char ptr;

  isneck = DMnextDiag(mydm);
  curr = DMgetDiagStart(mydm, 1, &size, &x, &y);
  curr->N = curr->O = gapstart;
  curr->M = 0;
  DMsetPtr(mydm, 0, 1, 1);
  //  printf("[%d %d]=%d %d %d\n",x,y,curr->M, curr->N, curr->O); 
  for (i = 2; i <= ndiags; i++) {
    isneck = DMnextDiag(mydm);
    if (!(i%10000))
      fprintf(stderr, "WORKING %d/%d\n", i/10000, ndiags/10000);
    curr = DMgetDiagStart(mydm, i, &size, &x, &y);

    pasts2 = DMgetElem(mydm, x-1, y);
    pasts1 = DMgetElem(mydm, x-1, y-1);
    for (j = 0; j < size; j++) {

      /***************************************************/
      pasts0 = pasts2;
      pasts2 = DMgetElem2(mydm, x, y-1, pasts2);

      s1 = pasts1->M;
      s2 = pasts1->N + ((ISCB(seq2[y-1]))?0:gapcont);
      s3 = pasts1->O + ((ISCB(seq1[x-1]))?0:gapcont);
      curr->M = matchscore (seq1[x-1], seq2[y-1]);
      if (s1 >= s2){
	if (s1 >= s3){ curr->M += s1; /*ptr = 0;*/ }
	else         { curr->M += s3; /*ptr = 2;*/ }
      }
      else {
	if (s2 >= s3){ curr->M += s2; /*ptr = 1;*/ }
	else         { curr->M += s3; /*ptr = 2;*/ }
      }

      s1 = curr->M + ((ISCB(seq2[y-1]))?0:gapstart);
      s2 = pasts0->N + ((ISCB(seq2[y-1]))?0:gapcont);
      if (s1 >= s2){ curr->N = s1; ptr = 0; }
      else         { curr->N = s2; ptr = 4; }
      
      s1 = curr->M + ((ISCB(seq1[x-1]))?0:gapstart);
      s2 = pasts2->O + ((ISCB(seq1[x-1]))?0:gapcont);
      if (s1 >= s2){ curr->O = s1; }
      else         { curr->O = s2; ptr |= 8; }
      
      s1 = curr->M;
      s2 = curr->N;
      s3 = curr->O;
      if (curr->M >= curr->N){
	if (curr->M < curr->O)
	  ptr |= 2;
      }
      else {
	if (curr->N >= curr->O)
	  ptr |= 1;
	else
	  ptr |= 2;
      }
      //ptr |= WEQ3(curr->M, curr->N, curr->O, MAX3(curr->M, curr->N, curr->O));
      //ptr = ptr | (WEQ2(curr->M+gapstart, pasts0->N+gapcont, curr->N) << 2); 
      //ptr = ptr | (WEQ2(curr->M+gapstart, pasts0->O+gapcont, curr->O) << 3);
      /***************************************************/
      /*
	curr->M = MAX3(pasts[1]->M, pasts[1]->N+gapcont, pasts[1]->O+gapcont); 
	curr->M += matchscore(seq1[x-1], seq2[y-1]);
	curr->N = MAX2(curr->M+gapstart, pasts[0]->N+gapcont); 
	curr->O = MAX2(curr->M+gapstart, pasts[2]->O+gapcont); 
	ptr = WEQ3(curr->M, curr->N, curr->O, MAX3(curr->M, curr->N, curr->O));
	ptr = ptr | (WEQ2(curr->M+gapstart, pasts[0]->N+gapcont, curr->N) << 2); 
	ptr = ptr | (WEQ2(curr->M+gapstart, pasts[0]->O+gapcont, curr->O) << 3);
      */

      DMsetPtr(mydm, ptr, x, y);
      curr++; x++; y--;

      pasts1 = DMgetElem2(mydm, x-1, y-1, pasts1);
    }
    if ((i < ndiags - 2) && isneck) {
      saveNeck(mydm, seq1, seq2, i);
    }
  }
  mydm->currneck++;
  a = getChain(mydm, seq1, seq2, mydm->d1, mydm->d2, 0);
  curr--;
  a->score = MAX3(curr->M, curr->N, curr->O);
  //  printf("here! %d\n", a);
  joinAligns(a);
  return a;
}

int printBinAlign(char* seq1, char* seq2, align* myalign) {
  int s1=1, s2=1, c;
  char lets[256];
  char left, right;
  //  fprintf(stderr,"kuku\n");
  for (c = 0; c < 256; c++)
    lets[c] = -1;
  lets['A'] = 1;  lets['C'] = 2;  lets['T'] = 3;  lets['G'] = 4; lets['N'] = 5; lets['.'] = 0;
  for (c = 1; c < myalign->algnlen; c++) {
    left=right=0;
    if (myalign->algn[c] != DELETION)
      left = lets[seq1[s1++]];
    if (myalign->algn[c] != INSERTION)
      right = lets[seq2[s2++]];
    right = right | (left << 4);
    putc(right, outfile);
  }
  fclose(outfile);
}

int printTextAlign(char* seq1, char* seq2, align* myalign) {
  int s1=1, s2=1, c, k;
  int nm=0, nga=0, ngb=0, nlets=0;
  int hasst=0;
  for (c = 1; c < myalign->algnlen; c = c + 60) {
    for (k = c; (k < (c + 60)) && (k < myalign->algnlen); k++) {
      if (myalign->algn[k] != DELETION)
	fprintf(outfile, "%c", seq1[s1++]);
      else {
	fprintf(outfile,"-");
	if (hasst)
	  nga++;
      }
    } 
    fprintf(outfile,"\n");
    for (k = c; (k < (c + 60)) && (k < myalign->algnlen); k++) {
      if (myalign->algn[k] == 1) {
	fprintf(outfile, ":");
	nm++; 
	nlets++;
	hasst = 1; 
      }
      else {
	fprintf(outfile, " ");
	if (hasst) nlets++;
      }
    } 
    fprintf(outfile, "\n");
    for (k = c; (k < (c + 60)) && (k < myalign->algnlen); k++) {
      if (myalign->algn[k] != INSERTION)
	fprintf(outfile, "%c", seq2[s2++]);
      else {
	fprintf(outfile, "-");
	if (hasst)
	  ngb++;
      }
    } 
    fprintf(outfile, "\n\n");
  }
  fprintf(outfile,"score = %d, nmatches = %d, nga=%d, ngb=%d nletters=%d, perc = %f\n",
	 myalign->score,nm,nga,ngb,nlets,(float)nm/(float)nlets);
  fprintf(outfile,"\n");
}

int printMFAAlign(char* seq1, char* seq2, align* myalign, char* n1, char* n2) {
  int s1=1, s2=1, c, k;
  int nm=0, nga=0, ngb=0, nlets=0;
  int hasst=0;
  fprintf(outfile,">%s\n", n1);
  for (c = 1; c < myalign->algnlen; c = c + 60) {
    for (k = c; (k < (c + 60)) && (k < myalign->algnlen); k++) {
      if (myalign->algn[k] != DELETION)
	fprintf(outfile, "%c", seq1[s1++]);
      else {
	fprintf(outfile,"-");
	if (hasst)
	  nga++;
      }
    } 
    fprintf(outfile,"\n");
  }
  fprintf(outfile,">%s\n", n2);
  for (c = 1; c < myalign->algnlen; c = c + 60) {
    for (k = c; (k < (c + 60)) && (k < myalign->algnlen); k++) {
      if (myalign->algn[k] != INSERTION)
	fprintf(outfile, "%c", seq2[s2++]);
      else {
	fprintf(outfile, "-");
	if (hasst)
	  ngb++;
      }
    } 
    fprintf(outfile, "\n");
  }
}

int printXMFAAlign(char* seq1, char* seq2, align* myalign, char* n1, char* n2) {
  int s1=1, s2=1, c, k;
  int nm=0, nga=0, ngb=0, nlets=0;
  int hasst=0;
  fprintf(outfile,">1:%d-%d + %s\n", s1start, s1end, n1);
  for (c = 1; c < myalign->algnlen; c = c + 60) {
    for (k = c; (k < (c + 60)) && (k < myalign->algnlen); k++) {
      if (myalign->algn[k] != DELETION)
	fprintf(outfile, "%c", seq1[s1++]);
      else {
	fprintf(outfile,"-");
	if (hasst)
	  nga++;
      }
    } 
    fprintf(outfile,"\n");
  }
  fprintf(outfile,">2:%d-%d + %s\n", s2start, s2end, n2);
  for (c = 1; c < myalign->algnlen; c = c + 60) {
    for (k = c; (k < (c + 60)) && (k < myalign->algnlen); k++) {
      if (myalign->algn[k] != INSERTION)
	fprintf(outfile, "%c", seq2[s2++]);
      else {
	fprintf(outfile, "-");
	if (hasst)
	  ngb++;
      }
    } 
    fprintf(outfile, "\n");
  }
}









