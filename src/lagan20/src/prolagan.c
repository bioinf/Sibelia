#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>

#include "skiplist.h"
#include "multial.h"
#include "filebuffer.h"

#define VER_NUM "1.1"
#define MIN2(x,y)   ( (x) >= (y) ? (y) : (x) )
#define MAX2(x,y)   ( (x) >= (y) ? (x) : (y) )

// Global variables

static int nested = 0;
static int postir = 0;
static int lazy = 0;
static int notree = 1;
static int verbose = 0;
static int numseqs = 0;
static int itertimes = 1;
static int cutoffmatch = 12;
static int translate = 0;
static int extend = 1;
static int fastreject = 0;
static int gapfreechunks = 0;

static align *simaligns[MAX_SEQ];
static char* lagan_dir;

static align *profile1 = 0;
static align *profile2 = 0;

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


void usage(void) {
  printf("mlagan seqfile_1 seqfile_2 [... seqfile_%d] [-parameters]\n\n",
	 MAX_SEQ);
  printf("-lazy : uses lazy mode\n");
  printf("-translate : use translated anchors\n");
  //  printf("-ext : extend the anchors\n");   This is now default
  printf("-fastreject : use fast rejection (tuned for human/mouse or closer)\n");
  //  printf("-gfc : find gap free chunks as anchors\n");   This is currently broken
  printf("-verbose : give debug output\n");
  printf("-tree \"(...)\" : runs with given phylogenetic tree\n");
  printf("-out \"filename\": outputs to filename\n");
  printf("-version : prints version info\n");
}

seq* readfile(FILE* input) {
  int seqstart=0;
  int seqend=0; 
  char* res = (char*) malloc(sizeof(char)*2);
  int ressize = 2, numread=1; //N at 1st letter
  char temp[256];
  seq* myseq = (seq*) malloc(sizeof(seq));
  char currchar;

  res[0] = 'N';
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
        fprintf(stderr, "Warning: %c converted to 'N'\n", currchar, alpha);
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

  if (seqstart > 0) {
    res = &res[seqstart-1];
    res[seqend-seqstart+1] = 0;
    numread = seqend-seqstart+1;
  }

  myseq->lets = res;
  myseq->numlets = numread;
  //  printf("read: %d lets\n",numread);
  return myseq;
}

int starts_with(char *str, char *word) {
  int len;
  char *first_word;

  len = strlen(str);
  first_word = (char *)malloc((len + 1) * sizeof(char));
  sscanf(str, "%s", first_word);
  return !strcmp(word, first_word);
}

align* findAlignByName(align *aligns[], char *name) {
  int i=0;
  // printf("findAlignByName: %s\n", name);
  while(i<numseqs) {
    if (starts_with(aligns[i]->seqs[0]->name, name)) {
      return(aligns[i]);
    }
    i++;
  }
  fprintf(stderr, "alignment not found for: %s", name);
  exit(2);
  return NULL;
}

int kk = 0;

// Profile stuff start

// replaces the sequence of same name with replacer, returning which was
// replaced or -1 if none.

int getSeqNumber(align* res, seq* replacer) {
  int i;
  for (i=0; i < res->numseq; i++) {
    if (!strcmp(res->seqs[i]->name, replacer->name)) {
      res->seqs[i] = replacer;
      return i;
    }
  }
  return -1;
}

void appendAlignProfile(align *res, seq* seqwgaps) {
  int i,j,k;
  res->seqs[res->numseq] = seqwgaps;
  for (i=1; i < res->algnlen; i++) {
    if (seqwgaps->lets[i] != '-') {
      k=strchr(alpha,seqwgaps->lets[i])-alpha;
      if (k < 4) {
	res->cnts[k][i]++;
      }
      res->algn[i] |= (1 << res->numseq);
      if (i > 0 && seqwgaps->lets[i-1] == '-')
	res->cnts[CNTS_GE][i]++;
    }
    else if (i > 0) {
      if (i > 0 && seqwgaps->lets[i-1] != '-') {
	res->cnts[CNTS_GS][i]++;      
      }
      else 
	res->cnts[CNTS_GC][i]++;      
      res->algn[i] |= (0 << res->numseq);
    }
  }
  res->numseq++;
}

align* readProfile(FileBuffer with_gaps) {
  int i,j;  
  seq* myseq;
  align* res = (align*) malloc (sizeof(align));
  res->score = 0;
  res->nextalign = 0;
  res->dirty = 0;
  res->numseq = 0;
  res->algnlen = -1;
  res->index = 32;
  
  while ( myseq = FileRead( with_gaps,0,0,VER_MLAGAN )) {
    //    fprintf(stdout, "seq: %s\n", myseq->lets);
    if (res->algnlen < 0) {
      res->algnlen = myseq->numlets;
      res->algn = (long long int*) malloc((res->algnlen+1) * sizeof(long long int));
      assert (res->algn);
      for (j=0; j<CNTS_LEN; j++) {
	res->cnts[j] = (char*) malloc((res->algnlen+1) * sizeof(char));    
	assert (res->cnts[j]);
      }
      for (i=0; i<= res->algnlen;i++) {
	for (j=0; j<CNTS_LEN; j++)
	  res->cnts[j][i] = 0; 
	res->algn[i] = 0;
      }
    }
    if ( res->algnlen != myseq->numlets) {
      fprintf (stderr, "Lengths screwed up!!!\n");
      exit(1);
    }
    appendAlignProfile(res, myseq);
  }
  if (verbose) {
    fprintf(stdout, "LOADED RES\n");
    printTextAlign(stdout,res);
  }
  return res;
}


// Profile stuff end


void printHLL(hll *myres) {
  fprintf(stderr, "into %d\n", ++kk);
  fflush(stderr);
  while(myres) {

    fprintf(stderr, "(%d %d)=(%d %d) %f\n", 
	   myres->seq1start, myres->seq1end,
	   myres->seq2start, myres->seq2end, myres->score);    
    fflush(stderr);
    myres=myres->next;
  }
}

hll* getAnchsFromFile(char *fname, FileBuffer f1, FileBuffer f2) {
  FILE *ancfile;
  hll *myres = 0, *tt = 0, *first = 0;
  char buff[256];
  int i=0, j=0;

  //  printf("getHLLFromNames: %s, %s\n", name1, name2);

  sprintf(buff, "%s.anchors", fname);
  ancfile=fopen(buff, "r");
  if(ancfile==NULL) {
    fprintf(stderr, "anchor file not found:: %s.anchors\n",
	   fname);
    exit(2);
  }

  while (!feof(ancfile)) {
    if (!fgets(buff, 256, ancfile)) {
      break;
    }
    tt = (hll*) malloc(sizeof(hll));
    sscanf(buff, "(%d %d)=(%d %d) %f", &tt->seq1start, &tt->seq1end,
           &tt->seq2start, &tt->seq2end, &tt->score);
    tt->next = myres;
    i++;
    myres = tt;
  }
  if (fastreject) {
    f1->startpos = MAX2(f1->startpos, myres->seq1end);
    f2->startpos = MAX2(f2->startpos, myres->seq2end);
    for (tt = myres; tt->next->next; tt = tt->next) {
      j++;
    }
    f1->endpos = MIN2(f1->endpos, tt->next->seq1start);
    f2->endpos = MIN2(f2->endpos, tt->next->seq2start);
    //    fprintf (stderr, "%d %d %d %d %d\n", j, f1->startpos, f1->endpos, f2->startpos, f2->endpos);
    myres = myres->next;
    tt->next = 0;
  }
  fprintf(stderr,"read %d anchs\n", i);
  fclose(ancfile);
  return myres;
}



hll* generateAnchors( FileBuffer a1, FileBuffer a2) {
  char buff[256];
  char fname[80];
  char *name1, *name2;
  char *endpnt;
  int diff1, diff2;
  align* temp;
  hll* res;
  char flip = 0;
  int retstat;

  name1 = strrchr (a1->filename, '/');
  if (!name1) name1 = a1->filename;
  else name1++;
  name2 = strrchr (a2->filename, '/');
  if (!name2) name2 = a2->filename;
  else name2++;

  endpnt = strchr ( name1, '.');
  diff1 = (endpnt)? endpnt - name1: strlen(name1);
  endpnt = strchr ( name2, '.');
  diff2 = (endpnt)? endpnt - name2: strlen(name2);
  strncpy (fname, name1, diff1);
  strncpy (fname+diff1, name2, diff2);
  fname[diff1+diff2] = 0;

  sprintf(buff, "%s/rechaos.pl %s %s -out %s.anchors %s %s %s %s %s\n",
          lagan_dir,
	  a1->filename,
	  a2->filename,
	  fname,
	  (extend ? "-ext" : ""),
	  (translate ? "-translate" : ""),
	  (fastreject ? "-fastreject" : ""),
	  (gapfreechunks ? "-gfc" : ""),
	  (lazy ? "-lazy" : ""));

  retstat = system(buff) >> 8;
  if (fastreject && (retstat == 3)) {
    return 0;
  }
  else if (retstat) {
    fprintf (stderr, "Error from rechaos\n");
    exit (1);
  }
  res = getAnchsFromFile(fname, a1, a2);
  return res;
}


void printFASTASeq(FILE *outfile, seq *myseq) {
  int i;
  //  printf("kva\n");
  if (!outfile)
    outfile = stdout;

  fprintf(outfile, ">%s\n", myseq->name);
  //  printf("kva2\n");
  for(i=0; i<myseq->numlets; i++)
    fprintf(outfile, "%c", myseq->rptr[i]);
  //  printf("kva %d\n",i);
  fprintf(outfile, "\n");
  
  if (outfile!=stdout) fclose(outfile);
}


hll* findBestChain(hptr* array, int arrsize) {
  sklst* skipper = makeSkLst();
  sle* help;
  int i;
  hll* t;
  for (i = 0; i < arrsize; i++) {
    if (array[i].isstart) {
      help = SLfind(skipper, array[i].myhll->seq2start);
      if (help->myelem) {
	array[i].myhll->bkptr = help->myelem;
	array[i].myhll->scoreSoFar = ((hll*)help->myelem)->scoreSoFar + array[i].myhll->score;
      }
      else {
	array[i].myhll->bkptr = 0;
	array[i].myhll->scoreSoFar = array[i].myhll->score;
      }
    }
    else {
      help = SLfind(skipper, array[i].myhll->seq2end);
      if (help->myelem && (array[i].myhll->scoreSoFar <= ((hll*)help->myelem)->scoreSoFar))
	continue;
      SLinsertAfter(skipper, help, array[i].myhll->seq2end, array[i].myhll);
      help = help->next[0];
      while (help->next[0] && 
	     ((hll*)help->myelem)->scoreSoFar >= ((hll*)help->next[0]->myelem)->scoreSoFar)
	SLremove(skipper, help->next[0]);
    }
  }
  t= (hll*)SLgetLast(skipper)->myelem;
  delSkLst(skipper);
  return t;
}


hll* remakeHLL(hll* bestPtr) { 
  int len;
  hll *res=0;
  hll *temp, *t2, *t3;
  int i, bestscore=-1;
  for (temp = bestPtr; temp; temp = temp->bkptr) {
    temp->next=res;
    temp->dirty = 1;
    res=temp;    
  }
  
  return res;
}


hll* reanchorHLL(hll* mylist) {

  hll *temp, *best, *t2;
  int numhits=0, i=0;
  hptr* myptrs;

  temp=mylist;
  while (temp) { numhits++; temp->dirty = 1; temp=temp->next; }

  myptrs = (hptr*) malloc (sizeof(hptr) * numhits *2);
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
  temp=best;
  while (temp) { temp->dirty = 0; temp=temp->bkptr; }
  temp=mylist;
  while (temp) { t2 = temp; temp=temp->next; if (t2->dirty) free(t2); }

  best = remakeHLL(best);
  //  printf("newbest\n");
  //  printHLL(best);
  free (myptrs);
  return best;
}


void orderAligns(align *a1, align *a2,
		 align **first, align **second,
		 int *index, int *hllindex) {
  int a1index, a2index;

  a1index = a1->index; 
  a2index = a2->index;
  
  if (a1index > a2index) {    
    *first = a2;
    *second = a1;
    *index = a2index;
    *hllindex = a1index;
  } else {
    *first = a1;
    *second = a2;
    *index = a1index;
    *hllindex = a2index;
  }
}


void doRemapHLLs(align *aligns[], align *uni, int *index, int hllindex) {
  int i, mapi, done=0;

  // take all hlls into first, and into the second and remap them

  for(mapi=*index; !done; mapi=hllindex)  {

    for (i=0; i<mapi; i++) {
      if (aligns[i]->hlls[mapi] != NULL && i != *index) {
	// remap them into i
	//	fprintf(stderr, "\n called1 %d %d(%d)\n", i, mapi, *index);
	aligns[i]->hlls[mapi] = remapHLLs(aligns[i]->hlls[mapi],
					  1, uni, 
					  (mapi!=*index));
      }
    }
    for (i=mapi+1; i<numseqs; i++) {
      if (aligns[mapi]->hlls[i] != NULL && i != hllindex) {
	// remap them into first or second
	//	fprintf(stderr, "\n called2 %d %d(%d)\n", mapi, i,*index);
	aligns[mapi]->hlls[i] = remapHLLs(aligns[mapi]->hlls[i],
					  0, uni,
					  (mapi!=*index));
      }
    }
    if (mapi==hllindex) done=1;
  }

  // free memory?  what's that?
  //  aligns[*index] = result;
  //  aligns[hllindex] = result;


}

void doReanchorHLLs(align *aligns[],
		 int *index, int hllindex) {
  int i;

  // for each pair of hlls from (i to first) and (i to second)

  for(i=0; i<*index; i++) {
    aligns[i]->hlls[*index] = 
      reanchorHLL(mergeHLLs(aligns[i]->hlls[*index], 0, 
			    aligns[i]->hlls[hllindex], 0));

    //    if (verbose) {
    //  printf("aligns[%d]->hlls[%d]\n",i ,*index);
    //    printHLL(aligns[i]->hlls[*index]);
    //   }
    aligns[i]->hlls[hllindex] = 0;
  }
  for(i=*index+1; i<hllindex; i++) {
    aligns[*index]->hlls[i] = 
      reanchorHLL(mergeHLLs(aligns[*index]->hlls[i], 0, 
			    aligns[i]->hlls[hllindex], 1));
    //  if (verbose) {
    //  printf("aligns[%d]->hlls[%d]\n",*index ,i);
    //    printHLL(aligns[*index]->hlls[i]);
    //  }
    aligns[i]->hlls[hllindex] = 0;
  }
  for(i=hllindex+1; i<numseqs; i++) {
    aligns[*index]->hlls[i] =  
      reanchorHLL(mergeHLLs(aligns[*index]->hlls[i], 0, 
			    aligns[hllindex]->hlls[i], 0));
    // if (verbose) {
    //  printf("aligns[%d]->hlls[%d]\n", *index, i);
    //    printHLL(aligns[*index]->hlls[i]);
    // }
    aligns[hllindex]->hlls[i] = 0;
  }
}


align* processAnchors(align *aligns[], align *a1, align *a2, int *index) {
  int hllindex;
  align *first, *second, *result, *uni;

  result = (align*) malloc(sizeof(align));
  
  assert (result);
  result->score = -1;
  result->numseq = a1->numseq + a2->numseq;
  result->algnlen = -1;
  result->nextalign = 0;
  result->dirty = 0;

  orderAligns(a1, a2, &first, &second, index, &hllindex);

  if (verbose)
    printHLL(aligns[first->index]->hlls[hllindex]);  

  //  result = makeAlign(first, second, aligns[first->index]->hlls[hllindex], &uni);
  result->index = *index;

  doReanchorHLLs(aligns, index, hllindex);

  fprintf(stderr,"done reanchor, leaving processAnchors\n");
  return(result);
}

align* processAlign(align *aligns[], align *a1, align *a2, int *index) {
  int hllindex;
  align *first, *second, *result, *uni;

  fprintf(stderr, "into processalign\n");

  orderAligns(a1, a2, &first, &second, index, &hllindex);

  if (verbose)
    printHLL(aligns[first->index]->hlls[hllindex]);  
  
  fprintf(stderr, "about to make\n");
  result = makeAlign(first, second, aligns[first->index]->hlls[hllindex], &uni);
  fprintf(stderr, "done make\n");
  result->index = *index;
  return(result);
}


align* iterativeImprovement (align *current, align *rpntree[], int length) {
  int converged = 0;
  int i=0, oldscore, cutoff;
  seq *removed;
  align *readd, *old, *new;
  hll* anchs, *tt;
  if (current->numseq <= 2)
    return current;
  //  printf("iterative improvement!\n");

  cutoff = cutoffmatch * 100;
  fprintf(stderr, "cutoff = %d\n", cutoff);
  while (!converged) {

    // Throw out a sequence.  Calling code in multial.
    removed = current->seqs[0];
    new = findAlignByName(simaligns, removed->name);
    old = current;
    anchs = getAnchsFromAlign(current, 0, cutoff);
    current = removeSeq(current, 0);
    free (old);

    // Re-align this thrown-out sequence to the remaining alignment.

    current = makeAlign (current, new, anchs, &old);
    if (verbose) {
      printf("improved:\n");
      printHLL(anchs);  
      printTextAlign(stdout, current);  
    }
    while (anchs) {
      tt = anchs;
      anchs = anchs->next;
      free (tt);
    }
    free (old);

    i++;
    if (i==numseqs*itertimes) converged = 1;
  }
  return current;
}



int treeToRPN(char *treestr, align *stack[MAX_SEQ*2], int *depth) {

  int i=0; int j, k; 
  char buffer[256];

  while (treestr[i]!='(') { i++; } i++;

  while ((treestr[i] != ')') && (treestr[i] != '\0')) { 
    //    printf("%d: %s\n", *depth, treestr+i);

  
    if (treestr[i]=='(') {
      i += treeToRPN(treestr+i, stack, depth);
    }  
    else if (isalnum(treestr[i])) {
      k = 0;
      // push alignment
      while((!isspace(treestr[i])) && (treestr[i]!='(') && (treestr[i]!=')')) { 
	buffer[k++] = treestr[i++];
      }
      buffer[k] = 0;
      stack[(*depth)++]=findAlignByName(simaligns, buffer);
      //      printf("pushed: %s\n", stack[*depth-1]->seqs[0]->name);
    }
    else if (treestr[i]==')')
      // (*depth)++;
      break;
    else { i++; }

  }

  if (treestr[i]==')') {
    (*depth)++; //null is '+'
    return i+1;
  }
 if (treestr[i] == '\0') { 
   fprintf(stderr, "ERROR parsing tree, depth %d, %d chars read", *depth, i);
   exit(1);
 }
}

align* procStack(align* rpntree[MAX_SEQ*2], int length, align *myaligns[]) {
  align* stack[MAX_SEQ];
  int i = 0, sp = 0;
  int index=0;

  while (i < (length-1)) {

    if (rpntree[i]) {
      stack[sp++] = rpntree[i];
    }
    else {
      stack[sp-2] = processAnchors(myaligns, stack[sp-2], stack[sp-1], &index);
      stack[--sp] = 0;      
      //      if(verbose) printTextAlign(stdout, stack[sp-1]);  
    }
    i++;
  }
  if (rpntree[i]) {
    fprintf(stderr,"Unexpeceted error\n");
  }
  else {
    stack[sp-2] = processAlign(myaligns, profile1, profile2, &index);
    stack[--sp] = 0;      
    if(verbose) printTextAlign(stdout, stack[sp-1]);  
  }

  return stack[sp-1];
}


void graphCollapsal (align *simaligns[]) {
  
  // for now...
  
  fprintf(stderr, "Please specify a phylogenetic tree, using [-tree]\n");
  exit(1);
}

int parseParameters(int argc, char** argv, FileBuffer *files, char **treestr) {

  int i=1;

  FileBuffer fb;

  if (argc < 3) {
    if (argc == 2)
      if (!strcmp(argv[1], "-version") || !strcmp(argv[1], "-Version")) {
        fprintf(stderr, "PROLAGAN version %s\n", VER_NUM);
        exit(0);
      }
    usage();
    return 1;
  }
  while((argv[i][0]!='-')) {

    // Read in sequence files
   
    //    printf("sequence %d: %s\n", i, argv[i]);

    if (!(files[numseqs++] = FileOpen(argv[i]))) {
      fprintf(stderr, "couldnt open dbase file %s\n",argv[i]);
      usage();
      return 2;
    }

    //    seqs[numseqs] = FileRead(seqfile, 0, 0, VER_MLAGAN);
    //    seqs[numseqs]->filename = argv[i];    
    //    numseqs++;


    if(++i>=argc) break;
  }

  //  printf("\n");

  while (i<argc) {
   
    // printf("parameters: %s\n", argv[i]);

    if (!(strcmp(argv[i], "-nested") || 
	  strcmp(argv[i], "-nopost") || 
	  strcmp(argv[i], "-postir") || 
	  strcmp(argv[i], "-fastreject") || 
	  strcmp(argv[i], "-gfc") || 
	  strcmp(argv[i], "-lazy") || 
	  strcmp(argv[i], "-verbose") || 
	  strcmp(argv[i], "-out") ||
	  strcmp(argv[i], "-translate") ||
	  strcmp(argv[i], "-ext") ||
	  strcmp(argv[i], "-match") || strcmp(argv[i], "-mismatch") ||
	  strcmp(argv[i], "-pro1") || strcmp(argv[i], "-pro2") ||
	  strcmp(argv[i], "-gapstart") || strcmp(argv[i], "-gapend") ||
	  strcmp(argv[i], "-gapcont") || strcmp(argv[i], "-gapperseq") ||
	  strcmp(argv[i], "-overlap") || strcmp(argv[i], "-glwidth") ||
	  strcmp(argv[i], "-tree"))) {
      fprintf(stderr, "unrecognized parameter: %s\n", argv[i]);
      usage();
      return 1;
    }
    if (!strcmp(argv[i], "-nested")) { 
      nested = 1; 
    }

    if (!strcmp(argv[i], "-translate")) { 
      translate = 1; 
    }

    if (!strcmp(argv[i], "-ext")) {  //default, do not use
      extend = 1; 
    }


    if (!strcmp(argv[i], "-verbose")) { 
      verbose = 1; 
    }

    if (!strcmp(argv[i], "-postir")) { 
      postir = 1; 
    }
    if (!strcmp(argv[i], "-lazy")) { 
      lazy = 1; 
    }
    if (!strcmp(argv[i], "-fastreject")) { 
      fastreject = 1; 
    }
    if (!strcmp(argv[i], "-gfc")) {  //Broken, do not use
      gapfreechunks = 1; 
    }

    if (!strcmp(argv[i], "-out")) {
      i++;
      if ((i>=argc) || (argv[i][0]=='-')) {
	fprintf(stderr, "missing parameter specification for [-out].\n");
	return 1;
      }
      fprintf(stderr, "outputting to: %s\n", argv[i]);
      outfile = fopen(argv[i], "w");
      if (outfile==NULL) {
	fprintf(stderr, "error with output file...\n");
	exit(2);
      }
    }

    if (!strcmp(argv[i], "-tree")) {
      i++;
      if ((i>=argc) || (argv[i][0]=='-')) {
	fprintf(stderr, "missing parameter specification for [-tree].\n");
	return 1;
      }
      notree = 0;
      *treestr = argv[i];
      fprintf(stderr, "using given phylogenetic tree:\n%s\n", *treestr); 
    }

    if (!strcmp(argv[i], "-gapperseq")) {
      i++;
      if (i>=argc) {
	fprintf(stderr, "missing parameter specification for [-gapperseq].\n");
	return 1;
      }
      gapperseq = atoi(argv[i]);
      fprintf(stderr, "using gapperseq score: %d\n", gapperseq); 
    }
    if (!strcmp(argv[i], "-overlap")) {
      i++;
      if (i>=argc) {
	fprintf(stderr, "missing parameter specification for [-overlap].\n");
	return 1;
      }
      overlap = atoi(argv[i]);
      fprintf(stderr, "using overlap value: %d\n", overlap); 
    }
    if (!strcmp(argv[i], "-glwidth")) {
      i++;
      if (i>=argc) {
	fprintf(stderr, "missing parameter specification for [-glwidth].\n");
	return 1;
      }
      glwidth = atoi(argv[i]);
      fprintf(stderr, "using glwidth value: %d\n", glwidth); 
    }

    if (!strcmp(argv[i], "-pro1")) {
      i++;
      if (i>=argc) {
	fprintf(stderr, "missing filename for [-pro1].\n");
	return 1;
      }
      fb = FileOpen (argv[i]);
      profile1 = readProfile(fb);
      fprintf(stderr, "Profile1 is: %s\n", argv[i]); 
    }

    if (!strcmp(argv[i], "-pro2")) {
      i++;
      if (i>=argc) {
	fprintf(stderr, "missing filename for [-pro2].\n");
	return 1;
      }
      fb = FileOpen (argv[i]);
      profile2 = readProfile(fb);
      fprintf(stderr, "Profile2 is: %s\n", argv[i]); 
    }

    i++;
  }

  //  setScores(gapstart, gapcont, gapend, gapperseq, overlap, glwidth);

  return 0;
}

hll* updateAnchorPos(hll* myhll, FileBuffer f1, FileBuffer f2) {
  hll *res, *temp, *prev=0;
  res = myhll;
  fprintf (stderr, "Updating anchs...\n");
  for ( ; myhll; myhll = myhll->next) {
    myhll->seq1start -= (f1->startpos-1);
    myhll->seq1end -= (f1->startpos-1);
    myhll->seq2start -= (f2->startpos-1);
    myhll->seq2end -= (f2->startpos-1);
  }
  while (res && (res->seq1start < 0 || res->seq2start < 0)) {
    //    fprintf (stderr, "first..\n");
    temp = res;
    //    fprintf(stderr, "Tossed %d %d(%d %d)\n", temp->seq1end, temp->seq2end,
    //    	    f1->endpos, f2->endpos);    
    res = res->next;
    free(temp);
  }
  temp = res;
  while (temp && temp->seq1end < (f1->endpos-f1->startpos) && temp->seq2end < (f2->endpos-f2->startpos)) {
    //    fprintf (stderr, "second...\n");
    //       fprintf(stderr, "Kept %d %d(%d %d)\n", temp->seq1end, temp->seq2end,
    //       	    f1->endpos-f1->startpos, f2->endpos-f2->startpos);
    prev = temp;
    temp = temp->next;
  }
  if (prev) {
    temp = prev;
    prev = prev->next;
    temp->next = 0;
  }
  else if (temp == res) {
    res = 0;
  }
  else {
    //    fprintf (stderr, "returning %d\n", res);
    return res;
  }
  while ( prev ) {
    //    fprintf (stderr, "third...\n");
    //        fprintf(stderr, "Tossed %d %d(%d %d)\n", temp->seq1end, temp->seq2end,
    //        	    f1->endpos, f2->endpos);
    temp = prev; 
    prev = prev->next;
    free(temp);
  }
  return res;
}

int connectedGraph(hll* graph[MAX_SEQ][MAX_SEQ], int numseqs) {
  int M[MAX_SEQ][MAX_SEQ];
  int i, j, k;

  for (i = 0; i < numseqs - 1; i++){
    for (j = i + 1; j < numseqs; j++){
      M[i][j] = M[j][i] = (graph[i][j] != NULL);
    }
  }

  for (k = 0; k < numseqs; k++)
    for (i = 0; i < numseqs; i++)
      for (j = 0; j < numseqs; j++)
	if (M[i][k] && M[k][j]) M[i][j] = 1;

  k = 1;
  for (i = 0; k && i < numseqs; i++)
    k = M[0][i];

  return k;
}


int main(int argc, char** argv) {
  FileBuffer seqfile;
  seq **seqs;
  int i = 1, j = 1, x, y;
  int pro1cnt=0, pro2cnt=0;
  int pro1lst[MAX_SEQ], pro2lst[MAX_SEQ];
  int pro1ptr[MAX_SEQ], pro2ptr[MAX_SEQ];
  char command[256];

  char *treestr = NULL;
  align *stack[MAX_SEQ*2];
  align *final;
  align *myaligns[MAX_SEQ];
  hll* table[MAX_SEQ][MAX_SEQ];
  FileBuffer files[MAX_SEQ];

  outfile = stdout;
  lagan_dir = getenv ("LAGAN_DIR");
  if (!lagan_dir) {
    fprintf(stderr, "Environment variable LAGAN_DIR not set\n");
    exit(1);
  }

  buildcache();
  initLib();

  seqs = (seq**) malloc((argc-1)*sizeof(seq*));


  if (parseParameters(argc, argv, files, &treestr)) return 1;

  gapstart += gapcont;


  // Take all sequences and make simple alignments

  for (i=0; i<numseqs; i++) {
    seqs[i] = FileRead(files[i], 0, 0, VER_MLAGAN);
    seqs[i]->index = i+1;
    myaligns[i]=simaligns[i]=mkSimAlign(seqs[i]);
    simaligns[i]->index = i;
    x = getSeqNumber(profile1, seqs[i]);
    y = getSeqNumber(profile2, seqs[i]);
    if (x < 0 && y < 0) {
      fprintf(stderr, "Sequence %s not found in either profile!!!\n", seqs[i]->name);
      exit(1);
    }
    if (x >= 0 && y >= 0) {
      fprintf(stderr, "Sequence %s found in both profiles!!!\n", seqs[i]->name);
      exit(1);
    }
    if (x >= 0) {
      fprintf(stderr, "Sequence %s[%d/%d] in 1st profile\n", seqs[i]->name, i, numseqs);
      if (profile1->index > i) {
	profile1->index = i;
      }
      pro1lst[pro1cnt++] = i;
      pro1ptr[i] = x;
      pro2ptr[i] = -1;
    }
    if (y >= 0) {
      fprintf(stderr, "Sequence %s[%d/%d] in 2nd profile\n", seqs[i]->name, i, numseqs);
      if (profile2->index > i) {
	profile2->index = i;
      }
      pro2lst[pro2cnt++] = i;
      pro1ptr[i] = -1;
      pro2ptr[i] = y;
    }
  } 


  // Find all pairwise anchors.
  fprintf(stderr,"pro1cnt = %d, pro2cnt = %d\n", pro1cnt, pro2cnt);
  for (i=0; i<(numseqs-1); i++) {
    for (j=i+1; j<numseqs; j++) {
      simaligns[i]->hlls[j]=0;
    }
  }
  for (i=0; i< pro1cnt; i++) {
    for (j=0; j< pro2cnt; j++) {
      if (pro1lst[i] < pro2lst[j]) {
	simaligns[pro1lst[i]]->hlls[pro2lst[j]] = generateAnchors(files[pro1lst[i]], files[pro2lst[j]]);
	simaligns[pro1lst[i]]->hlls[pro2lst[j]] = remapHLLs(simaligns[pro1lst[i]]->hlls[pro2lst[j]],
							    0, profile1, pro1ptr[pro1lst[i]]);
	simaligns[pro1lst[i]]->hlls[pro2lst[j]] = remapHLLs(simaligns[pro1lst[i]]->hlls[pro2lst[j]],
							    1, profile2, pro2ptr[pro2lst[j]]);
      }
      else {
	simaligns[pro2lst[j]]->hlls[pro1lst[i]] = generateAnchors(files[pro2lst[j]], files[pro1lst[i]]);
	simaligns[pro2lst[j]]->hlls[pro1lst[i]] = remapHLLs(simaligns[pro2lst[j]]->hlls[pro1lst[i]],
							    0, profile2, pro2ptr[pro2lst[j]]);
	simaligns[pro2lst[j]]->hlls[pro1lst[i]] = remapHLLs(simaligns[pro2lst[j]]->hlls[pro1lst[i]],
							    1, profile1, pro1ptr[pro1lst[j]]);
      }
    }
  }

  //  printf("\n");

  for (i=0; i<MAX_SEQ*2; i++) {
    stack[i] = NULL;
  }


  /*
  for (i=0; i<(numseqs-1); i++) {
    for (j=i+1; j<numseqs; j++) {
      printf("Sanity Check: simaligns[%d]->hlls[%d].score=%g\n",
	     i,j,
	     simaligns[i]->hlls[j]==NULL ? 0 : simaligns[i]->hlls[j]->score);
    }
  }
  */

  // Processall closest pairs 

  if (notree) { // Not yet implemented
    graphCollapsal(myaligns);
  }
  else {

    fprintf(stderr, "\n****************************\n");
    fprintf(stderr, "gs: %d; ge: %d;\n", gapstart, gapend);
    fprintf(stderr, "gc: %d; gp: %d\n", gapcont, gapperseq);
    //fprintf(stderr, "match: %d; mismatch: %d\n", match, mismatch);
    fprintf(stderr, "overlap: %d; glwidth: %d\n", overlap, glwidth);
    fprintf(stderr, "\n****************************\n");

    i = 0;
    treeToRPN(treestr, stack, &i);
    final = procStack(stack, i, myaligns);
  }


  // Ouput end result.
  fprintf(stderr, "final alignment... \n");
  if (fastreject) {
    printXMFAAlign(outfile, final);
  }
  else {
    printFASTAAlign(outfile, final);
  }
  if (outfile != stdout) fclose (outfile);


  fprintf(stderr, "mlagan -- end.\n");
  return 0;
}














