#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>

#include "skiplist.h"
#include "multial.h"
#include "filebuffer.h"

#define VER_NUM "2.0"
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
  printf("-nested : runs improvement in a nested fashion\n");  
  printf("-postir : incorporates the final improvement phase\n");
  printf("-lazy : uses lazy mode\n");
  printf("-translate : use translated anchors\n");
  //  printf("-ext : extend the anchors\n");   This is now default
  printf("-fastreject : use fast rejection (tuned for human/mouse or closer)\n");
  //  printf("-gfc : find gap free chunks as anchors\n");   This is currently broken
  printf("-verbose : give debug output\n");
  printf("-tree \"(...)\" : runs with given phylogenetic tree\n");
  printf("-out \"filename\": outputs to filename\n");
  printf("-nucmatrixfile \"filename\": uses given substitution matrix instead of $LAGAN_DIR/nucmatrix.txt\n");
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
  return strcmp(word, first_word);
}

align* findAlignByName(align *aligns[], char *name) {
  int i=0;
  // printf("findAlignByName: %s\n", name);
  while(i<numseqs) {
    if (strstr(aligns[i]->seqs[0]->name, name)) {
      return(aligns[i]);
    }
    i++;
  }
  fprintf(stderr, "alignment not found for: %s", name);
  exit(2);
  return NULL;
}

int kk = 0;

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


align* processAlign(align *aligns[], align *a1, align *a2, int *index) {
  int hllindex;
  align *first, *second, *result, *uni;

  orderAligns(a1, a2, &first, &second, index, &hllindex);

  //  if (verbose
    //    printHLL(aligns[first->index]->hlls[hllindex]);  

  result = makeAlign(first, second, aligns[first->index]->hlls[hllindex], &uni);
  result->index = *index;

  freeHLLs(aligns[first->index]->hlls[hllindex]);
  aligns[first->index]->hlls[hllindex] = 0;    
  

  doRemapHLLs(aligns, uni, index, hllindex);

  doReanchorHLLs(aligns, index, hllindex);

  // if the constituent alignments were not simple alignments, free them
  freeAlign(uni); uni = 0;
  if (first->numseq > 1){ freeAlign(first); first = 0; }
  if (second->numseq > 1){ freeAlign(second); second = 0; }

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

  while (i < length) {

    if (rpntree[i]) {
      stack[sp++] = rpntree[i];
    }
    else {
      stack[sp-2] = processAlign(myaligns, stack[sp-2], stack[sp-1], &index);
      stack[--sp] = 0;      
      if(verbose) printTextAlign(stdout, stack[sp-1]);  
    }

    if (nested) {
      iterativeImprovement(stack[sp-1], rpntree, i);
    }

    i++;
  }
  return stack[sp-1];
}


char* buildTree (align *simalign[], float distances[MAX_SEQ][MAX_SEQ]) {
  char *names[MAX_SEQ];
  int namelens[MAX_SEQ];
  float max;
  int mli, mlj;
  int i, j;
  char *result, *temp;

  //  fprintf (stderr, "into build\n");

  for (i=0; i< numseqs; i++) {
    namelens[i] = strlen(simalign[i]->seqs[0]->name);
    names[i] = (char*) malloc ((namelens[i]+1) * sizeof (char));
    sscanf (simalign[i]->seqs[0]->name,"%s",names[i]); 
  }
  
  do {
    max = -1;
    for (i=0; i<(numseqs-1); i++) {
      for (j=i+1; j<numseqs; j++) {
	if (distances[i][j] > max) {
	  max = distances[i][j];
	  mli = i;
	  mlj = j;
	}
      }
    }
    if (max < 0)
      break;
    //    fprintf (stderr, "join! %d %d (score %f)\n", mli, mlj, distances[mli][mlj]);
    temp = (char*) malloc ((namelens[mli] + namelens[mlj] +4)* sizeof(char));
    sprintf(temp, "(%s %s)", names[mli], names[mlj]);

    //    fprintf (stderr, "%d(%d)+%d(%d)+3=%d(really %d)\n", namelens[mli],strlen(names[mli]),
    //	     namelens[mlj], strlen(names[mlj]), strlen(temp), namelens[mli]+namelens[mlj]+3);

    //    fprintf (stderr, "malloc gave %x\n", temp);
    //    fprintf (stderr, "new = %s\n", temp);
    //    fprintf (stderr, "done free1 %x\n", names[mli]);
    free (names[mli]);
    //    fprintf (stderr, "done free2 %x\n", names[mlj]);
    free (names[mlj]);
    names[mlj] = 0;
    names[mli] = result = temp;
    namelens[mli] = namelens[mli] + namelens[mlj] + 3;
    distances[mli][mlj] = -1;
    //    fprintf (stderr, "done concat\n");
    for (i=0; i < mli; i++) {
      //      fprintf (stderr, "h1\n");
      if (distances[i][mli] >= 0)
	distances[i][mli] = (distances[i][mli] + distances[i][mlj]) / 2;
      distances[i][mlj] = -1;
    }
    for (i=mli+1; i < mlj; i++) {
      //      fprintf (stderr, "h2\n");
      if (distances[mli][i] >= 0) 
	distances[mli][i] = (distances[mli][i] + distances[i][mlj]) / 2;
      distances[i][mlj] = -1;
    }
    for (i=mlj+1; i < numseqs; i++) {
      //      fprintf (stderr, "h3\n");
      if (distances[mli][i] >= 0) 
	distances[mli][i] = (distances[mli][i] + distances[mlj][i]) / 2;
      distances[mlj][i] = -1;
    }
    //    fprintf (stderr, "end of loop\n");
  } while (max >= 0);

  for (i=0; i< numseqs; i++) {
    if (names[i] != result)
      free (names[i]);
  }
  fprintf (stderr, "We built the tree: \"%s\"\n", result);
  return result;
}


char* graphCollapsal (align *simaligns[]) {
  float distances[MAX_SEQ][MAX_SEQ];
  int i, j;
  float sum = 0, length = 0;
  float score = 0, count = 0;
  hll* temp;

  for (i=0; i< MAX_SEQ; i++)
    for (j=0; j< MAX_SEQ; j++)
      distances[i][j] = -1;
  
  for (i=0; i<(numseqs-1); i++) {
    for (j=i+1; j<numseqs; j++) {
      sum = 0; count = 0;
      length = 0; score = 0;
      temp = simaligns[i]->hlls[j];
      while (temp) {
	sum += temp->score;
	length += (temp->seq1end - temp->seq1start);
	score += temp->score/(temp->seq1end - temp->seq1start);
	count += 1;
	temp = temp->next;
      }
      if (count != 0 && sum > 0) {
	//distances[i][j] = score/count;
	distances[i][j] = sum/length;
	//MIN2(simaligns[i]->seqs[0]->numsiglets, simaligns[j]->seqs[0]->numsiglets);
	fprintf (stderr, "Similarity %s and %s = %f\n",
		 simaligns[i]->seqs[0]->name, simaligns[j]->seqs[0]->name, distances[i][j]);
      }
      else 
	distances[i][j] = 0;
    }
  }
  return buildTree (simaligns, distances);
}

int parseParameters(int argc, char** argv, FileBuffer *files, char **treestr) {

  int i=1;

  if (argc < 3) {
    if (argc == 2)
      if (!strcmp(argv[1], "-version") || !strcmp(argv[1], "-Version")) {
        fprintf(stderr, "MLAGAN version %s\n", VER_NUM);
        exit(0);
      }
    usage();
    return 1;
  }
  while((argv[i][0]!='-')) {

    // Read in sequence files.
   
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
	  strcmp(argv[i], "-ext") || 	  strcmp(argv[i], "-scorematrix") ||
	  strcmp(argv[i], "-match") || strcmp(argv[i], "-mismatch") ||
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

    if (!strcmp(argv[i], "-nucmatrixfile")) {
      i++;
      if (i>=argc) {
	fprintf(stderr, "missing parameter specification for [-scorematrix.\n");
	return 1;
      }
      nucmatrixfile = argv[i];
      fprintf(stderr, "using nucmatrixfile value: %s\n", nucmatrixfile); 
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

  for (i=0; i<(numseqs-1); i++) {
    for (j=i+1; j<numseqs; j++) {
      table[i][j] = generateAnchors(files[i], files[j]);
    }
  }

  if (fastreject && !connectedGraph(table, numseqs)) {
    if (outfile != stdout)
      fclose (outfile);
    exit (0);
  }

  if (fastreject) {
    for (i=0; i<numseqs; i++) {
      for (j=i+1; j<numseqs; j++) {
	if (table[i][j])
	  table[i][j] = updateAnchorPos(table[i][j], files[i], files[j]);
	else
	  fprintf (stderr, "hmm\n");
      } 
    }
  }

  if (fastreject && !connectedGraph(table, numseqs)) {
    if (outfile != stdout)
      fclose (outfile);
    exit (0);
  }

  gapstart += gapcont;


  // Take all sequences and make simple alignments

  for (i=0; i<numseqs; i++) {
    if (fastreject) {
      if (files[i]->startpos > files[i]->endpos) {
	if (outfile != stdout)
	  fclose (outfile);
	exit (0);
      }
      seqs[i] = FileRead(files[i], 1, 0, VER_MLAGAN);
      


    }
    else 
      seqs[i] = FileRead(files[i], 0, 0, VER_MLAGAN);
    seqs[i]->index = i+1;
    myaligns[i]=simaligns[i]=mkSimAlign(seqs[i]);
    simaligns[i]->index = i;
  }


  // Find all pairwise anchors.

  for (i=0; i<(numseqs-1); i++) {
    for (j=i+1; j<numseqs; j++) {
      simaligns[i]->hlls[j]=table[i][j];
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

  fprintf(stderr, "\n****************************\n");
  fprintf(stderr, "gs: %d; ge: %d;\n", gapstart, gapend);
  fprintf(stderr, "gc: %d; gp: %d\n", gapcont, gapperseq);
  //fprintf(stderr, "match: %d; mismatch: %d\n", match, mismatch);
  fprintf(stderr, "overlap: %d; glwidth: %d\n", overlap, glwidth);
  fprintf(stderr, "\n****************************\n");
  
  if (notree) {
    treestr = graphCollapsal(myaligns);
  }

  //REMOVE the next line once debugged!!!
  //  exit(2);
  //End of remove

  i = 0;
  treeToRPN(treestr, stack, &i);
  
  final = procStack(stack, i, myaligns);
  

  if (postir) {
    final = iterativeImprovement(final, stack, i);
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














