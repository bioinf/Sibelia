#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define MAX_SEQ 31
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

#define CNTS_LEN 6
#define CNTS_A 0
#define CNTS_T 1
#define CNTS_C 2
#define CNTS_G 3
#define CNTS_N 4
#define CNTS_GAP 5

double logs[MAX_SEQ+1];
double maxentr;
char* alpha = "ATCGN-";
int s1shift = 0, s2shift = 0;

typedef struct pair_ints {
  int s;
  int e;
} pair;

typedef struct align_res {
  char* names[MAX_SEQ];
  int algnlen;
  int numseq;
  int* algn;
  char* cnts[CNTS_LEN];
} align;

int cntlets(FILE* input) {
  int numread=0;
  char temp[256];
  char currchar = '~';

  if (feof(input))
    return 0;
  fgets(temp, 255, input);
  if (temp[0] != '>') {
    fprintf(stderr, "File is not in FASTA format!!\n");
    exit(1);
  }
  while ((currchar != '>') && (currchar != EOF)) {
    currchar = fgetc(input);
    if (!isspace(currchar)) {
      currchar = toupper(currchar);
      numread++;
    }
  }
  rewind(input);
  return numread-1;
}

int readseq(FILE* input, align* myal, int seqnum, int checksum) {
  int numread=0, help;
  char temp[256];
  char currchar;

  if (feof(input))
    return 0;
  fgets(temp, 255, input);
  if (temp[0] != '>') {
    fprintf(stderr, "File is not in FASTA format!!\n");
    exit(1);
  }
  myal->names[seqnum] = (char*) malloc((strlen(temp))*sizeof(char));
  strcpy(myal->names[seqnum], temp+1);
  *(strchr(myal->names[seqnum], '\n')) = 0;

  currchar = fgetc(input);
  while (numread <= checksum &&(currchar != '>') && (currchar != EOF)) {
    if (!isspace(currchar)) {
      currchar = toupper(currchar);
      if (!strchr(alpha, currchar)) {
	//	fprintf(stderr, "WARNING %c converted to N\n", currchar, alpha);
	currchar = 'N';
      }
      help = strchr(alpha, currchar)-alpha;
      myal->cnts[help][numread]++;
      if (help != CNTS_GAP) {
	myal->algn[numread] |= (1 << seqnum);
      }
      numread++;
    }
    currchar = fgetc(input);
  }
  if (currchar == '>')
    ungetc(currchar, input);
  if (numread != checksum) {
    fprintf(stderr, "Sequence (%s) of different lengths (%d v. %d)!!\n", 
	    myal->names[seqnum], numread, checksum);
    exit(1);
  }
  return 1;
}


align* readMultial(FILE* alfile) {
  int letcnt = cntlets(alfile), i, j;
  align* res = (align*)malloc (sizeof(align));
  res->algn = (int*) malloc (sizeof(int)* letcnt);
  for (j=0; j<CNTS_LEN; j++)
    res->cnts[j] = (char*) malloc (sizeof(char)* letcnt);
  for (i=0; i<letcnt; i++) {
    res->algn[i] = 0;
    for (j=0; j<CNTS_LEN; j++)
      res->cnts[j][i] = 0;
  }
  i = 0;
  while (readseq(alfile, res, i++, letcnt)) 
    ;

  res->numseq = i-1;
  res->algnlen = letcnt;
  return res;
}

inline int getScore (align* a, int i){
  return
    ((a->cnts[0][i] * (a->cnts[0][i] - 1)) +
     (a->cnts[1][i] * (a->cnts[1][i] - 1)) +
     (a->cnts[2][i] * (a->cnts[2][i] - 1)) +
     (a->cnts[3][i] * (a->cnts[3][i] - 1))) / 2;
}

void skipto (align *myal, int trgt, int *i, int* pos){
  int j;

  while (*i < trgt){
    for (j = 0; j < myal->numseq; j++)
      pos[j] += (myal->algn[*i] & (1 << j)) > 0;
    (*i)++;
  }
}

void print (align *myal, int *first, int *last, int len){
  int *start, *end, i, j, s = 0, e = 0;

  start = (int *) malloc (sizeof (int) * myal->numseq); assert (start);
  end = (int *) malloc (sizeof (int) * myal->numseq); assert (end);

  for (i = 0; i < myal->numseq; i++) start[i] = end[i] = 0;

  for (i = 0; i < len; i++){
    skipto (myal, first[i], &s, start);
    skipto (myal, last[i], &e, end);

    printf ("(%d %d) --> ", first[i] + s1shift, last[i] + s1shift);
    if (myal->numseq == 2){
      printf ("(%d %d)%s", start[0] + s1shift, end[0] + s1shift, (0 == myal->numseq - 1) ? "\n" : ", ");
      printf ("(%d %d)%s", start[1] + s2shift, end[1] + s2shift, (1 == myal->numseq - 1) ? "\n" : ", ");
    }
    else {
      for (j = 0; j < myal->numseq; j++){
	printf ("(%d %d)%s", start[0], end[0], (j == myal->numseq - 1) ? "\n" : ", ");
      }
    }

    // this is a hack -- can't handle multiple seq's
    /*
    for (j = 0; j < myal->numseq; j++){
      printf ("(%d %d)%s", start[j], end[j], (j == myal->numseq - 1) ? "\n" : ", ");
    }
    */
  }

  free (start);
  free (end);
}

void analyze (align *myal, int cutoff, int window){
  int *first, *last, size = 1, len = 0, i, score, count = 0;
  int runstart = -1, numpairs = myal->numseq * (myal->numseq - 1) / 2;

  window = MIN (window, myal->algnlen);
  first = (int *) malloc (size * sizeof (int)); assert (first);
  last = (int *) malloc (size * sizeof (int)); assert (last);

  score = 0;
  for (i = 0; i < window; i++)
    score += getScore (myal, i);

  if (score * 100 >= window * numpairs * cutoff) runstart = 0;
  for (i = 1; i <= myal->algnlen - window; i++){
    score += getScore (myal, i + window - 1) - getScore (myal, i - 1);

    if (score * 100 >= window * numpairs * cutoff){
      if (runstart == -1){
	if (len > 0 && last[len - 1] >= i)
	  runstart = first[--len];
	else
	  runstart = i;
      }
    }
    else if (runstart >= 0){
      first[len] = runstart;
      last[len++] = i + window - 1;
      runstart = -1;
	
      if (len == size){
	size *= 2;

	first = (int *) realloc (first, sizeof (int) * size); assert (first);
	last = (int *) realloc (last, sizeof (int) * size); assert (last);
      }
    }
  }

  if (runstart >= 0){
    first[len] = runstart;
    last[len++] = myal->algnlen - 1;
  }

  for (i = 0; i < len; i++){
    count += last[i] - first[i];
  }

  printf ("%d\n", count);
  print (myal, first, last, len);

  free (first);
  free (last);
}

int main(int argc, char** argv) {
  FILE *alignfile;
  align* myal;
  int i;

  if (argc != 4 && argc != 7) {
    fprintf(stderr, "usage:\ncstat multi_fasta_file cutoff window_size [-shift s1shift s2shift]\n");
    exit(1);    
  }
  if (!(alignfile = fopen(argv[1],"r"))) {
    fprintf(stderr, "couldnt open alignment file %s\n",argv[1]);
    return 2;
  }

  if (argc == 7){
    s1shift = atoi (argv[5]);
    s2shift = atoi (argv[6]);
  }

  myal = readMultial(alignfile);
  analyze (myal, atoi (argv[2]), atoi (argv[3]));
}
