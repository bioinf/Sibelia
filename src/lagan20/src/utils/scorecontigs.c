#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define MAX_SEQ 1024
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

#define CNTS_LEN 6
#define CNTS_A 0
#define CNTS_T 1
#define CNTS_C 2
#define CNTS_G 3
#define CNTS_N 4
#define CNTS_GAP 5

#define STATE_NULL 0
#define STATE_MATCH 1
#define STATE_MISMATCH 2
#define STATE_GAP 3
#define CACHE_SIZE 1000

int PEN_0_MIS, PEN_0_MTC, PEN_0_GAP;
int PEN_1_MIS, PEN_1_MTC, PEN_1_GAP;
int PEN_TO_0, PEN_TO_1;

char* alpha = "ATCGN-.";
double scoreMatch = 12;
double scoreMismatch = -4;
double scoreGapOpen = -80;
double cache[CACHE_SIZE];

typedef struct align_res {
  char *names[MAX_SEQ];
  int algnlen;
  int numseqs;
  char *data[MAX_SEQ];
} align;

typedef struct rangelist_res {
  int seqlen;
  int *score;
} rangelist;

int cntlets(FILE* input, int lettersonly) {
  int numread=0;
  char temp[1024];
  char currchar = '~';

  rewind (input);
  if (feof(input))
    return 0;
  fgets(temp, 1024, input);
  if (temp[0] != '>') {
    fprintf(stderr, "File is not in FASTA format!!\n");
    exit(1);
  }
  currchar = fgetc(input);
  while ((currchar != '>') && !feof (input)) {

    if (!isspace(currchar)) {
      currchar = toupper(currchar);
      if (!lettersonly || isalpha (currchar)){
	numread++;
      }
    }
    currchar = fgetc(input);
  }

  rewind(input);
  return numread;
}

int readseq (FILE *input, align *res){
  int numread = 0;
  char temp[1024], currchar, *write;

  if (feof (input)) return 0;
  fgets (temp, 1024, input);
  if (temp[0] != '>'){
    fprintf (stderr, "scorealign: File is not in FASTA format!!\n");
    exit (1);
  }
  res->names[res->numseqs] = (char*) malloc((strlen(temp))*sizeof(char));
  strcpy(res->names[res->numseqs], temp+1);
  *(strchr(res->names[res->numseqs], '\n')) = 0;

  write = res->data[res->numseqs] = (char *) malloc (sizeof (char) * res->algnlen); assert (write);

  currchar = fgetc (input);
  while (numread <= res->algnlen && (currchar != '>') && !feof (input)){
    if (!isspace (currchar)){
      currchar = toupper (currchar);
      if (!strchr(alpha, currchar)) currchar = 'N';
      write[numread++] = currchar;
    }
    currchar = fgetc (input);
  }

  if (currchar == '>'){
    ungetc (currchar, input);
  }

  if (numread != res->algnlen) {
    fprintf (stderr, "Sequence (%s) of different lengths (%d v. %d)!!\n", 
	     res->names[res->numseqs], numread, res->algnlen);
    exit(1);
  }
  return 1;
}
  
align *readMultial (char *filename){
  FILE *alfile;
  align *res;

  if (!(alfile = fopen (filename, "r"))){
    fprintf (stderr, "scorecontigs: couldn't open alignment file: %s\n", filename);
    exit (1);
  }

  res = (align *) malloc (sizeof (align)); assert (res);
  res->algnlen = cntlets (alfile, 0);
  res->numseqs = 0;
  
  while (readseq (alfile, res)) res->numseqs++;
  
  assert (res->numseqs == 2);
      
  fclose (alfile);

  return res;
}

inline int getstate (char c, char d){
  if (c == '-' || d == '-') return 2;
  if (c == 'N' || d == 'N') return 3;
  return c == d;
}

rangelist *getranges (char *filename, int offs){
  FILE *file;
  align *myal = readMultial (filename);
  rangelist *r = (rangelist *) malloc (sizeof (rangelist));
  int *scores[2], i, j, k, l, m, state, from0, from1, herescore;
  int *states, len, used, tot;
  char *traceback[2];
  
  assert (r);

  file = fopen (filename, "r"); assert (file);
  r->seqlen = cntlets (file, 1);
  len = cntlets (file, 0);
  for (i = 0; i < 2; i++){
    scores[i] = (int *) malloc (sizeof (int) * len); assert (scores[i]);
    traceback[i] = (char *) malloc (sizeof (char) * len); assert (traceback[i]);
  }

  for (i = 0; i < len; i++){
    state = getstate (myal->data[0][i], myal->data[1][i]);
    assert (i >= 0 && i < myal->algnlen);
    
    if (i <= 5){
      scores[0][i] = scores[1][i] = 0;
      traceback[0][i] = traceback[1][i] = 0;
    }
    else {

      // go to state 0
      herescore = (state == 0 ? PEN_0_MIS : (state == 1 ? PEN_0_MTC : (state == 2 ? PEN_0_GAP : 0)));
      from0 = scores[0][i-1] + herescore;
      from1 = scores[1][i-1] + herescore + PEN_TO_0;      
      if (from0 > from1){ scores[0][i] = from0; traceback[0][i] = 0; }
      else              { scores[0][i] = from1; traceback[0][i] = 1; }

      // go to state 1
      herescore = (state == 0 ? PEN_1_MIS : (state == 1 ? PEN_1_MTC : (state == 2 ? PEN_1_GAP : 0)));
      from0 = scores[0][i-1] + herescore + PEN_TO_1;
      from1 = scores[1][i-1] + herescore;      
      if (from0 > from1){ scores[1][i] = from0; traceback[1][i] = 0; }
      else              { scores[1][i] = from1; traceback[1][i] = 1; }
    }
  }

  states = (int *) malloc (sizeof (int) * len); assert (states);
  states[len - 1] = (scores[0][len - 1] > scores[1][len - 1]) ? 0 : 1;
  for (i = len - 2; i >= 0; i--) states[i] = traceback[states[i+1]][i+1];
  r->score = (int *) malloc (sizeof (int) * r->seqlen); assert (r->score);

  k = tot = used = 0;
  for (i = 0; i < len; i++){

    if (!states[i]){
      if (isalpha (myal->data[0][i])){
	r->score[k] = 0;
	k++;
      }
      continue;
    }

    used = 1;
    herescore = l = 0;
    
    for (j = i; j < len && states[j]; j++){
      if (isalpha (myal->data[0][j])) l++;
      state = getstate (myal->data[0][j], myal->data[1][j]);
      herescore += (state == 0 ? PEN_1_MIS : (state == 1 ? PEN_1_MTC : (state == 2 ? PEN_1_GAP : 0)));
    }
    tot += herescore;
    herescore /= l;

    //    fprintf (stderr, "%s: (%d %d) %d %d\n", filename, k + offs, k + l + offs, herescore, r->seqlen);
    for (m = k; m < k + l; m++) r->score[m] = herescore;

    k += l;
    i = j - 1;
  }

  //  printf ("%d\n", tot);

  free (states);

  for (i = 0; i < 2; i++){
    free (scores[i]);
    free (traceback[i]);
  }
  
  if (!used){
    free (r->score);
    free (r);
    return NULL;
  }
 
  return r;
}

inline int getdata (rangelist **ranges, int *offs, int j, int i){
  i -= offs[j];
  if (i >= 0 && i < ranges[j]->seqlen)
    return ranges[j]->score[i];
  return 0;
}


inline int match (rangelist **ranges, int numContigs, int i, int j, int *offs){
  int k;
  for (k = 0; k < numContigs; k++)
    if ((getdata (ranges, offs, k, i) != 0) != (getdata (ranges, offs, k, j) != 0)) return 0;
  return 1;
}

inline int allzeroes (rangelist **ranges, int numContigs, int pos, int *offs){
  int i;

  for (i = 0; i < numContigs; i++)
    if (getdata (ranges, offs, i, pos) != 0) return 0;
  return 1;
}

inline void print (int start, int end, int *score, int numContigs){
  int j;

  printf ("(%7d %7d)", start, end);
  for (j = 0; j < numContigs; j++) printf (" %7d", score[j]);
  printf ("\n");
}

void printRanges (rangelist **ranges, int numContigs, int seqLen, int *offs){
  int i, j, start = 0, end;
  int *score = (int *) malloc (sizeof (int) * numContigs);
  int *pattern = (int *) malloc (sizeof (int) * numContigs);

  assert (score);
  assert (pattern);
  
  printf ("numContigs = %d\n", numContigs);
  printf ("seqLen = %d\n", seqLen);

  for (i = 0; i < numContigs; i++) score[i] = 0;
  for (i = 0; i <= seqLen; i++)
    if (!allzeroes (ranges, numContigs, i, offs)) break;
  if (i > 0) print (0, i - 1, score, numContigs);

  start = end = i;
  while (i <= seqLen){
    if (i != seqLen && match (ranges, numContigs, start, i, offs)){
      end = i;
      for (j = 0; j < numContigs; j++){
	score[j] += getdata (ranges, offs, j, i);
      }
    }
    else if (i == seqLen || !allzeroes (ranges, numContigs, i, offs)){
      print (start, end, score, numContigs);
      for (j = 0; j < numContigs; j++) score[j] = 0;
      if (end < i - 1) print (end + 1, i - 1, score, numContigs);
      start = end = i;
    }
    i++;
  }

  free (score);
  free (pattern);
}

inline double scoregap (int gaplen){
  if (gaplen == 0) return 0;
  //return (gaplen - 1) * -1 - 50;
  return (log (gaplen) / log (10) + 1) * scoreGapOpen;
}

double scorealign (align *myal, int a, int b){
  int i, gaplen = 0;
  double score = 0;
  double best = 0;
  char c, d;


  // compensate for lagan bug
  for (i = 10; i < myal->algnlen; i++){
    c = myal->data[a][i]; d = myal->data[b][i];
    if (c == '-' && d == '-') continue;
    if (c == '-' || d == '-') gaplen++;
    else {
      if (gaplen != i){
	if (gaplen < CACHE_SIZE)
	  score += cache[gaplen];
	else
	  score += scoregap (gaplen);
      }
      gaplen = 0;
      if (c == d) score += scoreMatch;
      else score += scoreMismatch;
      if (score > best) best = score;
      if (score < 0) score = 0;
    }
  }

  return best;
}

void analyze (align *myal){

  int i, j, k;
  double score = 0;

  for (i = 0; i < CACHE_SIZE; i++) cache[i] = scoregap (i);

  for (i = 0; i < myal->numseqs; i++)
    for (j = i + 1; j < myal->numseqs; j++)
      score += scorealign (myal, i, j);

  printf ("%d\n", (int) score);
}
  
int main(int argc, char** argv) {
  FILE *filelist, *cfile;
  char contignames[MAX_SEQ][1024];
  rangelist *ranges[MAX_SEQ];
  int numseqs, i, j;
  int offs1[MAX_SEQ], offs2[MAX_SEQ], off[MAX_SEQ], num[MAX_SEQ];

  if (argc != 5) {
    fprintf(stderr, "Usage:\n\nscorecontigs file_list fasta_file contig_list cons_rate\n");
    exit (1);    
  }

  PEN_1_MIS = -(25 * atoi(argv[4])) / (101 - atoi (argv[4]));
  PEN_1_MTC = 25;
  PEN_1_GAP = PEN_1_MIS / 2;
  PEN_0_MIS = 0;
  PEN_0_MTC = 0;
  PEN_0_GAP = 0;
  PEN_TO_0 = -250; //-300;
  PEN_TO_1 = -350; //-400;

  if (!(filelist = fopen (argv[1], "r"))) {
    fprintf(stderr, "scorecontigs: Couldn't open alignment file: %s\n", argv[1]);
    exit (1);
  }

  numseqs = 0;
  while (!feof (filelist)){
    if (fscanf (filelist, "%d %d %d %s\n", &(num[numseqs]), &(offs1[numseqs]), &(offs2[numseqs]), &(contignames[numseqs])) == 4){
      numseqs++;
    }
  }
  fclose (filelist);

  if (numseqs == 0){
    fprintf (stderr, "scorecontigs: No contigs found.\n");
    exit (1);
  }

  cfile = fopen (argv[3], "w"); assert (cfile);
  j = 0;
  for (i = 0; i < numseqs; i++){
    ranges[j] = getranges (contignames[i], offs1[i]);
    if (ranges[j]){
      fprintf (cfile, "%d %d %d %s\n", num[i], offs1[i], offs2[i], contignames[i]);
      off[j] = offs1[i];
      j++;
    }
  }
  fclose (cfile);

  filelist = fopen (argv[2], "r"); assert (filelist);
  printRanges (ranges, j, cntlets (filelist, 1), off);
  fclose (filelist);
}
