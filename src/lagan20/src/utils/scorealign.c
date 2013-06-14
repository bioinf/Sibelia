#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>

#define NUCLEOTIDE_MATRIX_FILE "nucmatrix.txt"
#define COLUMNS 60

int cons_rate = 0;
int doibounds = 0, doubounds = 0, leftbound, rightbound, pairseqlen;
int doregions = 0, docropxmfa = 0;
char **seqs;
int *seqid, *seqstart, *seqend;
char *seqdir, **seqcomment;
int numseqs, seqlen = -1;
int matchscore[256][256];
int gapopen = -1500, gapcont = -50;

inline int min (int a, int b){
  if (a < b) return a;
  return b;
}

inline int max (int a, int b){
  if (a > b) return a;
  return b;
}

inline int scoreMatch (char c, char d){
  if (c == '-' && d == '-') return 0;
  if (c == '-' || d == '-') return gapcont;
  return matchscore[(unsigned char) c][(unsigned char) d];
}

int conv2seqcoords (int pos, int i, int j){
  int alignpos = -1, pairpos = -1; 
  
  while (pairpos < pos && alignpos < seqlen){
    alignpos++;
    if (seqs[i][alignpos] != '-' || seqs[j][alignpos] != '-') pairpos++;
    if (alignpos >= seqlen){
      printf ("%d %d %d %d", pairpos, pos, alignpos, seqlen);
    }
    assert (alignpos < seqlen);
  }
  
  return alignpos+1;
}

#define CN 0
#define NC 1

int scorePair (char *seq1, char *seq2, int seqindex1, int seqindex2){
  int score[2][2];
  char *dad[2], *state;
  int i, j, CNscore, NCscore, left = pairseqlen, right = 1;
  
  for (i = 0; i < 2; i++){
    dad[i] = (char *) malloc (sizeof (char) * pairseqlen); assert (dad[i]);
    dad[i][0] = -1;
    score[i][0] = 0;
  }
  state = (char *) malloc (sizeof (char) * pairseqlen); assert (state);

  j = 0;
  for (i = 0; i < pairseqlen; i++){
    CNscore = score[CN][j];
    NCscore = score[NC][j] + gapopen;
    if (CNscore > NCscore){ score[CN][!j] = CNscore; dad[CN][i] = CN; }
    else                  { score[CN][!j] = NCscore; dad[CN][i] = NC; }
    score[CN][!j] += scoreMatch (seq1[i], seq2[i]);

    CNscore = score[CN][j] + gapopen;
    NCscore = score[NC][j];
    if (CNscore > NCscore){ score[NC][!j] = CNscore; dad[NC][i] = CN; }
    else                  { score[NC][!j] = NCscore; dad[NC][i] = NC; }

    j = !j;
  }

  i = pairseqlen - 1;
  j = (score[CN][j] > score[NC][j]) ? CN : NC;
  
  while (i >= 0){
    state[i] = j;
    assert (j == CN || j == NC);
    j = dad[j][i];
    i--;
  }

  j = 0;
  CNscore = 0;
  for (i = 0; i < pairseqlen; i++){
    if (state[i] == CN){
      if (!CNscore){
	CNscore = 1;
	if (doregions) printf ("Conserved region: %d ", i+1);
	left = min (left, i+1);
      }
      else if (i == pairseqlen - 1){
	if (doregions) printf ("%d\n", i+1);       
	right = max (right, i+1);
      }
      j++;
    }
    else if (CNscore){
      CNscore = 0;
      if (doregions) printf ("%d\n", i);
      right = max (right, i);
    }
  }

  if (j > 0){
    left = conv2seqcoords(left-1, seqindex1, seqindex2);
    right = conv2seqcoords(right-1, seqindex1, seqindex2);
    
    if (doibounds){
      leftbound = max (leftbound, left);
      rightbound = min (rightbound, right);
    }
    else if (doubounds){
      leftbound = min (leftbound, left);
      rightbound = max (rightbound, right);
    }
  }
  else {
    leftbound = 1;
    rightbound = seqlen;
  }
    
  for (i = 0; i < 2; i++) free (dad[i]);
  free (state);

  return j;
}

void project (char *orig1, char *orig2, char *dest1, char *dest2, int *length){
  int i, j;

  j = 0;
  for (i = 0; i < *length; i++){
    if (orig1[i] != '-' || orig2[i] != '-'){
      dest1[j] = orig1[i];
      dest2[j] = orig2[i];
      j++;
    }
  }
  *length = j;
}

int countleft (int pos, int i){
  int j, k;

  k = 0;			       
  for (j = 0; j < pos; j++)
    if (seqs[i][j] != '-') k++;

  return k;
}

int countright (int pos, int i){
  int j, k;

  k = 0;			       
  for (j = seqlen - 1; j > pos; j--)
    if (seqs[i][j] != '-') k++;

  return k;
}

void printXMFA (int score){
  int i, j, k;

  if (leftbound  > rightbound) {
    return;
  }

  if (seqid[0] == -1){
    for (i = 0; i < numseqs; i++){
      seqid[i] = i+1;
      seqstart[i] = 1;
      seqend[i] = countleft (seqlen, i);
      seqdir[i] = '+';
      strcpy (seqcomment[i], "");
    }
  }

  for (i = 0; i < numseqs; i++){
    if (seqcomment[i][strlen(seqcomment[i]) - 1] == '\n')
      seqcomment[i][strlen(seqcomment[i]) - 1] = '\0';

    printf (">%d:%d-%d %c %s\n", seqid[i],
	    seqstart[i] + countleft (leftbound-1, i), seqend[i] - countright(rightbound-1, i),
	    seqdir[i], seqcomment[i]);
    
    k = 0;
    for (j = leftbound - 1; j <= rightbound - 1; j++){
      printf ("%c", seqs[i][j]);
      k++;
      if (k % COLUMNS == 0) printf("\n");
    }
    if (k % COLUMNS != 0) printf("\n");
  }
  printf ("= score=%d\n", score);
}

void scoreAlign (){
  int i, j;
  int score = 0;
  char *u, *v;

  for (i = 0; i < numseqs - 1; i++){
    for (j = i + 1; j < numseqs; j++){
      pairseqlen = seqlen;
      u = (char *) malloc (sizeof (char) * seqlen); assert (u);
      v = (char *) malloc (sizeof (char) * seqlen); assert (v);
      project (seqs[i], seqs[j], u, v, &pairseqlen);
      score += scorePair (u, v, i, j);
      free (u);
      free (v);
    }
  }

  if (!doregions){
    if (doibounds || doubounds)
      if (docropxmfa){
	printXMFA(score);
      }
      else
	printf ("score=%d start=%d end=%d\n", score, leftbound, rightbound);
    else 
      printf ("%d\n", score);
  }
}

inline int issymbol (char ch){
  return ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'N' || ch == '.' || ch == '-';
}

void extractXMFAinfo (char *line, int *si, int *ss, int *se, char *sd, char **sc){
  int numread;

  *sc = malloc (sizeof (char) * 1024);  
  numread = sscanf (line, ">%d:%d-%d %c %s", si, ss, se, sd, *sc);

  if (numread < 4){
    *si = *ss = *se = -1;
    *sd = '~';
    strcpy (*sc, "");
  }
  else if (numread < 5){
    strcpy (*sc, "");
  }
}

char *getSequence (FILE *file, int *si, int *ss, int *se, char *sd, char **sc){
  int charsread = 0;
  int bufsize = 1;
  char *buffer;
  char prevch = '~';
  char line[1024];

  if (feof (file)) return NULL;
  fgets (line, 1024, file);
  if (line[0] == '='){
    return NULL;
  }

  extractXMFAinfo (line, si, ss, se, sd, sc);

  buffer = (char *) malloc (sizeof (char) * bufsize); assert (buffer);

  while (!feof (file)){
    buffer[charsread] = toupper (fgetc (file));

    if (buffer[charsread] == '>' || buffer[charsread] == '='){
      ungetc (buffer[charsread], file);
      break;
    }

    if (issymbol (buffer[charsread]))
      charsread++;
    
    if (charsread == bufsize){
      bufsize *= 2;
      buffer = (char *) realloc (buffer, sizeof (char) * bufsize);
    }
    
    prevch = buffer[charsread];
  }

  if (charsread == 0){
    free (buffer);
    return NULL;
  }

  if (seqlen == -1)
    seqlen = charsread;
  else {
    assert (seqlen == charsread);
  }

  return buffer;
}

int getSequences (FILE *file){
  char *newseq, sd, *sc;
  int i, si, ss, se;
  
  seqlen = -1;
  numseqs = 0;

  seqs = (char **) malloc (sizeof (char *) * 0);
  seqid = (int *) malloc (sizeof (int) * 0);
  seqstart = (int *) malloc (sizeof (int) * 0);
  seqend = (int *) malloc (sizeof (int) * 0);
  seqdir = (char *) malloc (sizeof (char) * 0);
  seqcomment = (char **) malloc (sizeof (char *) * 0);

  while (newseq = getSequence (file, &si, &ss, &se, &sd, &sc)){
    numseqs++;

    seqs = (char **) realloc (seqs, sizeof (char *) * numseqs);
    seqid = (int *) realloc (seqid, sizeof (int) * numseqs);
    seqstart = (int *) realloc (seqstart, sizeof (int) * numseqs);
    seqend = (int *) realloc (seqend, sizeof (int) * numseqs);
    seqdir = (char *) realloc (seqdir, sizeof (char) * numseqs);
    seqcomment = (char **) realloc (seqcomment, sizeof (char *) * numseqs);

    seqs[numseqs - 1] = newseq;
    seqid[numseqs - 1] = si;
    seqstart[numseqs - 1] = ss;
    seqend[numseqs - 1] = se;
    seqdir[numseqs - 1] = sd;
    seqcomment[numseqs - 1] = sc;
  }

  if (numseqs > 0) return 1;

  free (seqs);
  free (seqid);
  free (seqstart);
  free (seqend);
  free (seqdir);
  free (seqcomment);

  return 0;
}

int processSequences (FILE *file){
  int i, j;

  if (getSequences (file)){
    if (doibounds){
      leftbound = 0;
      rightbound = 1000000000;
    }
    else if (doubounds){
      leftbound = 1000000000;
      rightbound = 0;
    }

    scoreAlign();

    for (i = 0; i < numseqs; i++) free (seqs[i]);
    free (seqs);
    free (seqid);
    free (seqstart);
    free (seqend);
    free (seqdir);
    for (i = 0; i < numseqs; i++) free (seqcomment[i]);
    free (seqcomment);

    return 1;
  }
  return 0;
}

void calculateScoreMatrix(){
  char *alpha = "ATCG";
  int i, j;

  double p_ij = (double) cons_rate / 100.0;
  double match = log (p_ij / 0.25);
  double mismatch = log ((1 - p_ij) / 0.75);

  for (i = 0; i < strlen (alpha); i++){
    for (j = 0; j < strlen (alpha); j++){
      matchscore[(unsigned char) alpha[i]][(unsigned char) alpha[j]] = 
	(i == j) ? (int)(match * 100) : (int)(mismatch * 100);
    }
  }
  gapopen = (int)(-40 * match * 100);
}

void readScoreMatrix (char *filename){
  FILE *file;
  int i, j, k, numlets = 0;
  char lets[256], line[1024];  
  char *lagan_dir;

  lagan_dir = getenv ("LAGAN_DIR");
  if (!lagan_dir){
    fprintf (stderr, "Error: $LAGAN_DIR not set.\n");
    exit (1);
  }

  sprintf (line, "%s/%s", lagan_dir, filename);
  fprintf (stderr, "%s\n", line);

  file = fopen (line, "r"); assert (file);

  fgets (line, 1024, file);
  for (i = 0; i < strlen (line); i++){
    if (!isspace (line[i])){
      lets[numlets++] = line[i];
    }
  }

  for (i = 0; i < numlets; i++){
    fscanf (file, "%1s", &(line[0]));
    for (j = 0; j < numlets; j++){
      fscanf (file, "%d", &k);
      matchscore[(unsigned char) line[0]][(unsigned char) lets[j]] = k;
    }
  }

  fscanf (file, "%d%d", &gapopen, &gapcont);
  fclose (file);
}

void processFile (char *filename){
  FILE *file;
  int i, j;

  for (i = 0; i < 256; i++)
    for (j = 0; j < 256; j++)
      matchscore[i][j] = 0;

  if (cons_rate >= 0)
    calculateScoreMatrix();
  else
    readScoreMatrix (NUCLEOTIDE_MATRIX_FILE);

  file = fopen (filename, "r"); assert (file);
  while (!feof (file)){
    processSequences (file);
  }
  fclose (file);
}

int main (int argc, char **argv){
  int i;

  if (argc < 3 || argc > 6){
    // [-bounds seqidx]
    fprintf (stderr, "Usage: scorealign mfa_file cons_rate [-regions] [-ibounds | -ubounds [-cropxmfa]]\n");
    exit (1);
  }

  cons_rate = atoi (argv[2]);
  for (i = 3; i < argc; i++){
    if (strcmp (argv[i], "-cropxmfa") == 0)
      docropxmfa = 1;
    else if (strcmp (argv[i], "-ibounds") == 0)
      doibounds = 1;
    else if (strcmp (argv[i], "-ubounds") == 0)
      doubounds = 1;
    else if (strcmp (argv[i], "-regions") == 0)
      doregions = 1;
  }

  if (docropxmfa) assert (doibounds || doubounds);
  
  processFile (argv[1]);
  return 0;
}
