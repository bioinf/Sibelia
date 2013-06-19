#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void Add_Tick(char *line, int count, int length);
void Print_Lines(char *line1, char *line2, char *ticks1, char *ticks2,
  char *match);
int Usage(void);

char MyName[1024];

int main(int argc, char **argv) {
  FILE *infile = NULL;
  char bases[] = {'-', 'A', 'C', 'T', 'G', 'N'};
  char *seq1, *seq2;
  int seqsize=1, numread=0;
  int bp, base1, base2, i;
  seq1 = (char*) malloc(sizeof(char));
  seq2 = (char*) malloc(sizeof(char));
// parse my command line and open input file(s)

  if (argc < 2) return Usage();

  if ((strcmp(argv[1], "-") != 0) &&
      ((infile = fopen(argv[1], "r")) == NULL))
    return Usage();

  if (infile == NULL)
    infile = stdin;

  while (!feof(infile)) {
    if ((bp = getc(infile)) == EOF) {  // get next char
      break;
    }
    // decode bp char
    base1 = bp >> 4;
    base2 = bp & 0xf;
    seq1[numread] = bases[base1];
    seq2[numread] = bases[base2];
    numread++;
    if (numread >= seqsize) {
      seq1 = (char*) realloc(seq1, sizeof(char)* (seqsize *2));
      seq2 = (char*) realloc(seq2, sizeof(char)* (seqsize *2));
      seqsize *= 2;
    }
  }

  printf(">seq1");
  for (i = 0; i < numread; i++) {
    if (!(i%60))
      printf("\n");
    printf("%c", seq1[i]);
  }
  printf("\n>seq2");
  for (i = 0; i < numread; i++) {
    if (!(i%60))
      printf("\n");
    printf("%c", seq2[i]);
  }

  return 0;
}

int Usage() {
  fprintf(stderr, " \
Usage: %s { - | alignment_file }]\n",
    MyName);
  return 1;
}
