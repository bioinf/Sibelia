#include <stdio.h>
#include <assert.h>

#define INTMAX (100000000)
#define INTMIN (-INTMAX)

int max (int a, int b){ if (a > b) return a; return b; }
int min (int a, int b){ if (a < b) return a; return b; }

int main (int argc, char **argv){
  FILE *file;
  int seq1begin = INTMAX, seq1end = INTMIN, seq2begin = INTMAX, seq2end = INTMIN;
  int a, b, c, d, e = 0;

  file = fopen (argv[1], "r"); assert (file);

  while (!feof (file)){
    if (fscanf (file, "(%d %d)=(%d %d) %*f\n", &a, &b, &c, &d) == 4){
      seq1begin = min (seq1begin, a);
      seq1end = max (seq1end, b);
      seq2begin = min (seq2begin, c);
      seq2end = max (seq2end, d);
      e++;
    }
  }

  fclose (file);

  if (!e)
    printf ("-1 -1 -1 -1\n");
  else
    printf ("%d %d %d %d\n", seq1begin, seq1end, seq2begin, seq2end);
}
