#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>

#define EXPAND 2

inline int max (int a, int b){ if (a > b) return a; return b; }
inline int min (int a, int b){ if (a < b) return a; return b; }

int getLength (char *filename){
  FILE *file;
  char buffer[1024], ch;
  int length = 0;

  file = fopen (filename, "r"); assert (file);
  fgets (buffer, 1024, file);
  while (!feof (file)){
    ch = fgetc (file);
    if (ch == '>') break;
    if (isalpha (ch) || ch == '.') length++;
  }
  fclose (file);

  return length;
}

int main (int argc, char **argv){
  FILE *file;
  int s1b, s1e, s2b, s2e, i;
  int S1B, S1E, S2B, S2E, ext, len1, len2;
  int m1b, m1e, m2b, m2e;
  float f;

  if (argc != 4){
    fprintf (stderr, "Usage:\n\ngetbounds anchfile seqfile1 seqfile2\n");
    exit (1);
  }

  file = fopen (argv[1], "r"); assert (file);
  len1 = getLength (argv[2]);
  len2 = getLength (argv[3]);

  m1b = m2b = 1000000000;
  m1e = m2e = -1000000000;
  while (!feof (file)){
    if (fscanf (file, "(%d %d)=(%d %d) %f\n", &s1b, &s1e, &s2b, &s2e, &f) == 5){
      m1b = min (m1b, s1b);
      m1e = max (m1e, s1e);
      m2b = min (m2b, s2b);
      m2e = max (m2e, s2e);
    }
  }
  m1e = len2 - m1e;
  m2e = len2 - m2e;
  fclose (file);
  file = fopen (argv[1], "r"); assert (file);

  i = 0;
  while (!feof (file)){
    if (fscanf (file, "(%d %d)=(%d %d) %f\n", &s1b, &s1e, &s2b, &s2e, &f) == 5){
      if (i == 0){
	S1B = max (s1b - m2b * EXPAND, 1);
	S1E = min (s1e + m2e * EXPAND, len1);
	S2B = max (s2b - m2b * EXPAND, 1);
	S2E = min (s2e + m2e * EXPAND, len2);
	i = 1;
      }
      else {
	S1B = min (S1B, max (s1b - m2b * EXPAND, 1));
	S1E = max (S1E, min (s1e + m2e * EXPAND, len1));
	S2B = min (S2B, max (s2b - m2b * EXPAND, 1));
	S2E = max (S2E, min (s2e + m2e * EXPAND, len2));
      }
    }
  }
  if (i == 0){
    S1B = 1;
    S1E = len1;
    S2B = 1;
    S2E = len2;
  }

  printf ("-s1 %d %d -s2 %d %d\n", S1B, S1E, 1, len2);
  
  fclose (file);
  return 0;
}

