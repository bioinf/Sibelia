#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

char* alpha = "ATCGN";

typedef struct Sequence {
  char* lets;
  int numlets;
  char* name;
  char* rptr;
} seq;

char comp(char c) {
  switch(c) {
  case 'A': return 'T'; 
  case 'T': return 'A'; 
  case 'C': return 'G';
  case 'G': return 'C'; 
  case 'N': return 'N';
  case 'a': return 't';
  case 't': return 'a'; 
  case 'c': return 'g';
  case 'g': return 'c'; 
  case 'n': return 'n';
  default: return c;
  }
}

int main (int argc, char **argv){
  char* res = (char*) malloc(sizeof(char));
  int ressize = 1, numread = 0, i;
  char temp[256];
  char currchar;

  if (feof(stdin))
    return 0;
  fgets(temp, 255, stdin);
  if (temp[0] != '>') {
    fprintf(stderr, "File is not in FASTA format!!\n");
    exit(1);
  }
  *(strchr(temp,'\n')) = 0;
  //  strcat (temp, "(-)");
  printf ("%s\n", temp);

  currchar = fgetc(stdin);
  while ((currchar != '>') && (currchar != EOF)) {
    if (!isspace(currchar)) {
      res[numread++] = comp (currchar);
      if (numread >= ressize) {
	res=(char*)realloc(res, sizeof(char)*(ressize*=2)); 
      }
    }
    currchar = fgetc(stdin);
  }
  res[numread]=0;
  i = 0;
  while (--numread >= 0){
    putchar (res[numread]);
    i++;
    if (i % 60 == 0){
      putchar ('\n');
      i = 0;
    }
  }
  if (i != 0) putchar ('\n');
  free (res);
  return 0;
}
