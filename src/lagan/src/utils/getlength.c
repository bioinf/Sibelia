#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#define BUF_SIZE 1024

int main (int argc, char **argv){
  FILE *file;
  char buffer[BUF_SIZE], ch;
  int length = 0, i, done = 0, nread;

  if (argc != 2){
    fprintf (stderr, "Usage:\n\ngetlength seqfile\n");
    exit (1);
  }

  file = fopen (argv[1], "r"); assert (file);
  fgets (buffer, BUF_SIZE, file);
  while (!feof (file) && !done){
    nread = fread (buffer, 1, BUF_SIZE, file);
    for (i = 0; i < nread; i++){
      ch = buffer[i];
      if (ch == '>'){
	done = 1;
	break;
      }
      if (((ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z')) || ch == '.' || ch == '-')
	length++;
    }    
  }
  fclose (file);

  printf ("%d\n", length);
  return 0;
}










