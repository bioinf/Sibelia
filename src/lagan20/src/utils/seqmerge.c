#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

int main (int argc, char** argv){
  FILE *file;
  int i, written = 0;
  char buffer[1024], ch;

  if (argc == 1){
    fprintf (stderr, "Usage:\n\nseqmerge fasta_file1 fasta_file2 ...\n");
    exit (1);
  }

  for (i = 1; i < argc; i++){
    file = fopen (argv[i], "r"); assert (file);
    fgets (buffer, 1024, file);
    if (i == 1) printf ("%s", buffer);
    
    while (!feof (file)){
      ch = fgetc (file);
      if (ch == '>') break;
      if (isalpha (ch) || ch == '.' || ch == '-'){
	printf ("%c", ch);
	written++;
	if (written % 60 == 0) printf ("\n");
      }
    }
    fclose (file);
  }
  if (written ^ 60 != 0) printf ("\n");
}












