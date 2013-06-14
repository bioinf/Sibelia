#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

int begin, finish, seqIdx, seqExt, seqlen, numseqs, seqlen2, numseqs2;
char name[1024], name2[1024], **seqs, **seqs2;

int getLength (char *filename){
  FILE *file;
  char buffer[1024], ch;
  int length = 0;

  file = fopen (filename, "r"); assert (file);
  fgets (buffer, 1024, file);
  while (!feof (file)){
    ch = fgetc (file);
    if (ch == '>') break;
    if (isalpha (ch) || ch == '.' || ch == '-') length++;
  }
  fclose (file);

  return length;
}

void readfile (char *filename, int *seqlen, int *numseqs, char *name, char ***seqs){
  FILE *file;
  char buffer[1024], ch;
  int i;

  *numseqs = 0;
  *seqlen = getLength (filename);
  strcpy (name, "");
  *seqs = (char **) malloc (sizeof (char *) * 1); assert (*seqs);
  (*seqs)[0] = (char *) malloc (sizeof (char) * (*seqlen));

  file = fopen (filename, "r"); assert (file);
  while (!feof (file)){
    i = 0;
    fgets (buffer, 1024, file);
    if (strlen (name) == 0) strcpy (name, buffer);
    if (feof (file)) break;
    (*numseqs)++;
    if (*numseqs > 1){
      *seqs = (char **) realloc (*seqs, sizeof (char *) * (*numseqs)); assert (*seqs);
      (*seqs)[*numseqs - 1] = (char *) malloc (sizeof (char) * (*seqlen)); assert ((*seqs)[*numseqs - 1]);
    }
        
    while (!feof (file)){
      ch = fgetc (file);
      if (ch == '>') break;
      if (isalpha (ch) || ch == '.' || ch == '-'){
	assert (i < (*seqlen));
	(*seqs)[*numseqs - 1][i] = ch;
	i++;
      }
    }
    if (ch == '>') ungetc (ch, file);
    assert (i == *seqlen);
  }
  fclose (file);
}

void print (void){
  int i = 0, pos = 0, pos2 = 0, written = 0, j = 0;

  while (pos <= finish && i < seqlen){
    if (isalpha (seqs[0][i])) pos++;
    if (isalpha (seqs[1][i])) pos2++;
    if (pos == finish){
      printf ("%d\n", pos2);
      break;
    }
    i++;
  }
}

int main (int argc, char** argv){
  int i;

  if (argc == 0){
    fprintf (stderr, "Usage:\n\ngetcontigpos multi_fasta_file finished_index\n");
    exit (1);
  }

  finish = atoi (strdup(argv[2]));

  readfile (argv[1], &seqlen, &numseqs, name, &seqs);
  print ();

  for (i = 0; i < numseqs; i++) free (seqs[i]);
  free (seqs);
}





