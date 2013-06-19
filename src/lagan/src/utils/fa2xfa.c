#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

int begin, finish, seqIdx, seqExt, seqlen, numseqs, seqlen2, numseqs2;
int rcflag = 0;
char name[1024], name2[1024], **seqs, **seqs2;

char comp(char a) {
  if (!rcflag) return a;
  switch (a) {
  case 'A': case 'a': return 'T';
  case 'T': case 't': return 'A';
  case 'C': case 'c': return 'G';
  case 'G': case 'g': return 'C';
  case 'N': case 'n': return 'N';
  }
  fprintf (stderr, "bad letter to RC %c\n",a);
  exit(2);
}

int getLength (char *filename){
  FILE *file;
  char buffer[1024], ch;
  int length = 0;

  file = fopen (filename, "r"); assert (file);
  fgets (buffer, 1024, file);
  while (!feof (file)){
    ch = fgetc (file);
    if (ch == '>') break;
    if (((ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z')) || ch == '.' || ch == '-') length++;
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
      ch = toupper(ch);
      if (((ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z')) || (ch == '.') || (ch == '-')){
//	assert (i < (*seqlen));
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
  int i = 0, pos = 0, written = 0, j = 0;
  assert (seqExt >= 0 && seqExt < numseqs);
  name[0] = ' ';

  printf (">%d:%d-%d %c %s", seqIdx, begin+1, finish, (rcflag)?'-':'+', name);

  for (i = begin; i < finish; i++) {
    printf ("%c", comp(seqs[seqExt][(rcflag)?(finish+begin-i-1):i]));
    written++;
    if (written % 60 == 0) printf ("\n");
  }
    if (written % 60 != 0) printf ("\n");
}

int main (int argc, char** argv){
  int i;

  if (argc != 5 && !(argc == 6 && strcmp (argv[5], "-rc") == 0)){
    fprintf (stderr, "Usage:\n\nfa2xfa fasta_file begin end seqid [-rc]\n");
    exit (1);
  }

  seqExt = 0;
  begin = atoi (argv[2])-1;
  finish = atoi (strdup(argv[3]));
  seqIdx = atoi (argv[4]);
  if (argc == 6)
    rcflag = 1;
  seqlen2 = 0;

  readfile (argv[1], &seqlen, &numseqs, name, &seqs);

  print ();

  for (i = 0; i < numseqs; i++) free (seqs[i]);
  free (seqs);
}





