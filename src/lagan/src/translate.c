#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fchaos.h"
#include "translate.h"
#include "assert.h"

char toPeptide (char* dnaword, char revcomp) {
  int i, j, sum=0, mask = 0;
  char *table = 
    "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
  if (revcomp) mask = 3; /* Hacking... */
  for (i = 0; i < 3; i++) {
    sum*=4;
    switch (dnaword[(i^mask)-!!revcomp]) {
    case 'a': case 'A': sum+=(0^mask); break; 
    case 'c': case 'C': sum+=(1^mask); break; 
    case 'g': case 'G': sum+=(2^mask); break; 
    case 't': case 'T': sum+=(3^mask); break; 
    case 'n': case 'N': return 'X'; 
    default: 
      fprintf(stderr, "%d = %c: bad letter in sequence\n",i,dnaword[i^mask]);
      exit(1);
    }
  }
  return table[sum];
}


seq* transSeq(seq* theseq, int frame) {
  char* res;
  seq* resseq = (seq*) malloc(sizeof(seq));
  char revcomp = 0;
  int i, numXs = 0;

  assert (resseq);


  if (frame < 0 || frame > 5) {
    fprintf(stderr, "Valid frame numbers are 1-6\n");
    exit(1);
  }
  if (frame > 2) revcomp = 1;
  
  frame = frame % 3;
  resseq->numlets = (theseq->numlets-frame)/3;
  
  res = (char*) malloc((resseq->numlets+1)* sizeof(char));
  assert (res);

  /**
   * This was the error.
   */
  res[(theseq->numlets-frame)/3] = 0;
  for (i = 0;i < (theseq->numlets-frame)/3; i++) {
    res[i] = (!revcomp)?toPeptide(&theseq->lets[i*3+frame],0)
      :toPeptide(&theseq->lets[theseq->numlets-3*(i+1)-frame],1);
    if (res[i] == 'X') numXs++;
  } 

  resseq->numsiglets = resseq->numlets - numXs;  
  resseq->rptr = resseq->lets = res;
  resseq->name = (char*) malloc(strlen(theseq->name)+5);
  resseq->name[0] = 0;
  sprintf(resseq->name, "%s_f%c%d", theseq->name, (revcomp)?'-':'+', frame);
  return resseq;
}

/*
int main(int argc, char** argv) {
  printf("%s\n", transSeq(argv[1], strlen(argv[1]), 0));
  printf("%s\n", transSeq(argv[1], strlen(argv[1]), 1));
  printf("%s\n", transSeq(argv[1], strlen(argv[1]), 2));
  printf("%s\n", transSeq(argv[1], strlen(argv[1]), 3));
  printf("%s\n", transSeq(argv[1], strlen(argv[1]), 4));
  printf("%s\n", transSeq(argv[1], strlen(argv[1]), 5));
}
*/
