#include "global.h"
#include <stdlib.h>
#include <stdio.h>

extern int indeces[256];

#define MAX2(x,y)   ( (x) >= (y) ? (x) : (y) )
#define MAX3(x,y,z)  MAX2(MAX2(x,y),z)


int ismatch(char a, char b) {
  return indeces[a] == indeces[b];
}

int matchscore (char a, char b) {
  if (a == b)
    return 4;
  return -3;
}

void reverse (char* a, int length) {
  char lft;
  int i;
  for (i=0; i < length/2; i++) {
    lft = a[i];
    a[i] = a[length-i-1];
    a[length-i-1] = lft;
  }
}

align* global(char* seq1, int start1, int end1, char* seq2, int start2, 
	      int end2, int gapopen, int gapext) {

  int mm = end2 - start2 + 1, score;
  int i,j,k,c, temp, lastdiag=0;
  int*  M = (int*) malloc (sizeof(int) * (end1-start1+1) * (end2 - start2+1));
  int*  N = (int*) malloc (sizeof(int) * (end1-start1+1) * (end2 - start2+1));
  int*  O = (int*) malloc (sizeof(int) * (end1-start1+1) * (end2 - start2+1));
  align* result = (align*) malloc (sizeof(align));
  char* almt = (char*) malloc ( sizeof(char) * ((end1-start1)+(end2-start2)+2));

  M[mm*0+0] = matchscore(seq1[start1],seq2[start2]);
  N[mm*0+0] = -1*gapopen;
  O[mm*0+0] = -1*gapopen;
  for (i = 1; i <= end1-start1; i++) {
    O[mm*i+0] = O[mm*(i-1)+0]-gapext;
    N[mm*i+0] = 0;
    M[mm*i+0] = O[mm*(i-1)+0]+matchscore(seq1[start1+i],seq2[start2]);
  }
  for (j = 1; j <= end2-start2; j++) {
    N[mm*0+j] = N[mm*0 + (j-1)]-gapext;
    O[mm*0+j] = 0;
    M[mm*0+j] = N[mm*0+(j-1)]+matchscore(seq1[start1],seq2[start2+j]);
  }
  for ( k = 2; k <= end1-start1; k++) {
    for (i = k-1, j = 1; (i > 0) && (j <= end2-start2); i--, j++) {
      N[mm*i + j] = MAX2(M[mm*(i-1)+j] - gapopen, N[mm*(i-1)+j] - gapext);
      O[mm*i + j] = MAX2(M[mm*i+(j-1)] - gapopen, O[mm*i+(j-1)] - gapext);
      M[mm*i + j] = MAX3(M[mm*(i-1)+(j-1)],N[mm*(i-1)+(j-1)],O[mm*(i-1)+(j-1)]) +
	matchscore(seq1[start1+i], seq2[start2+j]);
    }
  } 
  for ( k = 1; k <= end2-start2; k++) {
    for (j = k, i = end1-start1; (i>0) && (j <= end2-start2); j++, i--) {
      N[mm*i + j] = MAX2(M[mm*(i-1)+j] - gapopen, N[mm*(i-1)+j] - gapext);
      O[mm*i + j] = MAX2(M[mm*i+(j-1)] - gapopen, O[mm*i+(j-1)] - gapext);
      M[mm*i + j] = MAX3(M[mm*(i-1)+(j-1)],N[mm*(i-1)+(j-1)],O[mm*(i-1)+(j-1)]) +
	matchscore(seq1[start1+i], seq2[start2+j]);
    }
  }
  i = end1-start1; 
  j = end2-start2;
  c = 0;
  result->score = MAX3 ( M[mm*(i)+(j)], 
			 N[mm*(i)+(j)],
			 O[mm*(i)+(j)]);
  
  while(i >= 0 &&  j >= 0) {
    if (!i) {
      almt[c++] = ismatch(seq1[start1], seq2[start2+j]);
      for ( j = j -1; j >=0; j--,c++) {
	lastdiag = 0;
	almt[c] = DELETION;
      }
    }
    else if (!j) {
      almt[c++] = ismatch(seq1[start1+i], seq2[start2]);
      for ( i = i -1; i >=0; i--,c++) {
	almt[c] = INSERTION;
	lastdiag = 0;
      }
    }
    else {
      if (!lastdiag) {
	M[mm*i+j] = M[mm*i+j] - gapopen;  
	N[mm*i+j] = N[mm*i+j] - gapext;  
	O[mm*i+j] = O[mm*i+j] - gapext;  
      }
     
      temp = MAX3 ( M[mm*(i)+(j)], 
		    N[mm*(i)+(j)],
		    O[mm*(i)+(j)]);
      if (temp == N[mm*(i)+(j)]) {
	lastdiag = 0;
	almt[c++] = INSERTION;
	i--;
      }
      else if (temp == O[mm*(i)+(j)]) {
	lastdiag = 0;
	almt[c++] = DELETION;
	j--;
      }
      else if (temp == M[mm*(i)+(j)]) {
	lastdiag = 1;
	almt[c++] = ismatch(seq1[start1+i], seq2[start2+j]);
	i--; j--;
      }
    }
  }
  free(M);
  free(N);
  free(O);
  result->algnlen = c;
  reverse(almt,c);
  result->algn = almt;
  return result;
}

int printalign(char* seq1, int start1, int end1, char* seq2, int start2, 
	       int end2, align* myalign) {
  int s1=start1, s2=start2, c, k;
  int nm=0, nga=0, ngb=0, nlets=0;
  int hasst=0;
  for (c = 0; c < myalign->algnlen; c = c + 60) {
    for (k = c; (k < (c + 60)) && (k < myalign->algnlen); k++) {
      if (myalign->algn[k] != DELETION)
	printf("%c", seq1[s1++]);
      else {
	printf("-");
	if (hasst)
	  nga++;
      }
    } 
    printf("\n");
    for (k = c; (k < (c + 60)) && (k < myalign->algnlen); k++) {
      if (myalign->algn[k] == 1) {
	printf(":");
	nm++; 
	nlets++;
	hasst = 1; 
      }
      else {
	printf(" ");
	if (hasst) nlets++;
      }
    } 
    printf("\n");
    for (k = c; (k < (c + 60)) && (k < myalign->algnlen); k++) {
      if (myalign->algn[k] != INSERTION)
	printf("%c", seq2[s2++]);
      else {
	printf("-");
	if (hasst)
	  ngb++;
      }
    } 
    printf("\n\n");
  }
  printf("score = %d, nmatches = %d, nga=%d, ngb=%d nletters=%d, perc = %f\n",
	 myalign->score,nm,nga,ngb,nlets,(float)nm/(float)nlets);
  printf("\n");
}




