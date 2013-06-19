#include "MultiSequence.h"
#include "SafeVector.h"
#include "Output.h"
#include <math.h>
#include <string.h>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <algorithm>

#define NUCLEOTIDE_MATRIX_FILE "nucmatrix.txt"
#define MAX_LINE_LENGTH 1024
#define CONS_RATE 70
#define INF 2000000000
#define CNTG_BRK_N 50

typedef SafeVector<int> vi;
typedef SafeVector<vi> vvi;
typedef SafeVector<vvi> vvvi;

MultiSequence seqs;
vvi matchScore (256, vi (256, 0));
vvi dad, score;
int gapopen, gapcont;
int NCtoNC = 0, NCtoCN = -1000, CNtoNC = -1000, CNtoCN = 0;

void readScoreMatrix (char *filename){
  FILE *file;
  int i, j, k, numlets = 0;
  char lets[256], line[1024];  
  char *lagan_dir;

  lagan_dir = getenv ("LAGAN_DIR");
  if (!lagan_dir){
    fprintf (stderr, "Error: $LAGAN_DIR not set.\n");
    exit (1);
  }

  sprintf (line, "%s/%s", lagan_dir, filename);
  fprintf (stderr, "%s\n", line);

  file = fopen (line, "r"); assert (file);

  fgets (line, 1024, file);
  for (i = 0; i < (int) strlen (line); i++){
    if (!isspace (line[i])){
      lets[numlets++] = line[i];
    }
  }

  for (i = 0; i < numlets; i++){
    fscanf (file, "%1s", &(line[0]));
    for (j = 0; j < numlets; j++){
      fscanf (file, "%d", &k);
      matchScore[(unsigned char) line[0]][(unsigned char) lets[j]] = k;
    }
  }

  fscanf (file, "%d%d", &gapopen, &gapcont);
  fclose (file);
}

void calculateScoreMatrix (int cons_rate){
  char *alpha = "ATCG";
  int i, j;

  for (int i = 0; i < 256; i++)
    for (int j = 0; j < 256; j++)
      matchScore[i][j] = 0;

  if (cons_rate == 0){
    readScoreMatrix (NUCLEOTIDE_MATRIX_FILE);
    return;
  }

  double p_ij = (double) cons_rate / 100.0;
  double match = log (p_ij / 0.25);
  double mismatch = log ((1 - p_ij) / 0.75);

  for (i = 0; i < (int) strlen (alpha); i++){
    for (j = 0; j < (int) strlen (alpha); j++){
      
      matchScore[(unsigned char) alpha[i]][(unsigned char) alpha[j]] = 
	(i == j) ? (int)(match * 100) : (int)(mismatch * 100);
    }
  }
  gapopen = (int)(-match * 750);
  gapcont = (int)(-match * 25);

  //  fprintf (stderr, "Using match=%d mismatch=%d gapopen=%d gapcont=%d...\n",
  //   (int)(match*100), (int)(mismatch*100), gapopen, gapcont);
}

#define NUM_STATES 2
#define NC 0
#define CN 1

void chooseBestOfTwo (int score1, int score2, int ptr1, int ptr2,
		      int &score, int &ptr){
  if (score1 >= score2){ score = score1; ptr = ptr1; }
  else                 { score = score2; ptr = ptr2; }
}

void chooseBestOfTwo (int score1, int score2, int &score){
  if (score1 >= score2){ score = score1; }
  else                 { score = score2; }
}

int scorePosition (char c, char d, int &isGap){
  if (c == '-' && d == '-') return 0;
  if (c == '-' || d == '-'){
    if (isGap) return gapcont;
    isGap = 1;
    return gapopen;
  }
  isGap = 0;
  return matchScore[(unsigned char) c][(unsigned char) d];
}

int rescoreRegion (Sequence &seq1, Sequence &seq2, int begin, int end){
  SafeVector<char>::iterator lets1 = seq1.getIterator();
  SafeVector<char>::iterator lets2 = seq2.getIterator();

  lets1 += begin - 1;
  lets2 += begin - 1;
  int isGap = 0;

  for (int i = 0; i < NUM_STATES; i++) score[i][begin-1] = dad[i][begin-1] = 0;

  for (int i = begin; i <= end; i++){
    chooseBestOfTwo (score[NC][i-1] + NCtoNC, score[CN][i-1] + CNtoNC, score[NC][i]);
    chooseBestOfTwo (score[NC][i-1] + NCtoCN, score[CN][i-1] + CNtoCN, score[CN][i]);
    score[CN][i] += scorePosition (*(++lets1), *(++lets2), isGap);
  }  
  
  chooseBestOfTwo (score[NC][end], score[CN][end], isGap);
  return isGap;
}

void getNucLabels (Sequence &seq1, Sequence &seq2, vi &nucLabels){
  SafeVector<char>::iterator lets1 = seq1.getIterator();
  SafeVector<char>::iterator lets2 = seq2.getIterator();
  int seqLen = seq1.getLength();
  int isGap = 0;

  nucLabels = vi (seqLen+1, 0);

  for (int i = 0; i < NUM_STATES; i++) score[i][0] = dad[i][0] = 0;

  for (int i = 1; i <= seqLen; i++){
    chooseBestOfTwo (score[NC][i-1] + NCtoNC, score[CN][i-1] + CNtoNC, NC, CN, score[NC][i], dad[NC][i]);
    chooseBestOfTwo (score[NC][i-1] + NCtoCN, score[CN][i-1] + CNtoCN, NC, CN, score[CN][i], dad[CN][i]);
    score[CN][i] += scorePosition (*(++lets1), *(++lets2), isGap);
  }

  chooseBestOfTwo (score[NC][seqLen], score[CN][seqLen], NC, CN, isGap, nucLabels[seqLen]);
  for (int i = seqLen - 1; i >= 1; i--){
    nucLabels[i] = dad[nucLabels[i+1]][i];
  }
}

int getSeqCoord (int seq, int pos){
  SafeVector<char>::iterator lets = seqs[seq].getIterator();
  int j = 0;
  
  for (int i = 1; i <= pos; i++)
    if (*(++lets) != '-') j++;
  
  return j;
}

void printCoordinates (int seq, int begin, int end){
  cout << seqs[seq].getID() << ":" << getSeqCoord(seq, begin) << "-" << getSeqCoord(seq, end) << " ";
}

int printRegion (int begin, int end){
  int score = 0;
  int numSeqs = seqs.getNumSeqs();

  for (int i = 0; i < numSeqs; i++){
    printCoordinates (i, begin, end);
    for (int j = i+1; j < numSeqs; j++){
      score += rescoreRegion (seqs[i], seqs[j], begin, end);
    }
  }
  cout << score << endl;
  return score;
}

void scoreAlign (){
  int numSeqs = seqs.getNumSeqs();
  int seqLen = seqs[0].getLength();
  vvvi nucLabels (numSeqs, vvi (numSeqs, vi()));

  for (int i = 0; i < numSeqs; i++){
    for (int j = i+1; j < numSeqs; j++){
      getNucLabels (seqs[i], seqs[j], nucLabels[i][j]);
    }
  }

  int begin = -1, end = -1, score = 0;
  for (int i = 1; i <= seqLen+1; i++){
    
    int conserved = 1;
    if (i == seqLen+1)
      conserved = 0;
    else {
      for (int j = 0; conserved && j < numSeqs; j++)
	for (int k = j+1; conserved && k < numSeqs; k++)
	  conserved = nucLabels[j][k][i];
    }

    if (conserved){
      if (begin == -1) 
	begin = i;
    }
    else {
      if (begin != -1){
	end = i-1;
	score += printRegion (begin, end);	
	begin = end = -1;
      }
    }    
  }

  cout << "= score=" << score << endl;
}

int countLets (SafeVector<char> &data){
  int ct = 0;
  for (int i = 0; i < (int) data.size(); i++){
    if (data[i] >= 'A' && data[i] <= 'Z' || data[i] >= 'a' && data[i] <= 'z')
      ct++;
  }
  return ct;
}

int findSplit (SafeVector<char> &data1, SafeVector<char> &data2, int overlap,
	       SafeVector<char> &data1a, SafeVector<char> &data2a){

  int offs1 = data1.size(), num1 = 0;
  for (int i = (int) data1.size() - 1; i >= 0; i--){
    if (overlap == 0) break;
    if (isalpha(data1[i])) num1++;
    if (num1 == overlap){
      offs1 = i;
      break;
    }
  }

  int offs2 = 0;
  num1 = 0;
  for (int i = 0; i < (int) data2.size(); i++){
    if (overlap == 0) break;
    if (isalpha(data2[i])) num1++;
    if (num1 == overlap){
      offs2 = i;
      break;
    }
  }

  SafeVector<int> score1 (overlap+1, 0);
  SafeVector<int> score2 (overlap+1, 0);

  int score = 0;
  for (int ct = 0,i=0; ct < overlap;i++){
    if (isalpha(data1[i+offs1])) ct++;
    score += (data1[i+offs1] == data1a[i+offs1]) ? 18 : -8;
    score1[ct] = score;
  }
  
  score = 0;
  for (int ct = 0,i=0; ct < overlap;i++){
    if (isalpha(data2[offs2-i])) ct++;
    score += (data2[offs2-i] == data2a[offs2-i]) ? 18 : -8;
    score2[ct] = score;
  }

  int j = 0, best = -1000000;
  for (int i = 0; i <= overlap; i++){
    if (score1[i] + score2[overlap-i] > best){
      best = score1[i] + score2[overlap-i];
      j = i;
    }
  }

  //  fprintf (stderr, "0 <= %d <= %d\n", j, overlap);
  
  return j;
}

template<class T>
int chopLeft (SafeVector<T> &data1, SafeVector<T> &data2, int num, bool inAlign){
  int num1 = 0, here = -1;

  if (inAlign)
    here = num - 1;
  else {
    for (int i = 0; i < (int) data1.size(); i++){
      if (num == 0) break;
      if (isalpha(data1[i])) num1++;
      if (num1 == num){
	here = i;
	break;
      }
    }
  }

  int chopped = here + 1;
  for (int i = here + 1; i < (int) data1.size(); i++){
    data1[i - chopped] = data1[i];
    data2[i - chopped] = data2[i];
  }

  data1.resize ((int) data1.size() - chopped);
  data2.resize ((int) data2.size() - chopped);

  return chopped;
}

template<class T>
int chopRight (SafeVector<T> &data1, SafeVector<T> &data2, int num, bool inAlign){
  int num1 = 0, here = data1.size();

  if (inAlign)
    here = data1.size() - num;
  else {
    for (int i = (int) data1.size() - 1; i >= 0; i--){
      if (num == 0) break;
      if (isalpha(data1[i])) num1++;
      if (num1 == num){
	here = i;
	break;
      }
    }
  }
    
  int ret = (int) data1.size() - here;
  data1.resize (here);
  data2.resize (here);

  return ret;
}

template<class T>
SafeVector<T> merge (SafeVector<T> &data1, SafeVector<T> &data2){
  SafeVector<T> temp;
  for (int i = 0; i < (int) data1.size(); i++) temp.push_back (data1[i]);
  for (int i = 0; i < (int) data2.size(); i++) temp.push_back (data2[i]);
  return temp;

}

int main (int argc, char **argv){
  FILE* outfile;
  
  if (argc < 2 || argc > 3){
    cerr << "Usage: Glue align.mfa \n" << endl;
    exit (1);
  }
  
  if (argc == 3) {
    if (!(outfile = fopen (argv[2], "w"))) {
      fprintf (stderr, "couldn't open %s for writing\n", argv[2]);
      exit(1);
    }

  }
  else outfile = stderr;

  //  calculateScoreMatrix (CONS_RATE);
  
  SafeVector<char> merged1, merged2;
  SafeVector<char> strand;
  SafeVector<int> merged1label, merged2label;
  int begin1 = 1, end1 = 1;

  ifstream data (argv[1]);
  int alignNum = 0;
  strand.push_back ('?'); // nothing for alignNum 0

  while (true){
    
    seqs = MultiSequence();
    seqs.addRawFromMFA (data);
    
    if (seqs.getNumSeqs() != 2) break;
    alignNum++;

    strand.push_back (seqs[1].getStrand());

    if (alignNum == 1){
      begin1 = seqs[0].getStartCoord();
      end1 = seqs[0].getEndCoord();
      merged1 = seqs[0].getData(); merged1label = SafeVector<int>((int) merged1.size(), 1);
      merged2 = seqs[1].getData(); merged2label = SafeVector<int>((int) merged2.size(), 1);
      continue;
    }

    int b1 = seqs[0].getStartCoord();
    int e1 = seqs[0].getEndCoord();

    SafeVector<char> seqs0;
    SafeVector<char> seqs1;

    seqs0 = seqs[0].getData();
    seqs1 = seqs[1].getData();

    SafeVector<int> seqs0label((int) seqs0.size(), alignNum);
    SafeVector<int> seqs1label((int) seqs1.size(), alignNum);

    int overlap = e1 - begin1 + 1;

    if (overlap > 0){
      int numLeft = findSplit (seqs0, merged1, overlap, seqs1, merged2);
      int numRight = overlap - numLeft;
      
      int choppedLeft = chopLeft (merged1, merged2, numLeft, false);
      int choppedRight = chopRight (seqs0, seqs1, numRight, false);

      chopLeft (merged1label, merged2label, choppedLeft, true);
      chopRight (seqs0label, seqs1label, choppedRight, true);
    }
    else if (overlap < 0){
      SafeVector<char> temp1 (-overlap, 'N');
      SafeVector<char> temp2 (-overlap, 'N');
      merged1 = merge (temp1, merged1);
      merged2 = merge (temp2, merged2);

      SafeVector<int> temp1label (-overlap, 0);
      SafeVector<int> temp2label (-overlap, 0);
      
      merged1label = merge (temp1label, merged1label);
      merged2label = merge (temp2label, merged2label);
    }

    merged1 = merge (seqs0, merged1);
    merged2 = merge (seqs1, merged2);

    merged1label = merge (seqs0label, merged1label);
    merged2label = merge (seqs1label, merged2label);

    //seqs[0].writeXMFAHeader(cerr);

    begin1 = b1;
    
    if (data.eof()) break;
    if (data.peek() == '=') data.ignore (MAX_LINE_LENGTH, '\n');
    if (data.eof()) break;
  }

  SafeVector<char> temp1 (begin1 - 1, 'N');
  SafeVector<char> temp2 (begin1 - 1, '-');

  for (int i = 0; i < min ((int) temp2.size(), CNTG_BRK_N); i++)
    temp2[i] = 'N';

  merged1 = merge (temp1, merged1);
  merged2 = merge (temp2, merged2);

  SafeVector<int> temp1label (begin1 - 1, 0);
  SafeVector<int> temp2label (begin1 - 1, 0);
  merged1label = merge (temp1label, merged1label);
  merged2label = merge (temp2label, merged2label);

  for (int i = 1; i <= alignNum; i++){
    int min1 = INF, max1 = 0, min2 = INF, max2 = 0;
    int pos1 = 0, pos2 = 0;
    for (int j = 0; j < (int) merged1label.size(); j++){
      if (isalpha(merged1[j])) pos1++;
      if (isalpha(merged2[j])) pos2++;
      
      if (merged1label[j] == i){
	min1 = min (min1, pos1);
	max1 = max (max1, pos1);
      }
      if (merged2label[j] == i){
	min2 = min (min2, pos2);
	max2 = max (max2, pos2);
      }
    }

    //[FASTA line for this contig in the original sequence file]
    //n baseFrom baseTo mergedFrom mergedTo startChop endChop {+,-} score secFrom secTo
    fprintf (outfile, "Align %d\n", i);
    if (min1 == INF)
      fprintf (outfile, "%d was cropped completely.\n", i);
    else
      fprintf (outfile, "%d %d %d 0 0 0 0 %c 0 %d %d\n", i, min1, max1, strand[i], min2, max2);
  }
  
  printMFA (cout, merged1, string ("first"), 60);
  printMFA (cout, merged2, string ("second"), 60);
}
