#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#define MAX_SEQS 63
#define MIN2(y,z)        ((y)<(z))?(y):(z)
#define MIN3(x,y,z)      MIN2((x),MIN2((y),(z)))
#define MIN4(w,x,y,z)    MIN2((w),MIN3((x),(y),(z)))


// Newick: (((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7):0.0;

// Takes a tree in newick format, builds an internal "tree" structure
// generates calls to other programs with correct weights



typedef struct sequence {
  char* seqname; 
  char* aligned;
  char* overlay;
  int alignlen;
  int overlaylen;
  int mynum;
} seq;


seq* allseqs[MAX_SEQS];
int numseqs;


char* dna_alpha = "ACGT";
char* valid_alpha = "ACGTN-";
char* DNA_PRINT;
char* DNA_LET;
char* NUM_ONES;

void init_consts() {
  int i;
  DNA_LET = (char*) malloc (sizeof(char) * 0x10);
  DNA_PRINT = (char*) malloc (sizeof(char) * 0x10);
  NUM_ONES = (char*) malloc (sizeof(char) * 0x10);

  for (i=0; i < 0x10; i++) {
    NUM_ONES[i] = DNA_LET[i] = DNA_PRINT[i] = -1;
  }

  DNA_LET[1] = 0;
  DNA_LET[2] = 1;
  DNA_LET[4] = 2;
  DNA_LET[8] = 3;
  DNA_PRINT[0] = 'N';
  DNA_PRINT[1] = 'A';
  DNA_PRINT[2] = 'C';
  DNA_PRINT[4] = 'G';
  DNA_PRINT[8] = 'T';
  DNA_PRINT[1|2] = 'M';
  DNA_PRINT[1|4] = 'R';
  DNA_PRINT[1|8] = 'W';
  DNA_PRINT[2|4] = 'S';
  DNA_PRINT[2|8] = 'Y';
  DNA_PRINT[4|8] = 'K';
  DNA_PRINT[1|2|4] = 'V';
  DNA_PRINT[1|2|8] = 'H';
  DNA_PRINT[1|4|8] = 'D';
  DNA_PRINT[2|4|8] = 'B';
  DNA_PRINT[1|2|4|8] = 'X';
  NUM_ONES[0] = 0;
  NUM_ONES[1] = 1;
  NUM_ONES[2] = 1;
  NUM_ONES[4] = 1;
  NUM_ONES[8] = 1;
  NUM_ONES[1|2] = 2;
  NUM_ONES[1|4] = 2;
  NUM_ONES[1|8] = 2;
  NUM_ONES[2|4] = 2;
  NUM_ONES[2|8] = 2;
  NUM_ONES[4|8] = 2;
  NUM_ONES[1|2|4] = 3;
  NUM_ONES[1|2|8] = 3;
  NUM_ONES[1|4|8] = 3;
  NUM_ONES[2|4|8] = 3;
  NUM_ONES[1|2|4|8] = 4;
}


seq* mk_seq() {
  seq* res = (seq*)malloc(sizeof(seq));
  res->seqname = 0;
  res->aligned = 0;
  res->overlay = 0;
  res->mynum = -1;
  return res;
}

int read_align(FILE* input, int target) {
  char* res = (char*) malloc(sizeof(char)*1);
  int i, ressize = 1, numread=0; 
  char temp[1024];
  char currchar, checkchar, *tt;

  if (feof(input)) {
    fprintf(stderr, "2COULDN'T READ ALIGNMENT\n");
    exit (2);
  }


  fgets(temp, 255, input);
  if (temp[0] != '>') {
    fprintf(stderr, "File is not in FASTA format!!\n");
    exit(1);
  }
  *(strchr(temp, '\n')) = 0;

  currchar = fgetc(input);

  while ((currchar != '>') && (currchar != EOF)) {
    if (!isspace(currchar)) {
      checkchar = toupper(currchar);
      if (!strchr(valid_alpha, checkchar)) {
	//        fprintf(stderr, "Warning: %d:%c skipped'\n", numread,currchar);
        currchar = 'N';
      }
      res[numread++] = currchar;
      if (numread >= ressize) {
        res=(char*)realloc(res, sizeof(char)*(ressize*=2));
      }
    }
    currchar = fgetc(input);
  }
  if (target >= 0) {
    allseqs[target]->seqname = malloc (strlen(temp)+1);
    strncpy(allseqs[target]->seqname, temp, strlen(temp)+1);
    allseqs[target]->aligned = res;
    allseqs[target]->alignlen = numread;
  }
  else {
    for (i = 0; i < numseqs; i++) {
      if (!strncmp(allseqs[i]->seqname, temp, strlen(temp))) {
	//	fprintf(stderr, "found %d\n",i);
	allseqs[i]->overlay = res;
	allseqs[i]->overlaylen = numread;
	break;
      }
    }

    if (i == numseqs) {
      fprintf(stderr, "seq %s not found!\n", temp);
      exit(2);
    }
  }
  if (currchar == '>') {
    ungetc(currchar, input);
    return 1;
  }
  return 0;
}

void read_align_file (char* filename) {

  FILE* input;
  if (!(input = fopen (filename, "r"))) {
    fprintf(stderr, "COULDN'T OPEN ALIGNMENT\n");
    exit (2);
  }
  while (read_align(input,numseqs++))
    ;
}


void read_sequences(int argc, char**argv) {
  char* filename;
  FILE* input;
  seq* myn;
  int i, j, kmer, breaker;
  int zz;

  for (i=2; i < argc; i++) {
    filename = argv[i];
    myn = 0;
    if (!(input = fopen (filename, "r"))) {
      fprintf(stderr, "COULDN'T OPEN SEQ %d %s\n",i,argv[i]);
      exit (2);
    }
    
    do {
      myn= allseqs[i-1];
      myn->mynum = i-1;
      zz = read_align(input,-1);
    } while (zz)
	;
  }
}
void overlayseq(int w) {
  int pos=0, i;
  for (i = 0; i < allseqs[w]->alignlen; i++) {
    if (allseqs[w]->aligned[i] != '-')
      allseqs[w]->aligned[i] = allseqs[w]->overlay[pos++];
  }
  fprintf(stderr, "check %d == %d\n",pos,allseqs[w]->overlaylen);
}


void overlay() {
  int i;
  for (i=0; i < numseqs; i++) {
    overlayseq(i);
  }
}

void printAlign() {
  int i,j;
  seq* a;
  for (j=0; j < numseqs; j++) {
    a = allseqs[j];
    fprintf(stdout, "%s", a->seqname);
    for (i=0; i < a->alignlen; i++) {
	if (!(i%60))
	  fprintf(stdout, "\n");
	//    fprintf(stdout, "%d:[%x]%c", i+1,a->aligned[i],DNA_PRINT[a->aligned[i]]);
	fprintf(stdout, "%c", a->aligned[i]);
      }
      fprintf(stdout, "\n");
  }
}


int main(int argc, char** argv) {
  char string_tree[16537]; //noone will ever need more :)))
  int moved, i;
  float ttree, test;
  
  //  fprintf(stderr, "Parsed tree\n");
  if (argc < 3) {
    fprintf(stderr, "Usage: overlay align.mfa seq1 [seq2].... > newalign.mfa\n");
    exit(2);
  }
  numseqs = 0;
  init_consts();


  for (i=0; i < MAX_SEQS; i++) {
    allseqs[i] = mk_seq();
  }

  //  ttree = get_outgroups(align_node, 0);
  //  fprintf(stdout, "ALIGN %s %s RES %s OUTS", align_node->lc->seqname, 
  //	  align_node->rc->seqname, align_node->seqname);
  //  for (i=0; i< numouts; i++) {
  //    fprintf(stdout, " %s %f", outgroups[i]->seqname, outdists[i]);
  //    test += outdists[i];
  //  }
  //  fprintf(stdout, "\n");
  
  read_align_file(argv[1]);
  read_sequences(argc, argv);
  overlay();
  printAlign();
  return 0;
}
