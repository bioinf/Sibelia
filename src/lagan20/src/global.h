#define INSERTION 2
#define DELETION 3

typedef struct align_res {
  int score;
  int algnlen;
  char* algn;
} align;

align* global(char* seq1, int start1, int end1, char* seq2, int start2, int end2,
	      int gapstart, int gapcont);

int printalign(char* seq1, int start1, int end1, char* seq2, int start2, int end2,
	      align* myalign);
