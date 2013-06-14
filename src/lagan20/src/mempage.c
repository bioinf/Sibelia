#include <stdio.h>
#include <stdlib.h>

#define MINPAGESIZE 1000000

typedef struct MemoryPage {
  char* memory;
  int size;
  int used;
  struct MemoryPage* next;
} mpage;


mpage* globalpage = 0;

void initMP(int pagesize) {
  mpage* newpage;
  if (pagesize < MINPAGESIZE)
    pagesize = MINPAGESIZE;

  newpage = (mpage*) malloc(sizeof(mpage));
  newpage->next = globalpage;
  globalpage = newpage;
  globalpage->memory = (char*) malloc (pagesize);
  globalpage->used = 0;
  globalpage->size = pagesize;
}

void* MPmalloc(int size) {
  void* tbr;
  if (globalpage->size - globalpage->used < size) {
    initMP(size);
  }
  tbr = globalpage->memory+ globalpage->used;
  globalpage->used += size;
  return tbr;
}

void* MPallfree() {
  mpage *n;
  while (globalpage) {
    free (globalpage->memory);
    n = globalpage;
    globalpage = globalpage->next;
    free(n);
  }
  initMP(0);
}

void* MPrealloc(void* prevptr, int prevsize, int newsize) {
  void* tbr = MPmalloc(newsize);
  memcpy(tbr, prevptr, prevsize);
  //  fprintf(stderr, "realloc returns %x instead of %x, (%d %d)\n", tbr, prevptr, prevsize, newsize);
  return tbr;
}
