#define MINPAGESIZE 256

typdef struct MemoryPage {
  void* memory;
  int size;
  int used;
  struct MemoryPage* next;
} mpage;


mpage globalpage;

void* initMP() {
  globalpage.memory = realloc (globalpage.memory, MINPAGESIZE);
  globalpage.used = 0;
  globalpage.size = MINPAGESIZE;
}

void* MPmalloc(int size) {
  void* tbr;
  while (globalpage.size - globalpage.used > size)
    globalpage.memory = realloc (globalpage.memory, (globalpage.size *=2));
  tbr = &(globalpage.memory[globalpage.used]);
  globalpage.used += size;
  return tbr;
}

void* MPallfree() {
  globalpage.memory = realloc (globalpage.memory, MINPAGESIZE);
  globalpage.used = 0;
  globalpage.size = MINPAGESIZE;
}


