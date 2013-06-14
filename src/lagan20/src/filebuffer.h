#ifndef __FILEBUFFER_H
#define __FILEBUFFER_H

#include <stdio.h>

#ifndef MULTIAL__FLAG
#include "fchaos.h"
#else
#include "multial.h"
#endif

#define BUFFER_SIZE 1048576
#define VER_FCHAOS 0
#define VER_ORDER 1
#define VER_MLAGAN 2

struct FileBufferImplementation {
  FILE *data;
  char* filename;
  char buffer[BUFFER_SIZE];
  char *head, *tail;
  int startpos, endpos;
  //  int pos, len;
};

typedef struct FileBufferImplementation *FileBuffer;

FileBuffer FileOpen (const char *path);
int FileEOF (FileBuffer buf);
void FileGetS (char *buffer, int length, FileBuffer buf);
char FilePeekC (FileBuffer buf);
void FilePopC (FileBuffer buf);
void FileClose (FileBuffer buf);
seq* FileRead (FileBuffer buf, int start, int end, int version);

#endif
