#include "filebuffer.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#ifdef CHAOS__FLAG
char* alphabet = "ATCGNPCMHDEKRQSILVFYWX*";
#else
char* alphabet = "ATCGN-.";
#endif

FileBuffer FileOpen (const char *path){
  FileBuffer buf;
  FILE *data = fopen (path, "r");
  if (!data) return NULL;
  buf = (FileBuffer) malloc (sizeof (struct FileBufferImplementation));
  if (!buf) return NULL;
  buf->filename = (char*) path;
  buf->head = NULL;
  buf->tail = NULL;
  buf->startpos = 0; //100000000;
  buf->endpos = 100000000; //0;
  //buf->pos = BUFFER_SIZE;
  //buf->len = BUFFER_SIZE;
  buf->data = data;
  return buf;  
}

void FileUpdate (FileBuffer buf){
  if (buf->head >= buf->tail){
    buf->tail = buf->buffer + fread (buf->buffer, sizeof(char), BUFFER_SIZE, buf->data);
    buf->head = buf->buffer;
  }
}

int FileEOF (FileBuffer buf){
  FileUpdate (buf);
  return buf->head >= buf->tail && feof (buf->data);
}

void FileGetS (char *buffer, int length, FileBuffer buf){
  int a;

  for (a = 0; a < length && !FileEOF (buf); a++){
    buffer[a] = FilePeekC (buf);
    buf->head++;
    if (a + 1 < length && buffer[a] == '\n'){
      buffer[a + 1] = '\0';
      break;
    }
  }
}

char *FileGetLine (FileBuffer buf){
  int a = 0, length = 1;
  char *buffer = (char *) malloc (1 * sizeof(char));
  assert (buffer);

  while (!FileEOF (buf)){
    buffer[a] = FilePeekC (buf);
    buf->head++;
    if (buffer[a] == '\n'){
      buffer[a] = '\0';
      break;
    }
    a++;
    if (a == length){
      buffer = (char *) realloc (buffer, (length *= 2) * sizeof(char));
      assert (buffer);
    }
  }

  return buffer;
}

void FilePopC (FileBuffer buf){
  buf->head++;
}

char FilePeekC (FileBuffer buf){
  FileUpdate (buf);
  return *(buf->head);
  //  return buf->buffer[buf->pos];
}

void FileClose (FileBuffer buf){
  fclose (buf->data);
  free (buf);
}

seq* FileRead (FileBuffer buf, int start, int finish, int version){
  char* res = (char*) malloc(sizeof(char));
  int ressize = 1, numread = 0, i, numNs = 0;
  char *tempname, temp[256], currchar, *curr, *resend;
  seq* myseq = (seq*) malloc(sizeof(seq));


  if (FileEOF(buf))
    return 0;

  if (start == 1 && finish == 0) {
    start = buf->startpos;
    finish = buf->endpos;
    if (start == 0)
      start = 1;
  }

  tempname = FileGetLine (buf);
  if (tempname[0] != '>') {
    fprintf(stderr, "File is not in FASTA format!!\n");
    exit(1);
  }

  myseq->name = (char*) malloc((strlen(tempname))*sizeof(char));
  strcpy(myseq->name, tempname+1);
  if (strchr(myseq->name, '\n'))
    *(char *)(strchr(myseq->name, '\n')) = 0;

  free (tempname);

  for (i = 0; i < 256; i++){
    temp[i] = (strchr (alphabet, toupper ((char) i)) != 0) ?
      toupper((char) i) : 'N';
  }

  FileUpdate (buf);
  curr = res;
  resend = res + ressize;

  if (version == VER_ORDER || version == VER_MLAGAN){
    ressize = 2;
    numread = 1;
    if (version == VER_ORDER)
      res[0] = 0;
    else 
      res[0] = 'N';
    curr++;
  }

  while (buf->head < buf->tail || !feof (buf->data)){

    while (buf->head < buf->tail){
      currchar = *(buf->head);
      if (currchar == '>') goto outer;
      if (currchar != ' ' && currchar != '\n' && currchar != '\r' && 
	  currchar != '\t' && currchar != '\t' && currchar != '\v') {
	if (currchar == 'N') numNs++;
	*curr++ = temp[(int) currchar];
	if (curr >= resend) {
	  numread = curr - res;
	  res = (char *) realloc (res, sizeof(char) * (ressize *= 2));
	  curr = res + numread;
	  resend = res + ressize;
	}
      }
      buf->head++;
    }

    buf->tail = buf->buffer + fread (buf->buffer, sizeof(char), BUFFER_SIZE, buf->data);
    buf->head = buf->buffer;
  }
  
 outer:
  numread = curr - res;
  res[numread]=0;
  myseq->rptr = res;

  if (version == VER_FCHAOS){
    if (start > 0) {
      res[finish] = 0;
      res = &res[start-1];
      numread = finish-start+1;
    }
    myseq->numlets = numread;
  }
  else if (version == VER_ORDER){
    if (start > 0){
      res = &res[start-1];
      res[0] = 0;
      res[finish-start+2] = 0;
      numread = finish-start+2;
    }
    myseq->numlets = numread-1;
  }
  else if (version == VER_MLAGAN){
    if (start > 0 || finish > 0) {
      res[finish] = 0;
      res = &res[start-1];
      numread = finish-start+1;
    }
    myseq->numlets = numread;
    myseq->leftbound = start;
    myseq->rightbound = finish;
  }
  myseq->numsiglets = numread - numNs;
  myseq->lets = res;
  return myseq;
}
