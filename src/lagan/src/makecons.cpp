/**
 * @file
 *
 * [TODO]
 *
 * @author Mikhail Soloviev
 * @date 31.03.2006
 * @version 1.0
 *
 */

//#include <iostream>
//#include <string>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>

using namespace std;

#define fastaRowLength 50
#define bufSize 2000

typedef char* pchar;

int isArg(char* key,int argc, char* argv[]) {
	for (int i=0;i<argc;i++) {
		if (strcmp(key,argv[i])==0) return 1;
	}
	return 0;
}

char* getArg(char* key,int argc, char* argv[]) {
	for (int i=0;i<argc;i++) {
		if (strcmp(key,argv[i])==0 && i<argc-1) return argv[i+1];
	}
   	fprintf(stderr,"ERROR: Parameter for option '%s' not specified\n",key);
   	exit(1);
	return NULL;
}

int trim(char* s) {
	int i=strlen(s);
	while (i>0 && (s[i-1]=='\n' || s[i-1]=='\r')) s[--i]='\0';
	return i;
}

FILE* openFile(char* path,char* mode) {
	FILE *f=fopen(path,mode);
	if (f==NULL) {
    	printf("ERROR: Failed open file: %s\n",path);
    	exit(1);
  	}
  	return f;
}

char* loadSeq(FILE *f,char* annot,int& seqLen) {
	char* seq=NULL;
	char buf[bufSize];
	int bufLen=0;
	seqLen=0;
  	while (!feof(f)) {
		buf[0]='\0';
		fgets(buf,bufSize,f);
		bufLen=trim(buf);
		if (bufLen>0) {
			if (buf[0]=='>') {
				strcpy(annot,buf);
				break;
			}
			else {
				if (seqLen==0) seq=(char*)malloc(sizeof(char)*bufLen);
				else seq=(char*)realloc(seq,sizeof(char)*(seqLen+bufLen));
				memcpy(&seq[seqLen],buf,bufLen);
				seqLen+=bufLen;
			}
		}
	}
	return seq;
}

void writeSeq(FILE *f,char* seq,int len) {
	int j=0;
	for (int i=0;i<len;i++,seq++) {
		fputc(*seq,f);
		j++;
		if (j==fastaRowLength) {
			j=0;
			fputc('\n',f);
		}
	}
	if (j>0) fputc('\n',f);
}

/*
char* makeCons(char* seq1,char* seq2,int len) {
	char* cons=seq1;
	char ch=' ';
	for (int i=0;i<len;i++,seq1++,seq2++) {
		if (*seq1=='-') {
			*seq1=*seq2;
		}
		else if (toupper(*seq1)=='N') {
			if (*seq2!='-') *seq1=*seq2;
		}
		else if (toupper(*seq1)==toupper(*seq2)) {
			if (islower(*seq1)) *seq1=*seq2;
		}
		else {
			ch=(rand()&1)?*seq1:*seq2;
			if (isupper(*seq1) || isupper(*seq2)) *seq1=toupper(ch); else *seq1=ch;
		}
	}
	return cons;
}
*/

/*
void makeCons(char seq1[],char seq2[],char cons[],int len) {
	for (int i=0;i<len;i++) {
		if (seq1[i]=='-') {
			cons[i]=seq2[i];
		}
		else if (seq2[i]=='-') {
			cons[i]=seq1[i];
		}
		else if (toupper(seq1[i])=='N') {
			cons[i]=seq2[i];
		}
		else if (toupper(seq2[i])=='N') {
			cons[i]=seq1[i];
		}
		else if (toupper(seq1[i])==toupper(seq2[i])) {
			cons[i]=isupper(seq1[i])?seq1[i]:seq2[i];
		}
		else {
			cons[i]=(rand()&1)?seq1[i]:seq2[i];
			if (isupper(seq1[i]) || isupper(seq2[i])) cons[i]=toupper(cons[i]);
		}
	}
}
*/

char dna[]={'N','A','C','G','T'};

int findMaxLetter(int count[],char* letter) {
	int max=0;
	int index=0;
	for (int i=1;i<5;i++) if (count[i]>max) max=count[i];
	for (int i=1;i<5;i++) if (count[i]==max) letter[index++]=dna[i];
	return index;
}

char makeConsLetter(char letter[],int proto) {
	int count[5];
	char maxLetter[5];
	int maxNumber;
	for (int j=0;j<5;j++) count[j]=0;
	for (int i=0;i<proto;i++) count[letter[i]]++;
	if (count[1]==0 && count[2]==0 && count[3]==0 && count[4]==0) {
		return 'N';
	}
	else {
		maxNumber=findMaxLetter(count,maxLetter);
		return maxNumber==1?maxLetter[0]:maxLetter[rand()%maxNumber];
	}
}

void makeCons(char cons[],char** seq,int proto,int len) {
	char letter[proto];
	for (int i=0;i<len;i++) {
		for (int j=0;j<proto;j++) {
			switch (toupper(seq[j][i])) {
				case 'A': letter[j]=1; break;
				case 'C': letter[j]=2; break;
				case 'G': letter[j]=3; break;
				case 'T': letter[j]=4; break;
				default:  letter[j]=0; break;
			}
		}
		cons[i]=makeConsLetter(letter,proto);
	}
}

int main (int argc,char* argv[]) {

	pchar seq[100];
	pchar cons=NULL;
	int len=0;
	int proto=0;
	char annot[2000];
	char nextAnnot[2000];

	srand((int)time(NULL));

	FILE *out=openFile(getArg("-o",argc,argv),"w");
	FILE *in=openFile(getArg("-i",argc,argv),"r");
	proto=atoi(getArg("-p",argc,argv));

	cons=loadSeq(in,annot,len);

  	while (!feof(in)) {
		for (int i=0;i<proto;i++) seq[i]=loadSeq(in,nextAnnot,len);

		cons=(char*)malloc(sizeof(char)*len);
		makeCons(cons,seq,proto,len);

		fprintf(out,"%s\n",annot);
		writeSeq(out,cons,len);

		strcpy(annot,nextAnnot);
		for (int i=0;i<proto;i++) free(seq[i]);
		free(cons);
	}
  	fclose(in);
  	fclose(out);
	return 0;
}
