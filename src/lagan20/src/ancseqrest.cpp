/**
 * @file
 * Adds not aligned genome areas to ancestor FASTA file.
 *
 * Arguments:
 *
 * -b filename : block chunk mapping <br>
 * -g genomeindex : genome index, it refers to 2 files: genomeindex.ind and genomeindex.seq <br>
 * -n {1|2} : which genome is taken (1st or 2nd) from block chunk mapping <br>
 * -p proto : number of original species in genome
 * -o filename : ancestor fasta file, output sequence data to be appended here
 *
 * Block chunk mapping example:
 *
 * [TODO]
 *
 * Comment: [TODO].
 *
 *
 * @author Mikhail Soloviev
 * @date 23.05.2006
 * @version 1.0
 *
 */

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <stdlib.h>
#include <stdio.h>

using namespace std;

#define fastaRowLength 50
typedef char* pchar;

pchar seqData[100];
char seqStrand;

string itoa(int i) {
	char buf[20];
	sprintf(buf,"%d",i);
	return buf;
}

FILE* openFile(string path,char* mode) {
	FILE *f=fopen(path.c_str(),mode);
	if (f==NULL) {
    	fprintf(stderr,"ERROR: Failed open file: %s\n",path.c_str());
    	exit(1);
  	}
  	return f;
}

int isArg(char* key,int argc, char* argv[]) {
	for (int i=0;i<argc;i++) {
		if (strcmp(key,argv[i])==0) return 1;
	}
	return 0;
}

string getArg(char* key,int argc, char* argv[]) {
	for (int i=0;i<argc;i++) {
		if (strcmp(key,argv[i])==0 && i<argc-1) return argv[i+1];
	}
   	fprintf(stderr,"ERROR: Parameter for option '%s' not specified\n",key);
   	exit(1);
	return "";
}

string getArgAt(char* key,int index,int argc, char* argv[]) {
	for (int i=0;i<argc;i++) {
		if (strcmp(key,argv[i])==0 && i<argc-index) return argv[i+index];
	}
   	fprintf(stderr,"ERROR: Parameter for option '%s' not specified\n",key);
   	exit(1);
	return "";
}

struct Range {
	int start;
	int end;
	char strand;
};

struct Location {
	string genome;
	string name; // sequence name/id
	int start;
	int end;
	char strand;
};

struct ChunkMap {
	//int blockId;
	Location location[3];
};

vector<ChunkMap> chunkMap;

void loadChunkMap(string path) {
	char line[2000];
	char genome0[1000];
	char genome1[1000];
	char genome2[1000];
	char name0[1000];
	char name1[1000];
	char name2[1000];
	int tmp;
	FILE *in=openFile(path,"r");
  	while (!feof(in)) {
		line[0]='\0';
		fgets(line,2000,in);
		if (strlen(line)==0) continue;
		ChunkMap chunk;
		sscanf(line,"%s %s %d %d %c %s %s %d %d %c %s %s %d %d %c",
			genome0,name0,&chunk.location[0].start,&chunk.location[0].end,&chunk.location[0].strand,
			genome1,name1,&chunk.location[1].start,&chunk.location[1].end,&chunk.location[1].strand,
			genome2,name2,&chunk.location[2].start,&chunk.location[2].end,&chunk.location[2].strand);
		chunk.location[0].genome=genome0;
		chunk.location[1].genome=genome1;
		chunk.location[2].genome=genome2;
		chunk.location[0].name=name0;
		chunk.location[1].name=name1;
		chunk.location[2].name=name2;
		chunkMap.push_back(chunk);
	}
  	fclose(in);
}

void writeChunkSeq(FILE *out,string header,int start,int end,int protoStart,int protoEnd) {
	start--;
	end--;
	for (int p=protoStart;p<=protoEnd;p++) {
		fprintf(out,">%s\n",header.c_str());
		int j=0;
		for (int i=start;i<=end;i++) {
			fputc(seqData[p][i],out);
			j++;
			if (j==fastaRowLength) {
				j=0;
				fputc('\n',out);
			}
		}
		if (j>0) fputc('\n',out);
	}
}

void writeChunkGap(FILE *out,string header,int start,int end,int proto) {
	start--;
	end--;
	for (int p=1;p<=proto;p++) {
		fprintf(out,">%s\n",header.c_str());
		int j=0;
		for (int i=start;i<=end;i++) {
			fputc('-',out);
			j++;
			if (j==fastaRowLength) {
				j=0;
				fputc('\n',out);
			}
		}
		if (j>0) fputc('\n',out);
	}
}

Range noNext={0,0,'+'};

Range nextRange(int seqSize,Range prev) {
	Range next;
	prev.start--;
	prev.end--;
	next.start=prev.end+1;
	if (next.start>=seqSize) return noNext;
	while (seqData[1][next.start]=='*') {
		next.start++;
		if (next.start>=seqSize) return noNext;
	}
	next.end=next.start;
	while (next.end<seqSize && seqData[1][next.end+1]!='*') {
		next.end++;
	}
	next.start++;
	next.end++;
	return next;
}

void fillRange(int start,int end,int proto) {
	start--;
	end--;
	for (int p=1;p<=proto;p++) {
		for (int i=start;i<=end;i++) seqData[p][i]='*';
	}
}

void writeSeqRest(FILE *out,FILE *chunk,string ancestor,int seqSize,int& block,int genomeNumber,string descSeqName,int proto1,int proto2,string desc1,string desc2) {
	Range range=noNext;
	while ((range=nextRange(seqSize,range)).start!=0) {
		block++;
		string ancSeqName=ancestor+"-ancrest-"+itoa(genomeNumber)+"-"+itoa(block);
		if (genomeNumber==1) {
			writeChunkSeq(out,ancSeqName,range.start,range.end,1,proto1);
			writeChunkGap(out,ancSeqName,range.start,range.end,proto2);
			fprintf(chunk,"%s %s %d %d %c %s %s %d %d %c %s %s %d %d %c\n",
				ancestor.c_str(),ancSeqName.c_str(),1,(range.end-range.start+1),'+',
				desc1.c_str(),descSeqName.c_str(),range.start,range.end,seqStrand,
				desc2.c_str(),"-",0,0,'+');
		}
		else {
			writeChunkGap(out,ancSeqName,range.start,range.end,proto1);
			writeChunkSeq(out,ancSeqName,range.start,range.end,1,proto2);
			fprintf(chunk,"%s %s %d %d %c %s %s %d %d %c %s %s %d %d %c\n",
				ancestor.c_str(),ancSeqName.c_str(),1,(range.end-range.start+1),'+',
				desc1.c_str(),"-",0,0,'+',
				desc2.c_str(),descSeqName.c_str(),range.start,range.end,seqStrand);
		}
	}
	for (int i=1;i<=proto1+proto2;i++) free(seqData[i]);
}

struct FaRecord {
	string id;
	long offset;
	int length;
};

struct FaIndex {
	string id;
	FILE* file;
	map<string,FaRecord> record;
};

FaRecord readIndexRecord(FILE *ind) {
	FaRecord record;
	record.id="";
	char line[2000];
	char id[200];
	line[0]='\0';
	id[0]='\0';
	fgets(line,2000,ind);
	if (strlen(line)>0) {
		sscanf(line,"%s %ld %d",id,&record.offset,&record.length);
		record.id=id;
	}
  	return record;
}

FaIndex genomeIndex;

void openGenomeIndex(string genomePath) {
	genomeIndex.file=openFile(genomePath+".seq","r+");
	FILE *ind=openFile(genomePath+".ind","r");
  	while (!feof(ind)) {
		FaRecord record=readIndexRecord(ind);
		if (record.id.size()>0) genomeIndex.record[record.id]=record;
	}
  	fclose(ind);
}

char* readSeqBuf(FILE *seq,long offset,int length) {
	fseek(seq,offset,0);
	char* buf=(char*)malloc(length*sizeof(char));
	fread(buf,sizeof(char),length,seq);
	return buf;
}

void readGenomeSeq(string seqName,int& seqSize,int proto) {
	FILE *seq=genomeIndex.file;
	for (int i=1;i<=proto;i++) {
		string id=seqName+":"+itoa(i);
		FaRecord ind=genomeIndex.record[id];
		seqSize=ind.length;
		seqData[i]=readSeqBuf(seq,ind.offset,ind.length);
	}
}

int main (int argc,char* argv[]) {

	int block=0;
	string seqName="";
	string ancestor="";
	string desc1="";
	string desc2="";
	int seqSize=0;
	int proto=1;
	int proto1=1;
	int proto2=1;
	int genomeNumber=1;
	int first=1;

	FILE* out=openFile(getArg("-o",argc,argv),"w");
	FILE* chunk=openFile(getArg("-c",argc,argv),"w");
	loadChunkMap(getArg("-b",argc,argv));
	openGenomeIndex(getArg("-g",argc,argv));
	genomeNumber=atoi(getArg("-n",argc,argv).c_str());
	proto1=atoi(getArg("-p1",argc,argv).c_str());
	proto2=atoi(getArg("-p2",argc,argv).c_str());
	ancestor=getArg("-a",argc,argv);
	desc1=getArg("-d1",argc,argv);
	desc2=getArg("-d2",argc,argv);

	proto=genomeNumber==1?proto1:proto2;

	for (int i=0;i<chunkMap.size();i++) {
		Location loc=chunkMap[i].location[genomeNumber];
		seqStrand=loc.strand;
		if (loc.name=="-") continue;
		if (loc.name!=seqName) {
			if (!first) writeSeqRest(out,chunk,ancestor,seqSize,block,genomeNumber,seqName,proto1,proto2,desc1,desc2);
			else first=0;
			seqName=loc.name;
			seqStrand=loc.strand;
			readGenomeSeq(seqName,seqSize,proto);
		}
		fillRange(loc.start,loc.end,proto);
	}
	writeSeqRest(out,chunk,ancestor,seqSize,block,genomeNumber,seqName,proto1,proto2,desc1,desc2);

  	fclose(out);
  	fclose(chunk);
	return 0;
}
