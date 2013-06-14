/**
 * @file
 * Compiles ancestor FASTA file using ansestor generation script.
 *
 * Arguments:
 *
 * -i filename : ansestor generation script <br>
 * -g genome genomeindex : genome index, genomeindex refers to 2 files: genomeindex.ind and genomeindex.seq <br>
 * -a alignmentindex : alignment index, alignmentindex refers to 2 files: alignmentindex.ind and alignmentindex.seq <br>
 * -o filename : output -- ancestor fasta file
 *
 * Ansestor generation script example:
 *
 * [TODO]
 *
 * Comment: [TODO].
 *
 *
 * @author Mikhail Soloviev
 * @date 31.03.2006
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

#include "util.cpp"
#include "faindex.cpp"

#define fastaRowLength 50

void revComp(char* seq,char* rev,long size) {
	rev+=size-1;
	for (long i=0;i<size;i++) {
		*rev=comp(*seq);
		seq++;
		rev--;
	}
}

void appendSeq(FILE *out,string header,string path) {
	fprintf(out,">%s\n",header.c_str());
	char buf[fastaRowLength+1];
	FILE *in=openFile(path,"r");
  	while (!feof(in)) {
		buf[0]='\0';
		fgets(buf,fastaRowLength,in);
		if (strlen(buf)>0) fprintf(out,"%s\n",buf);
	}
	fclose(in);
}

typedef char* pchar;
typedef FILE* pfile;
typedef pfile* ppfile;

struct Range {
	int start;
	int end;
};

struct AlignLocation {
	string org;
	string name; // sequence name/id
	int start;
	int end;
	char strand;
};

struct AlignMap {
	string id;
	map<string,AlignLocation> location; // string: orgId
	char strand;
};

map<string,AlignMap> alignMap; // string: alignId

void loadAlignMap(string path) {
	char line[2000];
	char id[1000];
	char name1[1000];
	char name2[1000];
	char org0[1000];
	char org1[1000];
	char org2[1000];
	AlignLocation loc0;
	AlignLocation loc1;
	AlignLocation loc2;

	FILE *in=openFile(path,"r");
  	while (!feof(in)) {
		line[0]='\0';
		fgets(line,2000,in);
		if (strlen(line)==0) continue;
		AlignMap aMap;
		sscanf(line,"%s %s %d %d %c %s %s %d %d %c %s %s %d %d %c",
			org0,id,&loc0.start,&loc0.end,&loc0.strand,
			org1,name1,&loc1.start,&loc1.end,&loc1.strand,
			org2,name2,&loc2.start,&loc2.end,&loc2.strand);
		loc0.org="0";
		loc1.org=org1;
		loc2.org=org2;
		loc0.name=id;
		loc1.name=name1;
		loc2.name=name2;
		aMap.id=id;
		aMap.strand=loc2.strand;
		aMap.location[loc0.org]=loc0;
		aMap.location[loc1.org]=loc1;
		aMap.location[loc2.org]=loc2;
		alignMap[aMap.id]=aMap;
	}
  	fclose(in);
}

// direct cut calculation: genome -> align, receives relative coord., returns absolute coord.

int calcCutStartLetter(char* seq,int start,int end,int relCut) {
	if (relCut==0) return start;
	int j=0;
	for (int i=start;i<=end;i++) {
		if (seq[i]!='-') j++;
		if (j==relCut) return i;
	}
	return start;
}

int calcCutEndLetter(char* seq,int start,int end,int relCut) {
	if (relCut==0) return end;
	int j=0;
	for (int i=end;i>=start;i--) {
		if (seq[i]!='-') j++;
		if (j==relCut) return i;
	}
	return end;
}

// reverse cut calculation: align -> genome, receives absolute coord., returns relative coord.

int revCalcCutStartLetter(char* seq,int start,int end,int absCut) {
	if (absCut==0) return 0;
	int j=0;
	for (int i=start;(i<=end && i<absCut);i++) {
		if (seq[i]!='-') j++;
	}
	return j;
}

int revCalcCutEndLetter(char* seq,int start,int end,int absCut) {
	if (absCut==0) return 0;
	int j=0;
	for (int i=end;(i>=start && i>absCut);i--) {
		if (seq[i]!='-') j++;
	}
	return j;
}

char* readSeqBuf(FILE *seq,long offset,int length) {
	fseek(seq,offset,0);
	char* buf=(char*)malloc(length*sizeof(char));
	fread(buf,sizeof(char),length,seq);
	return buf;
}

void writeSeqBuf(FILE *out,char* buf,int length,int sameStrand) {
	if (sameStrand) {
		fwrite(buf,sizeof(char),length,out);
	}
	else {
		char* rev=(char*)malloc(length*sizeof(char));
		revComp(buf,rev,length);
		fwrite(rev,sizeof(char),length,out);
		free(rev);
	}
	free(buf);
}

void writeSeq(FILE *out,FILE *seq,long offset,int length,int sameStrand) {
	char* buf=readSeqBuf(seq,offset,length);
	writeSeqBuf(out,buf,length,sameStrand);
}


/*OLD
void writeSeqCut(FILE *out,FILE *seq,long offset,int length,int sameStrand,int cutStart,int cutEnd) {
	offset+=cutStart;
	length-=cutStart+cutEnd;
	writeSeq(out,seq,offset,length,sameStrand);
}
*/

/*OLD
Range writeSeqCutLetter(FILE *out,FILE *seq,long offset,int length,int sameStrand,int cutStart,int cutEnd) {
	char* buf=readSeqBuf(seq,offset,length);
	cutStart=cutStartLetter(buf,length,cutStart);
	cutEnd=cutEndLetter(buf,length,cutEnd);
	length-=cutStart+cutEnd;
	memmove(buf,&buf[cutStart],length);
	writeSeqBuf(out,buf,length,sameStrand);
	Range r;
	r.start=cutStart;
	r.end=cutEnd;
	return r;
}
*/

map<string,FaIndex> genomeIndex;

void openGenomeIndex(string genomeName,string protoNumber,string genomePath) {
	FaIndex index;
	index.id=genomeName;
	index.proto=atoi(protoNumber.c_str());
	index.file=openFile(genomePath+".seq","r+");
	FILE *ind=openFile(genomePath+".ind","r");
  	while (!feof(ind)) {
		FaRecord record=readIndexRecord(ind);
		if (record.id.size()>0) index.record[record.id]=record;
	}
  	fclose(ind);
	genomeIndex[index.id]=index;
}

AlignLocation writeGenomeSeq(pfile out[],string orgName,int orgProto,string seqName,int start,int end,char strand) {
	FILE *seq=genomeIndex[orgName].file;
	for (int p=1;p<=orgProto;p++) {
		string recId=seqName+":"+itoa(p);
		FaRecord ind=genomeIndex[orgName].record[recId];
		writeSeq(out[p-1],seq,ind.offset+start-1,end-start+1,strand=='+');
	}
	AlignLocation loc;
	loc.org=orgName;
	loc.name=seqName;
	loc.start=start;
	loc.end=end;
	// TODO check
	loc.strand='+';
	return loc;
}

AlignLocation writeGenomeGap(pfile out[],string orgName,int orgProto,string seqName,int start,int end) {
	int size=end-start+1;
	char* buf=(char*)malloc(size*sizeof(char));
	memset(buf,'-',size);
	for (int p=1;p<=orgProto;p++) {
		fwrite(buf,sizeof(char),size,out[p-1]);
	}
	free(buf);
	AlignLocation loc;
	loc.org=orgName;
	loc.name=seqName;
	loc.start=start;
	loc.end=end;
	// TODO check
	loc.strand='+';
	return loc;
}

FaIndex alignIndex;

void openAlignIndex(string path) {
	alignIndex.file=openFile(path+".seq","r+");
	FILE *ind=openFile(path+".ind","r");
  	while (!feof(ind)) {
		FaRecord record=readIndexRecord(ind);
		if (record.id.size()>0) alignIndex.record[record.id]=record;
	}
  	fclose(ind);
}

int writeAlignSeq(pfile out1[],int proto1,pfile out2[],int proto2,string alignId,string orgName,char strand) {
	FILE *seq=alignIndex.file;
	AlignLocation loc=alignMap[alignId].location[orgName];
	AlignLocation loc0=alignMap[alignId].location["0"];
	int start=loc0.start-1;
	int length=loc0.end-loc0.start+1;
	FaRecord ind;
	for (int p=1;p<=proto1;p++) {
		string recId=alignId+":"+itoa(p);
		ind=alignIndex.record[recId];
		writeSeq(out1[p-1],seq,ind.offset+start,length,strand==loc.strand);
	}
	for (int p=1;p<=proto2;p++) {
		string recId=alignId+":"+itoa(proto1+p);
		ind=alignIndex.record[recId];
		writeSeq(out2[p-1],seq,ind.offset+start,length,strand==loc.strand);
	}
	return length;
}

/* not used anymore
AlignLocation writeAlignSeqCut(FILE *out,string alignId,string orgIndex,string orgName,char strand,int cutAlignStart,int cutAlignEnd) {
	FILE *seq=alignIndex.file;
	FaRecord ind=alignIndex.record[alignId+":"+orgIndex];
	AlignLocation loc=alignMap[alignId].location[orgName];
	writeSeqCut(out,seq,ind.offset,ind.length,strand==loc.strand,cutAlignStart,cutAlignEnd);

	// TODO -- find it via cutAlignStart,cutAlignEnd -- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	//loc.start+=cutStart;
	//loc.end-=cutEnd;
	return loc;
}
*/

// TODO check implementation when start implementing overlapping, compare with writeAlignSeq

/* OLD
AlignLocation writeAlignSeqCutLetterAlign(FILE *out,string alignId,string orgIndex,string orgName,char strand,int cutAlignStart,int cutAlignEnd) {
	FILE *seq=alignIndex.file;
	FaRecord ind=alignIndex.record[alignId+":"+orgIndex];
	AlignLocation loc=alignMap[alignId].location[orgName];

	// TODO -- optimize by excluding double reading the same sequence

	writeSeqCut(out,seq,ind.offset,ind.length,strand==loc.strand,cutAlignStart,cutAlignEnd);
	char* buf=readSeqBuf(seq,ind.offset,ind.length);
	loc.start+=reCutStartLetter(buf,ind.length,cutAlignStart);
	loc.end-=reCutEndLetter(buf,ind.length,cutAlignEnd);
	free(buf);
	return loc;
}
*/

// TODO check implementation when start implementing overlapping, compare with writeAlignSeq

/* OLD
AlignLocation writeAlignSeqCutLetter(FILE *out,string alignId,string orgIndex,string orgName,char strand,int cutStart,int cutEnd,int& cutAlignStart,int& cutAlignEnd) {
	FILE *seq=alignIndex.file;
	FaRecord ind=alignIndex.record[alignId+":"+orgIndex];
	AlignLocation loc=alignMap[alignId].location[orgName];
	Range r=writeSeqCutLetter(out,seq,ind.offset,ind.length,strand==loc.strand,cutStart,cutEnd);
	cutAlignStart=r.start;
	cutAlignEnd=r.end;
	loc.start+=cutStart;
	loc.end-=cutEnd;
	return loc;
}
*/

Range calcCutRangeLetter(char* seqBuf,int start,int end,int cutStartLength,int cutEndLength) {
	Range r;
	r.start=calcCutStartLetter(seqBuf,start,end,cutStartLength);
	r.end=calcCutEndLetter(seqBuf,start,end,cutEndLength);
	return r;
}

char* makeCons(string alignId,int protoStart,int protoEnd) {
	FILE *seqFile=alignIndex.file;
	char* cons=NULL;
	for (int p=protoStart;p<=protoEnd;p++) {
		string recId=alignId+":"+itoa(p);
		FaRecord ind=alignIndex.record[recId];
		char* buf=readSeqBuf(seqFile,ind.offset,ind.length);
		if (p==protoStart) {
			cons=(char*)malloc(ind.length*sizeof(char));
			memcpy(cons,buf,ind.length);
		}
		else {
			for (int i=0;i<ind.length;i++) if (buf[i]!='-') cons[i]=buf[i];
		}
		free(buf);
	}
	return cons;
}

AlignLocation writeAlignSeqCutLetter(pfile out[],int protoStart,int protoEnd,string alignId,string orgName,char strand,int cutStart,int cutEnd,int& cutAlignStart,int& cutAlignEnd) {
	FILE *seq=alignIndex.file;
	AlignLocation loc=alignMap[alignId].location[orgName];
	AlignLocation loc0=alignMap[alignId].location["0"];
	int start=loc0.start-1;
	int end=loc0.end-1;
	int length=loc0.end-loc0.start+1;
	FaRecord ind;
	char* cons=makeCons(alignId,protoStart,protoEnd);
	Range r=calcCutRangeLetter(cons,start,end,cutStart,cutEnd);
	for (int p=protoStart;p<=protoEnd;p++) {
		string recId=alignId+":"+itoa(p);
		ind=alignIndex.record[recId];
		writeSeq(out[p-1],seq,ind.offset+r.start,(r.end-r.start+1),strand==loc.strand);
	}
	cutAlignStart=r.start;
	cutAlignEnd=r.end;
	loc.start+=cutStart;
	loc.end-=cutEnd;
	free(cons);
	return loc;
}

AlignLocation writeAlignSeqCutLetterAlign(pfile out[],int protoStart,int protoEnd,string alignId,string orgName,char strand,int cutAlignStart,int cutAlignEnd) {
	FILE *seq=alignIndex.file;
	AlignLocation loc=alignMap[alignId].location[orgName];
	AlignLocation loc0=alignMap[alignId].location["0"];
	int start=loc0.start-1;
	int end=loc0.end-1;
	int length=loc0.end-loc0.start+1;
	FaRecord ind;
	char* cons=makeCons(alignId,protoStart,protoEnd);
	for (int p=protoStart;p<=protoEnd;p++) {
		string recId=alignId+":"+itoa(p);
		ind=alignIndex.record[recId];
		writeSeq(out[p-1],seq,ind.offset+cutAlignStart,(cutAlignEnd-cutAlignStart+1),strand==loc.strand);
	}
	loc.start+=revCalcCutStartLetter(cons,start,end,cutAlignStart);
	loc.end-=revCalcCutEndLetter(cons,start,end,cutAlignEnd);
	free(cons);
	return loc;
}

struct Command {
	char operation;
	string orgName;
	string seqName;
	string alignId1;
	string alignId2;
	int start;
	int end;
	int over1;
	int over2;
	char strand;
};

vector<Command> command;

void loadCommand(string path) {
	char line[1000];
	char orgName[100];
	char seqName[100];
	char alignId1[100];
	char alignId2[100];
	char operation;

	FILE *in=openFile(path,"r");
  	while (!feof(in)) {
		line[0]='\0';
		fgets(line,1000,in);
		if (strlen(line)==0) continue;
		Command com;
		operation=' ';
		orgName[100]='\0';
		seqName[100]='\0';
		alignId1[100]='\0';
		alignId2[100]='\0';
		com.over1=0;
		com.over2=0;
		sscanf(line,"%c ",&operation);
		if (operation=='g') {
			sscanf(line,"%c %s %s %d %d %c",&operation,orgName,seqName,&com.start,&com.end,&com.strand);
		}
		else if (operation=='s') {
			sscanf(line,"%c %s %s %c",&operation,alignId1,orgName,&com.strand);
		}
		else if (operation=='o') {
			sscanf(line,"%c %s %s %s %c %d %d",&operation,alignId1,alignId2,orgName,&com.strand,&com.over1,&com.over2);
		}
		else if (operation=='d') {
			sscanf(line,"%c %s %s %s %c",&operation,alignId1,alignId2,orgName,&com.strand);
		}
		else if (operation=='e') {
		}
		com.operation=operation;
		com.orgName=orgName;
		com.seqName=seqName;
		com.alignId1=alignId1;
		com.alignId2=alignId2;
		command.push_back(com);
	}
  	fclose(in);
}

void writeChunkLocation(FILE* blockChunk,AlignLocation loc) {
	fprintf(blockChunk,"%s %s %d %d %c",loc.org.c_str(),loc.name.c_str(),loc.start,loc.end,loc.strand);
}

void writeChunk(FILE* blockChunk,AlignMap chunk,string org[]) {
	writeChunkLocation(blockChunk,chunk.location[org[0]]);
	fprintf(blockChunk," ");
	writeChunkLocation(blockChunk,chunk.location[org[1]]);
	fprintf(blockChunk," ");
	writeChunkLocation(blockChunk,chunk.location[org[2]]);
	fprintf(blockChunk,"\n");
}

void openTmp(pfile tmp[],string outPath,int size,int offset) {
	for (int i=0;i<size;i++) {
		tmp[i]=openFile(outPath+"."+itoa(offset+i)+".tmp","w");
	}
}

void closeTmp(pfile tmp[],int size) {
	for (int i=0;i<size;i++) {
		fclose(tmp[i]);
	}
}

int main (int argc,char* argv[]) {

	string org[3];
	string ancOrg;

	map<string,ppfile> outtmp;
	map<string,string> other;
	map<string,string> orgIndex;
	map<string,int> proto;
	map<string,int> protoStart;

	AlignMap chunk;
	string header;
	int block=1;
	int multi=0;
	int start=0;
	int end=0;
	int ancProto=0;
	int	ancEnd=0;

	int cutAlignStart=0;
	int cutAlignEnd=0;

	string outPath=getArg("-o",argc,argv);
	FILE* out=openFile(outPath,"w");

	FILE* blockChunk=openFile(getArg("-b",argc,argv),"w");

	org[1]=getArg("-g1",argc,argv);
	org[2]=getArg("-g2",argc,argv);

	proto[org[1]]=atoi(getArgAt("-g1",2,argc,argv).c_str());
	proto[org[2]]=atoi(getArgAt("-g2",2,argc,argv).c_str());

	protoStart[org[1]]=1;
	protoStart[org[2]]=proto[org[1]]+1;

	ancProto=proto[org[1]]+proto[org[2]];

	loadAlignMap(getArg("-c",argc,argv));

	openAlignIndex(getArg("-a",argc,argv));

	openGenomeIndex(getArgAt("-g1",1,argc,argv),getArgAt("-g1",2,argc,argv),getArgAt("-g1",3,argc,argv));
	openGenomeIndex(getArgAt("-g2",1,argc,argv),getArgAt("-g2",2,argc,argv),getArgAt("-g2",3,argc,argv));

	ancOrg=org[1]+"_"+org[2];
	org[0]=ancOrg;

	chunk.location[org[0]].org=org[0];
	chunk.location[org[1]].org=org[1];
	chunk.location[org[2]].org=org[2];

	header=ancOrg+"-anc"+itoa(block);
	chunk.location[org[0]].name=header;
	chunk.location[org[0]].start=0;
	chunk.location[org[0]].end=0;

	other[org[1]]=org[2];
	other[org[2]]=org[1];

	orgIndex[org[1]]="1";
	orgIndex[org[2]]="2";

	pfile tmp1[proto[org[1]]];
	pfile tmp2[proto[org[2]]];
	outtmp[org[1]]=tmp1;
	outtmp[org[2]]=tmp2;

	openTmp(outtmp[org[1]],outPath,proto[org[1]],1);
	openTmp(outtmp[org[2]],outPath,proto[org[2]],proto[org[1]]+1);

	loadCommand(getArg("-i",argc,argv));

	// TODO: check and implement if necessary linking between s,d,o,g
	// in the same block, currently only d & o is linked

	for (int i=0;i<command.size();i++) {
		Command com=command[i];

		if (com.operation=='g') {
			multi=0;
			chunk.location[org[0]].start=ancEnd+1;
			ancEnd+=com.end-com.start+1;
			chunk.location[org[0]].end=ancEnd;
			chunk.location[com.orgName]=writeGenomeSeq(outtmp[com.orgName],com.orgName,proto[com.orgName],com.seqName,com.start,com.end,com.strand);
			chunk.location[other[com.orgName]]=writeGenomeGap(outtmp[other[com.orgName]],other[com.orgName],proto[other[com.orgName]],"-",com.start,com.end);
			writeChunk(blockChunk,chunk,org);
		}
		else if (com.operation=='s') {
			multi=0;
			chunk.location[org[0]].start=ancEnd+1;
			ancEnd+=writeAlignSeq(outtmp[org[1]],proto[org[1]],outtmp[org[2]],proto[org[2]],com.alignId1,com.orgName,com.strand);
			chunk.location[org[0]].end=ancEnd;
			chunk.location[com.orgName]=alignMap[com.alignId1].location[com.orgName];
			chunk.location[other[com.orgName]]=alignMap[com.alignId1].location[other[com.orgName]];
			writeChunk(blockChunk,chunk,org);
		}
		else if (com.operation=='d') {
			if (multi==0) multi=1; else multi=2;

			// align. 1
			if (multi==1) {
				chunk.location[org[0]].start=ancEnd+1;
				ancEnd+=writeAlignSeq(outtmp[org[1]],proto[org[1]],outtmp[org[2]],proto[org[2]],com.alignId1,com.orgName,com.strand);
				chunk.location[org[0]].end=ancEnd;
				chunk.location[com.orgName]=alignMap[com.alignId1].location[com.orgName];
				chunk.location[other[com.orgName]]=alignMap[com.alignId1].location[other[com.orgName]];
				writeChunk(blockChunk,chunk,org);
			}
			// genome between
			AlignLocation loc1=alignMap[com.alignId1].location[com.orgName];
			AlignLocation loc2=alignMap[com.alignId2].location[com.orgName];
			// TODO check possible overlap
			if (com.strand=='+') {
				start=loc1.end-1;
				end=loc2.start-1;
			}
			else {
				start=loc2.end+1;
				end=loc1.start-1;
			}
			// TODO -- currently it is assumed that seqName in the 1st and 2nd align. are the same -- check it !!!
			// see also the equivalent line below
			if (start<end) {
				chunk.location[org[0]].start=ancEnd+1;
				ancEnd+=end-start+1;
				chunk.location[org[0]].end=ancEnd;
				chunk.location[com.orgName]=writeGenomeSeq(outtmp[com.orgName],com.orgName,proto[com.orgName],loc1.name,start,end,com.strand);
				chunk.location[other[com.orgName]]=writeGenomeGap(outtmp[other[com.orgName]],other[com.orgName],proto[other[com.orgName]],"-",start,end);
				writeChunk(blockChunk,chunk,org);
			}
			else {
				printf("Warning: No gap between alignments %s and %s in %s (%d to %d)\n",
					com.alignId1.c_str(),com.alignId2.c_str(),com.orgName.c_str(),start,end);
			}
			// align. 2
			chunk.location[org[0]].start=ancEnd+1;
			ancEnd+=writeAlignSeq(outtmp[org[1]],proto[org[1]],outtmp[org[2]],proto[org[2]],com.alignId2,com.orgName,com.strand);
			chunk.location[org[0]].end=ancEnd;
			chunk.location[com.orgName]=alignMap[com.alignId2].location[com.orgName];
			chunk.location[other[com.orgName]]=alignMap[com.alignId2].location[other[com.orgName]];
			writeChunk(blockChunk,chunk,org);
		}
		// overlapping
		else if (com.operation=='o') {
			if (multi==0) multi=1; else multi=2;
			Command comNext=command[i+1];

			// align. 1
			if (multi==1) {
				if (com.strand=='+') {
					chunk.location[com.orgName]=writeAlignSeqCutLetter(outtmp[com.orgName],protoStart[com.orgName],proto[com.orgName],com.alignId1,com.orgName,com.strand,0,com.over1,cutAlignStart,cutAlignEnd);
					chunk.location[other[com.orgName]]=writeAlignSeqCutLetterAlign(outtmp[other[com.orgName]],protoStart[other[com.orgName]],proto[other[com.orgName]],com.alignId1,other[com.orgName],com.strand,cutAlignStart,cutAlignEnd);
					writeChunk(blockChunk,chunk,org);
				}
				else {
					chunk.location[com.orgName]=writeAlignSeqCutLetter(outtmp[com.orgName],protoStart[com.orgName],proto[com.orgName],com.alignId1,com.orgName,com.strand,com.over1,0,cutAlignStart,cutAlignEnd);
					chunk.location[other[com.orgName]]=writeAlignSeqCutLetterAlign(outtmp[other[com.orgName]],protoStart[other[com.orgName]],proto[other[com.orgName]],com.alignId1,other[com.orgName],com.strand,cutAlignStart,cutAlignEnd);
					writeChunk(blockChunk,chunk,org);
				}
			}
			// genome between
			AlignLocation loc1=alignMap[com.alignId1].location[com.orgName];
			AlignLocation loc2=alignMap[com.alignId2].location[com.orgName];
			// TODO check possible overlap
			if (com.strand=='+') {
				start=loc1.end-com.over1+1;
				end=loc2.start+com.over2-1;
			}
			else {
				start=loc2.end-com.over2+1;
				end=loc1.start+com.over1-1;
			}
			// TODO -- see TODO above
			if (start<end) {
				chunk.location[com.orgName]=writeGenomeSeq(outtmp[com.orgName],com.orgName,proto[com.orgName],loc1.name,start,end,com.strand);
				chunk.location[other[com.orgName]]=writeGenomeGap(outtmp[other[com.orgName]],other[com.orgName],proto[other[com.orgName]],"-",start,end);
				writeChunk(blockChunk,chunk,org);
			}
			else {
				printf("Warning: No gap between alignments %s and %s in %s (%d to %d)\n",
					com.alignId1.c_str(),com.alignId2.c_str(),com.orgName.c_str(),start,end);
			}
			// align. 2
			if (com.strand=='+') {
				chunk.location[com.orgName]=writeAlignSeqCutLetter(outtmp[com.orgName],protoStart[com.orgName],proto[com.orgName],com.alignId2,com.orgName,com.strand,com.over2,comNext.over1,cutAlignStart,cutAlignEnd);
				chunk.location[other[com.orgName]]=writeAlignSeqCutLetterAlign(outtmp[other[com.orgName]],protoStart[other[com.orgName]],proto[other[com.orgName]],com.alignId2,other[com.orgName],com.strand,cutAlignStart,cutAlignEnd);
				writeChunk(blockChunk,chunk,org);
			}
			else {
				chunk.location[com.orgName]=writeAlignSeqCutLetter(outtmp[com.orgName],protoStart[com.orgName],proto[com.orgName],com.alignId2,com.orgName,com.strand,comNext.over1,com.over2,cutAlignStart,cutAlignEnd);
				chunk.location[other[com.orgName]]=writeAlignSeqCutLetterAlign(outtmp[other[com.orgName]],protoStart[other[com.orgName]],proto[other[com.orgName]],com.alignId2,other[com.orgName],com.strand,cutAlignStart,cutAlignEnd);
				writeChunk(blockChunk,chunk,org);
			}
		}
		else if (com.operation=='e') {
			multi=0;
	  		closeTmp(outtmp[org[1]],proto[org[1]]);
  			closeTmp(outtmp[org[2]],proto[org[2]]);
			for (int i=1;i<=proto[org[1]];i++) appendSeq(out,header,outPath+"."+itoa(i)+".tmp");
			for (int i=1;i<=proto[org[2]];i++) appendSeq(out,header,outPath+"."+itoa(proto[org[1]]+i)+".tmp");
			openTmp(outtmp[org[1]],outPath,proto[org[1]],1);
			openTmp(outtmp[org[2]],outPath,proto[org[2]],proto[org[1]]+1);
			block++;
			header=ancOrg+"-anc"+itoa(block);
			chunk.location[org[0]].name=header;
			chunk.location[org[0]].start=0;
			chunk.location[org[0]].end=0;
			ancEnd=0;
		}
	}
	closeTmp(outtmp[org[1]],proto[org[1]]);
	closeTmp(outtmp[org[2]],proto[org[2]]);
  	fclose(out);
  	fclose(blockChunk);
	return 0;
}
