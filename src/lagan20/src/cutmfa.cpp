/**
 * @file
 * Cuts Multi-FASTA file into parts using coordinate ranges
 * produced by supermap.
 *
 * Arguments:
 *
 * -i filename : input fasta file (containing only 1 sequence) <br>
 * -o filename : output fasta file <br>
 * -c filename : alignments' coordinate ranges (supermap output data) <br>
 * -s number   : take prototype organism sequences starting with number <br>
 * -e number   : take prototype organism sequences ending with number <br>
 * -u number   : which alignment coordinate range to use -- first or second,
 *               correspondingly number can be 1 or 2 <br>
 * -g {0|1}    : allow gaps <br>
 *
 * Alignments' coordinate range example:
 *
 * mouse-ENm001 1 12433   rat-ENm001 400 28619 + (DM, 13 aligns) <br>
 * mouse-ENm001 7001 14975   rat-ENm001 1 15303 + (M1, 1 aligns) <br>
 * mouse-ENm001 12872 51014   rat-ENm001 6891 71164 + (DM, 106 aligns)
 *
 * Comment: Only the first 6 fields are read, the rest can be anything.
 *
 * Resulted output example:
 *
 * >mouse-ENm001 <br>
 * GGACTCGTCGCAGTGCCTTGT <br>
 * TTTACTGTGCACTTCGCCTGG <br>
 * ACTGTCTACGCCATGCTTGAT <br>
 *
 * Comment: FASTA header contains sequence name (mouse-ENm001).
 *
 * @author Mikhail Soloviev
 * @date 05.04.2006
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

// TODO refactor in classes and normal make project

#include "util.cpp"
#include "faindex.cpp"

FaIndex faIndex;

void writeSeqDirect(FILE *out,char* seq,int start,int end,int gapped,int masked) {
	start--;
	end--;
	int j=0;
	for (int i=start;i<=end;i++) {
		if (gapped || seq[i]!='-') {
			fputc(masked?mask(seq[i]):seq[i],out);
			j++;
			if (j==fastaRowLength) {
				j=0;
				fputc('\n',out);
			}
		}
	}
	if (j>0) fputc('\n',out);
}

void writeSeqRevComp(FILE *out,char* seq,int start,int end,int gapped,int masked) {
	start--;
	end--;
	int j=0;
	for (int i=end;i>=start;i--) {
		if (gapped || seq[i]!='-') {
			fputc(masked?mask(comp(seq[i])):comp(seq[i]),out);
			j++;
			if (j==fastaRowLength) {
				j=0;
				fputc('\n',out);
			}
		}
	}
	if (j>0) fputc('\n',out);
}

void writeSeq(FILE *out,char* seq,int start,int end,int direct,int gapped,int masked) {
	if (direct) writeSeqDirect(out,seq,start,end,gapped,masked);
	else writeSeqRevComp(out,seq,start,end,gapped,masked);
}

int main (int argc,char* argv[]) {
	char buf[bufSize];
	char name[bufSize];
	int start;
	int end;
	char name2[bufSize];
	int start2;
	int end2;
	int count=0;
	char strand;

	int gapped=1;
	int useOrg=1;
	int protoStart=1;
	int protoEnd=1;
	int masked=0;

	string id;
	char* seq;

	FILE *out=openFile(getArg("-o",argc,argv),"w");
	FILE *in=openFile(getArg("-c",argc,argv),"r");

	readFaIndex(faIndex,getArg("-i",argc,argv));
	useOrg=atoi(getArg("-u",argc,argv).c_str());
	gapped=atoi(getArg("-g",argc,argv).c_str());
	protoStart=atoi(getArg("-s",argc,argv).c_str());
	protoEnd=atoi(getArg("-e",argc,argv).c_str());
	masked=atoi(getArg("-m",argc,argv).c_str());

  	while (!feof(in)) {
		buf[0]='\0';
		fgets(buf,bufSize,in);
		if (strlen(buf)==0) continue;
		sscanf(buf,"%s %d %d %s %d %d %c ",name,&start,&end,name2,&start2,&end2,&strand);
		if (useOrg==2) {
			strcpy(name,name2);
			start=start2;
			end=end2;
		}
		for (int n=protoStart;n<=protoEnd;n++) {
			id=name;
			id=id+":"+itoa(n);
			seq=getFaIndexSeq(faIndex,id);
			fprintf(out,">%s\n",name);
			writeSeq(out,seq,start,end,(useOrg==2 && strand=='-'),gapped,masked);
			free(seq);
		}
	}
  	fclose(in);
  	fclose(out);
	return 0;
}
