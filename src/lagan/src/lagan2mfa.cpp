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

void writeSeq(FILE *f,char* seq,int start,int end) {
	start--;
	end--;
	int j=0;
	for (int i=start;i<=end;i++) {
		fputc(seq[i],f);
		j++;
		if (j==fastaRowLength) {
			j=0;
			fputc('\n',f);
		}
	}
	if (j>0) fputc('\n',f);
}


int main (int argc,char* argv[]) {
	char buf[bufSize];

	char org0[1000];
	char name0[1000];
	int start0;
	int end0;
	char strand0;

	char org1[1000];
	char name1[1000];
	int start1;
	int end1;
	char strand1;

	char org2[1000];
	char name2[1000];
	int start2;
	int end2;
	char strand2;

	int proto=1;

	string id;
	string name;
	char* seq;

	FILE *out=openFile(getArg("-o",argc,argv),"w");
	FILE *chunk=openFile(getArg("-c",argc,argv),"w");
	FILE *in=openFile(getArg("-m",argc,argv),"r");
	proto=atoi(getArg("-p",argc,argv).c_str());
	readFaIndex(faIndex,getArg("-i",argc,argv));

  	while (!feof(in)) {
		buf[0]='\0';
		fgets(buf,bufSize,in);
		if (strlen(buf)==0) continue;

		sscanf(buf,"%s %s %d %d %c %s %s %d %d %c %s %s %d %d %c",
			org0,name0,&start0,&end0,&strand0,org1,name1,&start1,&end1,&strand1,org2,name2,&start2,&end2,&strand2);

		name=org0;
		name=name+"-anc"+name0;

		for (int n=1;n<=proto;n++) {
			id=name0;
			id=id+":"+itoa(n);
			seq=getFaIndexSeq(faIndex,id);
			fprintf(out,">%s\n",name.c_str());
			writeSeq(out,seq,start0,end0);
			free(seq);
		}
		end0=end0-start0+1;
		start0=1;

		fprintf(chunk,"%s %s %d %d %c %s %s %d %d %c %s %s %d %d %c\n",org0,name.c_str(),start0,end0,strand0,org1,name1,start1,end1,strand1,org2,name2,start2,end2,strand2);
	}
  	fclose(in);
  	fclose(out);
  	fclose(chunk);
	return 0;
}
