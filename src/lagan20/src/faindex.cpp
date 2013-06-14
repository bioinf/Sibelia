struct FaRecord {
	string id;
	long offset;
	int length;
};

struct FaIndex {
	string id;
	int proto;
	FILE* file;
	map<string,FaRecord> record;
};

FaRecord readIndexRecord(FILE *ind) {
	FaRecord record;
	record.id="";
	char line[1000];
	char id[100];
	line[0]='\0';
	id[0]='\0';
	fgets(line,1000,ind);
	if (strlen(line)>0) {
		sscanf(line,"%s %ld %d",id,&record.offset,&record.length);
		record.id=id;
	}
  	return record;
}

void readFaIndex(FaIndex& faIndex,string path) {
	faIndex.file=openFile(path+".seq","r+");
	FILE *ind=openFile(path+".ind","r");
  	while (!feof(ind)) {
		FaRecord record=readIndexRecord(ind);
		if (record.id.size()>0) faIndex.record[record.id]=record;
	}
  	fclose(ind);
}

char* getFaIndexSeq(FaIndex& faIndex,string seqId) {
	FaRecord ind=faIndex.record[seqId];
	fseek(faIndex.file,ind.offset,0);
	char* seq=(char*)malloc(ind.length*sizeof(char));
	fread(seq,sizeof(char),ind.length,faIndex.file);
	return seq;
}

char* getMFaIndexSeq(FaIndex& faIndex,string seqId,int protoIndex) {
	char protoId[20];
	sprintf(protoId,"%d",protoIndex);
	string id=seqId+":"+protoId;
	FaRecord ind=faIndex.record[id];
	fseek(faIndex.file,ind.offset,0);
	char* seq=(char*)malloc(ind.length*sizeof(char));
	fread(seq,sizeof(char),ind.length,faIndex.file);
	return seq;
}
