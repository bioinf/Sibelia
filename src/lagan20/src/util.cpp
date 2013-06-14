#define fastaRowLength 50
#define bufSize 2000

int trim(char* s) {
	int i=strlen(s);
	while (i>0 && (s[i-1]=='\n' || s[i-1]=='\r')) s[--i]='\0';
	return i;
}

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

char comp(char c) {
	switch(c) {
		case 'A': return 'T';
		case 'T': return 'A';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'N': return 'N';
		case 'a': return 't';
		case 't': return 'a';
		case 'c': return 'g';
		case 'g': return 'c';
		case 'n': return 'n';
		default: return c;
	}
}

char mask(char c) {
	return islower(c)?'N':c;
}
