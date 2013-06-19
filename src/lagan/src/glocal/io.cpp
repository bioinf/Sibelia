#include<structs.h>
#include<glocal.h>
#include<io.h>
#include<algorithm>

extern vector <Fragment> fragments;
extern vector <Point> startPoints;
extern vector <Point> endPoints;
extern Name allNames;

bool PointCompare(const Point &f1, const Point &f2) {
	if (f1.seq1 < f2.seq1) {
		return (f1.seq1 < f2.seq1);
	} else if (f1.seq1 == f2.seq1) {
		return (f1.seq2 < f2.seq2);
	} else {
		return (f1.seq1 < f2.seq1);
	}
}


//internal function that i dont need to care about.
char* rolltonum(char* str) {
	char *got1 = 0, *got2 = 0;
	long long int in = 0, i = 0;
	while (1) {
		if (str[i] == 0) { break; }

		if (str[i] == ';' && got1 && got2) { return got1; }

		if (isdigit(str[i])) {
			if (!in && (!i || isspace(str[i-1]))) {
				if (got1) {
					got2 = &str[i];
				} else {
					got1 = &str[i];
				}
				in = 1;
			}
		} else if (in && isspace(str[i])) {
			if (got2) {
				got1 = got2; got2 = 0; in = 0;
			}
			in = 0;
		} else {
			got1 = got2 = NULL;
		}
		i++;
	}
	return &str[i];
}


//reads one line of input at a time.
long long int getline(FILE *infile, hll *tt) {
	char temp[1024];
	char* help;
	long long int z;
	int h;
	fgets(temp, 1024, infile);
	sscanf(temp, "%s", tt->seq1Name);

	help = rolltonum(temp);
	z = sscanf(help, "%lld %lld;%n", &tt->seq1start, &tt->seq1end, &h);
	if (z < 2) { return 0; }

	sscanf(help+h, "%s", tt->seq2Name);
	help = rolltonum(help + h);

	if (sscanf(help, "%lld %lld; score = %f (%c)\n", &tt->seq2start, &tt->seq2end, &tt->score, &tt->strand)<3) {
		return 0;
	} else {
		return 1;
	}
}


void printFragment ( Fragment * curfrag ) {
	if (curfrag == NULL) {
		printf("done");
		return;
	}
    else if (curfrag->score == -1) {
        return;
    }

	// TODO: remove space after s2 and check supermap sorts and regexes
	printf("(%lld %lld)=(%lld %lld) %f %c [%f] s1:%s s2: %s\n",
		curfrag->seq1Start,
		curfrag->seq1End,
		curfrag->seq2Start-curfrag->base,
		curfrag->seq2End-curfrag->base,
		curfrag->score,
		(curfrag->strand==POSITIVE)?'+':'-',
		curfrag->totalScore,
		curfrag->seq1Name,
		curfrag->seq2Name
	);
}


void printAllFragments(long long int numFragments) {
	long long int i;
	for (i=0; i<numFragments; i++) {
		printFragment(&fragments[i]);
	}
	return;
}


// prints a chain upwards starting at the fragment called last.
long long int printChain(Fragment *current) {
	while (current) {
		printFragment(current);
		current = current->back;
	}
	return 0;
}


void swap(long long int *a, long long int *b) {
	long long int temp;
	temp = *a;
	*a = *b;
	*b = temp;
}


// initialises the parameters for a fragment.
// note the swap at the end of this function.
Fragment createFragment(hll *temp) {
	Fragment frag;
	frag.seq1Start = temp->seq1start;
	frag.seq1End = temp->seq1end;

	frag.seq2Start = temp->seq2start;

	frag.seq2End = temp->seq2end;

	strcpy(frag.seq1Name, temp->seq1Name);
	strcpy(frag.seq2Name, temp->seq2Name);

	if (temp->strand == '+') {
		frag.strand = POSITIVE;
	} else {
		frag.strand = NEGATIVE;
	}

	frag.score = temp->score;

	frag.back = NULL;

	frag.totalScore = -1;
	frag.deleted = FALSE;

	if (frag.seq1Start > frag.seq1End) {
		swap(&(frag.seq1Start), &(frag.seq1End));
	}
	return frag;
}


// reads the input file and returns the number of fragments read.
long long int readInput(char * fileName) {
	hll tempInput;
	FILE * fp;
	long long int i=0;
	char line[1024];

	unsigned long long int line_count = 0;

	fp = fopen(fileName, "r");

	if (!fp) {
		printf("SLAGAN: Error: Could not open file '%s'\n", fileName);
		exit(0);
	} else if (feof(fp)) {
		printf("SLAGAN: Error: Empty file %s\n", fileName);
		exit(0);
	}

	// Count the number of lines in the file
	while (fgets(line, 1023, fp)) {
		line_count++;
	}
	rewind(fp);

	fragments.reserve(line_count);

	while (!feof(fp)) {
		while (!feof(fp) && !getline(fp, &tempInput));
		if (feof(fp)) { break; }

		// ignoring the low scoring fragments ?
		if (tempInput.score < CUTOFF ) { continue; }

		//createfragment

		fragments.push_back(createFragment(&tempInput));
		i++;
	}

	return i;
}


void createPointLists(long long int numFragments) {
	long long int i;
	Point startPoint, endPoint;

	//SLAGANCHANGE:: Push -seq2,seq1 on the start list as well.

	for (i=0; i<numFragments; i++) {
		startPoint.seq1 = fragments[i].seq1Start;
		startPoint.seq2 = fragments[i].seq2Start;
		endPoint.seq1 = fragments[i].seq1End;
		endPoint.seq2 = fragments[i].seq2End;
		startPoint.frag = &fragments[i];
		endPoint.frag = &fragments[i];
		startPoints.push_back(startPoint);

		startPoint.seq2 = -fragments[i].seq2Start;
		startPoints.push_back(startPoint);
		endPoints.push_back(endPoint);
	}
	sort(startPoints.begin(), startPoints.end(), PointCompare);
	sort(endPoints.begin(), endPoints.end(), PointCompare);
}


void printPointLists(long long int numFragments) {
	long long int i;
	printf("StartPoint lists:\n");

	for (i=0; i<numFragments; i++) {
		printf(" seq1 :%lld seq2:%lld \n", startPoints[i].seq1, startPoints[i].seq2);
	}

	printf("EndPoint lists:\n");
	for (i=0; i<numFragments; i++) {
		printf(" seq1 :%lld seq2:%lld \n", endPoints[i].seq1, endPoints[i].seq2);
	}
	printf("End lists");
}


void findAllNames(long long int numFragments) {
	long long int i;
	long long int size;
	long long int numContigs=0;
	Name::iterator currName;

	for (i=0; i<numFragments; i++) {
		size = fragments[i].seq2Start>fragments[i].seq2End ? fragments[i].seq2Start : fragments[i].seq2End;

		currName = allNames.find(fragments[i].seq2Name);

		if (currName != allNames.end()) {
			if (currName->second < size) {
				currName->second = size;
			}
		} else {
			allNames[fragments[i].seq2Name] = size;
			numContigs ++;
		}
	}
	if (DEBUG) { fprintf(stderr, "The number of contigs is %lld",numContigs); }
}


void decideContigBase() {
	Name::iterator currName;
	long long int offset =0;
	long long int temp;

	for (currName=allNames.begin(); currName!=allNames.end(); currName++) {
		temp = currName->second;
		currName->second = offset;
		offset += (10 + temp);
	}
}


void storeIterators(long long int numFragments) {
	long long int i;

	for (i=0; i<numFragments; i++) {
		fragments[i].nameIter = allNames.find(fragments[i].seq2Name);
		fragments[i].seq2Start += (fragments[i].nameIter)->second;
		fragments[i].seq2End += (fragments[i].nameIter)->second;
		fragments[i].base = (fragments[i].nameIter)->second;
	}
}
