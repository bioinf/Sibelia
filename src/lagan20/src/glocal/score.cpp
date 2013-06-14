#include<structs.h>
#include<score.h>
#include<leftinfluence.h>
#include<rightinfluence.h>
#include<fstream>

extern vector<class Score*> scoreFunctions[1<<(UPSTRANDBITS+DOWNSTRANDBITS+RELPOSBITS)];


float Score::getScore(Fragment *up, Fragment * down) {
	long long int absSeq1,absSeq2,absDiagonal,absMin,absMax;

	absSeq1= Myabs((up->seq1End) - (down->seq1Start));
	absSeq2= Myabs((up->seq2End) - (down->seq2Start));

	absMin = Mymin(absSeq1,absSeq2);
	absMax=Mymax(absSeq1,absSeq2);

	absDiagonal = absMax-absMin;

	return absMin*(-minConstant) + absMax* (-maxConstant) + absDiagonal *(-diagConstant) -openConstant +up->totalScore;
}


ScoreInterface::ScoreInterface (float iopenConstant, float iminConstant, float imaxConstant, float idiagConstant) {
	openConstant = iopenConstant;
	minConstant = iminConstant;
	maxConstant = imaxConstant;
	diagConstant = idiagConstant;
}


Score::Score (float iopenConstant , float iminConstant ,float imaxConstant,float idiagConstant):ScoreInterface(iopenConstant,iminConstant, imaxConstant, idiagConstant) {

}


void initScoreFunctionPointers(char * scoreFileName) {
	ifstream SFP;
	char line[255];

	SFP.open(scoreFileName);

	if (!SFP.good()) {
		printf("The score file is invalid");
		exit(0);
	}

	while (1) {
		SFP.getline(line,255);
		if (line[0]=='\0') { break; }
		createScoreFunctionObjects(line);
	}
}

void createScoreFunctionObjects(char * line) {
	long long int i;
	long long int j;
	long long int rem[4];
	long long int remCases[MAXCASES],remObjects[MAXOBJECTS];
	long long int numCases;
	long long int numObjects;
	long long int cases [MAXCASES];
	float objects[MAXOBJECTS][4];
	char updir,downdir,relpos;

	Score * SFObjects[MAXOBJECTS];

	j=0;

	for (i=0; (unsigned)i<strlen(line); i++) {
		if (line[i]=='{' || line[i]=='}') {
			rem[j++]=i;
		}
	}

	//forming cases

	numCases=0;

	for (i=rem[0]; i<=rem[1]; i++) {
		if (line[i]=='{' ||line[i]=='}'||line[i]==';') {
			remCases[numCases++]=i;
		}
	}

	numCases--;

	for (i=0; i<numCases; i++) {
		sscanf(&line[remCases[i]+1],"%c %c %c",&updir,&relpos,&downdir);
		if (DEBUG) { fprintf(stderr,"\n%c %c %c",updir,downdir,relpos); }
		cases[i]= charToCase(updir)<<UPSTRANDSHIFT | charToCase(downdir)<<DOWNSTRANDSHIFT |charToCase(relpos)<<RELPOSSHIFT;
	}

	numObjects=0;
	for (i=rem[2]; i<=rem[3]; i++) {
		if (line[i]=='{' || line[i]=='}' || line[i]==';') {
			remObjects[numObjects++]=i;
		}
	}

	numObjects--;

	for (i=0; i<numObjects; i++) {
		sscanf(&line[remObjects[i]+1],"%f %f %f %f",&objects[i][0],&objects[i][1],&objects[i][2],&objects[i][3]);
		if (DEBUG) { fprintf(stderr,"\t%f %f %f %f\n",objects[i][0],objects[i][1],objects[i][2],objects[i][3]); }
		SFObjects[i] = new Score(objects[i][0],objects[i][2],objects[i][3],objects[i][1]);
	}

	for (i=0; i<numCases; i++) {
		for (j=0; j<numObjects; j++) {
			scoreFunctions[cases[i]].push_back(SFObjects[j]);
		}
	}
}


long long int charToCase(char in) {
	switch(in) {
		case '+': return POSITIVE;
		case '-': return NEGATIVE;
		case 'R': return RIGHT;
		case 'L': return LEFT;
		case 'U': return UNRELATED;

		default:
		{
			fprintf(stderr,"\n Unrecognisable character in score file");
			exit(0);
		}
	}
}


float scoreAll(Fragment * up, Fragment * down, long long int ret_case) {
	unsigned long long int i;
//  TODO TODO TODO
	float ret_score=NEGINF;
//	float ret_score = -99999999999;
	float temp_score;

	if (up->nameIter != down->nameIter) {
		if (ret_case >> RELPOSSHIFT != UNRELATED) {
			//MUKCHECK HOPE THIS WORKS
			return NEGINF;
		}
	}

	for (i=0; i<scoreFunctions[ret_case].size(); i++) {
		temp_score = scoreFunctions[ret_case][i]->getScore(up,down);

		if (temp_score > ret_score) {
			ret_score = temp_score;
		}
	}

	if (ret_score == NEGINF) {
		printf("Score function case not handled::%lld\n",ret_case);
		//exit(0);
	}
	return ret_score;
}


long long int Mymax(long long int a, long long int b) {
	return (a>=b ? a : b);
}


long long int Mymin(long long int a,long long int b) {
	return (a<=b ? a : b);
}


long long int Myabs(long long int a) {
	return (a<0 ? -a : a);
}


float fragmentSetScore(Fragment * current, Fragment *owner, LI *LeftInfluence, RI * RightInfluence, long long int rightInfluenceFlag) {
	/*SLAGANCHANGE change call to the score based on the Leftinfluence, this has to be passed i guess*/
	float tempScore;

	if (rightInfluenceFlag == 3) {
		tempScore = scoreAll(owner,current, current->strand << DOWNSTRANDSHIFT | owner->strand <<UPSTRANDSHIFT | UNRELATED<< RELPOSSHIFT);
        if (tempScore == NEGINF) { // TODO
            if (current->totalScore <= 0) {
                current->totalScore = current->score;
                current->back = owner;
            }
        } else 
        if (tempScore + current->score > current->totalScore) {
            current->totalScore = tempScore + current->score;
			current->back = owner;
		}
	} else if (rightInfluenceFlag == TRUE) {
		tempScore = scoreAll(owner,current,RightInfluence->scoreIndex);

        if (tempScore == NEGINF) { // TODO
            if (current->totalScore <= 0) {
                current->totalScore = current->score;
                current->back = owner;
            }
        } else
        if (tempScore + current->score > current->totalScore) {
			current->totalScore = tempScore + current->score;
			current->back = owner;
		}
	} else {
		tempScore = scoreAll(owner,current,LeftInfluence->scoreIndex);
        
        if (tempScore == NEGINF) { // TODO
            if (current->totalScore <= 0) {
                current->totalScore = current->score;
                current->back = owner;
            }
        } else
        if (tempScore + current->score > current->totalScore) {
			current->totalScore = tempScore + current->score;
			current->back = owner;
		}
	}

	return current->totalScore;
}
