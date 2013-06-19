#include<score.h>
#include<glocal.h>
#include<algorithm>

bool seq1StartCompare(const Fragment &f1, const Fragment &f2) {
	return f1.seq1Start < f2.seq1Start;
}

//vectors that would be needed globally
vector<Fragment> fragments;
vector<Point>startPoints;
vector<Point>endPoints;
long long int numFragments;
InterPoint inter;


/*SLAGANCHANGE This has to change*/

RI RI_regions[1<<(UPSTRANDBITS+DOWNSTRANDBITS+RELPOSBITS)];
LI LI_regions[1<<(UPSTRANDBITS+DOWNSTRANDBITS+RELPOSBITS)];

vector<class Score*> scoreFunctions[1<<(UPSTRANDBITS+DOWNSTRANDBITS+RELPOSBITS)];

Name allNames;


extern Fragment LI_dummy;
Fragment * unrelatedFrag;

Fragment *max_score_index;
float max_score;

int main(int, char **argv) {
	long long int nextEndRow,nextStartRow, nextInterPointRow;
	long long int i;
	Point intersectionPoint;

	numFragments = readInput(argv[1]);

	findAllNames( numFragments);
	decideContigBase();
	storeIterators(numFragments);

	initScoreFunctionPointers(argv[2]);
	unrelatedFrag = &LI_dummy;

	/*SLAGANCHANGE  need a LI, RI pointer array and init */
	/*SLAGANCHANGE:: Need score function init */

	if (DEBUG) { fprintf(stderr,"Numfrg::%lld",numFragments); }
	max_score_index=NULL;
	max_score =-INF;

	long long int break_flag =0;

	createPointLists(numFragments);
//	printFragmentsInPointListOrder(numFragments);
//	exit(0);

	//The initial Row upto which startPointHandler goes
	nextEndRow = endPoints[0].seq1;
	nextStartRow = startPoints[0].seq1;

	for (i=0;i<1<<TOTALSHIFT;i++) {
		initRI(&RI_regions[i],i);
		InitLI(&LI_regions[i],i);
	}

	if (DEBUG) { fprintf(stderr,"The number of regions was %lld",i); }

	while (1) {
		if (inter.begin()==inter.end()) {
			nextInterPointRow = INF;
			if (DEBUG) { fprintf(stderr,"\nORHERE"); }
		} else {
			intersectionPoint = (inter.begin())->first;
			nextInterPointRow = intersectionPoint.seq1;
			if (DEBUG) { fprintf(stderr,"\nHERE"); }
		}

		if (nextStartRow <= nextEndRow) {
			//CHANGE HERE
			if (nextStartRow<nextInterPointRow) {
				nextStartRow=startPointHandler();

				if (nextStartRow == INF) {
					//break;
					break_flag = 1;
				}
			} else {
				intersectionPointHandler();
			}
		} else {
			//CHANGE HERE
			if (nextEndRow<nextInterPointRow) {
				nextEndRow=endPointHandler();
				if (break_flag == 1) {
					break;
				}
			} else {
				intersectionPointHandler();
			}
		}
	}

	if (DEBUG) { fprintf(stderr,"\nMAX CHAIN\n"); }
	printChain(max_score_index);

	//fprintf(stderr,"\nALL\n");
	//printAllFragments(numFragments);
	return 0;
}


//Processes till the row number reaches the argument
long long int startPointHandler() {
	static long long int current=0;
	Fragment *owner;
	long long int current_seq1= startPoints[current].seq1;
	float current_score;
	if (DEBUG) { fprintf(stderr,"\nStart PointHandler"); }

	while (startPoints[current].seq1==current_seq1) {
		long long int upStrand,downStrand,relPos,possibleCase;

		downStrand = (startPoints[current].frag)->strand;

		relPos = startPoints[current].seq2 > 0 ? RIGHT:LEFT;

		upStrand = POSITIVE;
		possibleCase = downStrand << DOWNSTRANDSHIFT | upStrand <<UPSTRANDSHIFT | relPos<< RELPOSSHIFT;

		owner=LILookUpOwnerStart(&LI_regions[possibleCase],startPoints[current].frag);

		current_score = fragmentSetScore(startPoints[current].frag, owner, &LI_regions[possibleCase], NULL, FALSE);

		owner = lookUpOwnerStart(&RI_regions[possibleCase], startPoints[current].frag);

		current_score = fragmentSetScore(startPoints[current].frag, owner, NULL, &RI_regions[possibleCase], TRUE);

		upStrand = NEGATIVE;
		possibleCase = downStrand << DOWNSTRANDSHIFT | upStrand <<UPSTRANDSHIFT | relPos << RELPOSSHIFT;

		owner = lookUpOwnerStart(&RI_regions[possibleCase], startPoints[current].frag);

		current_score = fragmentSetScore(startPoints[current].frag, owner, NULL,&RI_regions[possibleCase], TRUE);
		if (DEBUG) { fprintf(stderr, "HI1"); }

		owner = LILookUpOwnerStart(&LI_regions[possibleCase],startPoints[current].frag);
		current_score = fragmentSetScore(startPoints[current].frag, owner, &LI_regions[possibleCase], NULL, FALSE);
		if (DEBUG) { fprintf(stderr, "HI2"); }

		current_score = fragmentSetScore(startPoints[current].frag, unrelatedFrag, NULL, NULL, 3);
		if (DEBUG) { fprintf(stderr, "HI3"); }

		if ((startPoints[current].frag)->back == NULL) {
			if (DEBUG) { fprintf(stderr, "\n The fragment did not chain!"); }
			// exit(1);
		} else if (DEBUG) {
			fprintf(stderr, "Score for the current fragment is::%f", startPoints[current].frag->totalScore);
			fprintf(stderr, "Score for the owner fragment is::%f", startPoints[current].frag->back->totalScore);
		}

		if (startPoints[current].frag->totalScore > max_score) {
			max_score = startPoints[current].frag->totalScore;
			max_score_index = startPoints[current].frag ;
		}

		current++;

		if (DEBUG) { fprintf(stderr,"\ncurrent fragment is %lld",current); }
		
		if (current>=2*numFragments) {
			return INF;
		}
	}

	return startPoints[current].seq1;
}


//takes as arguements the start row number and the end row number and processes all the rows
//This would usually have to find the case
long long int endPointHandler() {
	static long long int current=0;

	long long int current_seq1= endPoints[current].seq1;

	if (DEBUG) { fprintf(stderr,"\nEnd PointHandler"); }

	/*SLAGANCHANGE:: There is going to be a commit to 4 strucures depending on the strand, loop with continue*/
	/*SLAGANCHANGE:: find the best scoring fragment in the current row and update the best so far at the end*/

	while (endPoints[current].seq1 == current_seq1) {
		long long int upStrand, downStrand, relPos, possibleCase;

		//MUKFIXME: This sends the highest scoring one into the leftinfluence machinery

		while (current<2*numFragments-1 &&( endPoints[current].seq1== endPoints[current+1].seq1) && (endPoints[current+1].seq2 == endPoints[current].seq2)) {
			if ((endPoints[current].frag->totalScore) > (endPoints[current+1].frag->totalScore)) {
				Fragment * temp;

				temp=endPoints[current+1].frag;
				endPoints[current+1].frag=endPoints[current].frag;
				endPoints[current].frag =temp;
			}
			current++;
		}

		/*
		if( current>1 &&(endPoints[current].seq1== endPoints[current-1].seq1) && (endPoints[current-1].seq2 == endPoints[current].seq2))
		{
		current++;
		continue;
		}
		*/
		upStrand = endPoints[current].frag->strand;

		// This works because POSITIVE and NEGATIVE are 0 and 1
		// This works because LEFT and RIGHT are 0 and 1

		for (downStrand=0;downStrand<2;downStrand++) {
			for (relPos=0;relPos<2;relPos++) {
				possibleCase = downStrand << DOWNSTRANDSHIFT | upStrand <<UPSTRANDSHIFT | relPos<< RELPOSSHIFT;

				RICommitEndPoint(&RI_regions[possibleCase],endPoints[current].frag);
				LICommitPoint(&LI_regions[possibleCase],endPoints[current].frag);
			}
		}

		if (endPoints[current].frag->totalScore > unrelatedFrag->totalScore)
		unrelatedFrag = endPoints[current].frag;

		current++;
	}

	return endPoints[current].seq1;
}


void intersectionPointHandler() {
	long long int current_seq1;
	Point p,curr;

	p=inter.begin()->first;

	current_seq1=p.seq1;

	if (DEBUG) { fprintf(stderr,"\nIntersection PointHandler"); }
	do {
		// printState(&LI_regions[0]);
		HandleOneIntersectionPoint();

		//printState(&LI_regions[0]);
		p=inter.begin()->first;
		current_seq1=p.seq1;
	} while (current_seq1==curr.seq1);
}
