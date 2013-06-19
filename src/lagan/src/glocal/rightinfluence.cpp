#include <rightinfluence.h>

Fragment origin, end;

// Sets the first default owner of the whole region
void initRI(RI *RightInfluence, long long int scoreIndex) {
	RightInfluence->scoreIndex = scoreIndex;

	if (((scoreIndex >> RELPOSSHIFT) & 1) == LEFT) {
		RightInfluence->reflectFlag = TRUE;
	} else {
		RightInfluence->reflectFlag = FALSE;
	}

	// will lose to anyone
	origin.seq1End = 0; origin.seq2End = 0;
    origin.seq1Start = 0; origin.seq2Start = 0;

	// hack to aid winner selection
	origin.score = -1;
	end.score = -2;
	origin.totalScore = end.totalScore = 0;

	// will win against anyone
	end.seq1End = 0; end.seq2End = 0;
	end.seq1Start = 0; end.seq2Start = 0;

	origin.back = NULL;

    RightInfluence->act[-INF] = &origin;
    RightInfluence->act[+INF] = &end;
}


// Finds the owner in the current right influence region and returns the score using the appropriate score function
float lookUpScore(RI * RightInfluence, Fragment * current) {
	Fragment* owner;

	// find the owner of the region that you are in
	owner = lookUpOwnerStart(RightInfluence, current);

	// return the score using the appropriate score function
	return scoreAll(owner, current, RightInfluence->scoreIndex);
}


// Returns the owner of the region
Fragment * lookUpOwnerStart(RI * RightInfluence, Fragment * current) {
	Active::iterator ownerIterator;

	// find the owner of the region that you are in.
	ownerIterator = RightInfluence->act.upper_bound(current->getSeq2Start(RightInfluence->reflectFlag) - current->seq1Start);
	ownerIterator--;

	return (*ownerIterator).second;
}


Fragment * lookUpOwnerEnd(RI * RightInfluence, Fragment * current) {
	Active::iterator ownerIterator;

	// find the owner of the region that you are in.
	ownerIterator=RightInfluence->act.upper_bound(current->getSeq2End(RightInfluence->reflectFlag) - current->seq1End);
	ownerIterator--;

	return (*ownerIterator).second;
}


// Returns true if the first argument is the winner in their common region
long long int RIWinner(RI * RightInfluence, Fragment * first, Fragment * second) {
	Fragment dummy;

	//if the first frag is the origin or the second frag is the end then the first frag loses
	if (first->score==-1 || second->score==-2) { return FALSE; }

	//if the first frag is the end or the second frag is the origin then the first frag wins
	if (second->score==-1 || first->score==-2) { return TRUE; }

	dummy.seq1Start = Mymax(first->seq1End, second->seq1End) + 1;
	dummy.seq2Start = Mymax(first->getSeq2End(RightInfluence->reflectFlag), second->getSeq2End(RightInfluence->reflectFlag)) + 2;

	if (first->getSeq2End(RightInfluence->reflectFlag) > second->getSeq2End(RightInfluence->reflectFlag)) {
		dummy.nameIter = first->nameIter;
	} else {
		dummy.nameIter = second->nameIter;
	}

	if (scoreAll(first, &dummy, RightInfluence->scoreIndex) > scoreAll(second, &dummy, RightInfluence->scoreIndex)) {
		return TRUE;
	} else {
		return FALSE;
	}
}


long long int RICommitEndPoint(RI * RightInfluence, Fragment * current) {
	Fragment * owner;
	Fragment * temp;
	owner = lookUpOwnerEnd(RightInfluence, current);

	if (RIWinner(RightInfluence, owner, current)) { return 0; }
    
	owner = nextOnActive(RightInfluence, owner);
    
	while (1) {
		if (RIWinner(RightInfluence, current, owner)) {
			temp = owner;
			owner = nextOnActive(RightInfluence, owner);
			RightInfluence->act.erase(temp->getSeq2End(RightInfluence->reflectFlag)-temp->seq1End);
		} else {
			break;
		}
	}

    //inserting into the list of active owners
	RightInfluence->act[current->getSeq2End(RightInfluence->reflectFlag) - current->seq1End] = current;

int possibleCase = NEGATIVE << DOWNSTRANDSHIFT | NEGATIVE <<UPSTRANDSHIFT | LEFT << RELPOSSHIFT;
if (RightInfluence->scoreIndex == possibleCase) {
    Active::iterator j,i = RightInfluence->act.begin();
    i++;
    while(i != RightInfluence->act.end()) {
        //    if (i == NULL) { continue;}
        j = i;
        j++;
        if (j != RightInfluence->act.end()) {
            if ((*j).second->score == -2) { break;} // j is act.end (why does the check above fail?)
            if ((*i).second->totalScore > (*j).second->totalScore) {
                /*                fprintf(stdout,"Assertion failed in RICommitEndPoint: Cur frag:\n");
                printFragment(current);
                fprintf(stdout,"Cur orig owner:\n");
                printFragment(tempOwner);
                fprintf(stdout,"Cur frag diag: %lld\n", (current->getSeq2End(RightInfluence->reflectFlag) - current->seq1End));
                fprintf(stdout,"    Frag 1 in pair (j):\n    ");
                printFragment((*j).second);
                fprintf(stdout,"    Frag 2 in pair (i):\n    ");
                printFragment((*i).second);
                fprintf(stdout,"RI:\n");
                printActive(RightInfluence);
                assert (0);
                */
                break;
                //            assert(i->first->score >= j->first->score);
            }
         }
        i++;
    }
}
 

    return 1;
}


long long int diagonal(Fragment * current, RI * RightInfluence) {
	return (current->getSeq2End(RightInfluence->reflectFlag) - current->seq1End);
}


// Returns the successor on the active list
Fragment * nextOnActive(RI * RightInfluence, Fragment * current) {
	Active::iterator holder;
	long long int diagCurrent;

	diagCurrent = current->getSeq2End(RightInfluence->reflectFlag) - current->seq1End;

    //MUKMOD start
    if(current->score==-1)
        {
            diagCurrent = -INF;
            
        }

    if(current->score ==-2)
        {
            diagCurrent = INF;
        }
    //MUKMOD end
        

	holder = RightInfluence->act.upper_bound(diagCurrent);

	if (holder != RightInfluence->act.end()) {
		return (*holder).second;
	} else {
		return NULL;
	}
}


long long int printActive(RI * RightInfluence) {
	Active::iterator temp;
	long long int i = 0;
    fprintf(stdout, "Active RI:\n");
	for (temp = RightInfluence->act.begin(); temp != RightInfluence->act.end(); temp++) {
		fprintf(stdout, " %lld", (*temp).first);
        fprintf(stdout, ":sc=%f:totsc=%f;",((*temp).second)->score, ((*temp).second)->totalScore);
		i++;
	}
    fprintf(stdout, "\n");
	return i;
}
