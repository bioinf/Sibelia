#include<leftinfluence.h>

Fragment LI_dummy;

// Returns the fragment who is the owner of the region in which the current point is
Owner::iterator LILookUpOwnerIterator(LI * LeftInfluence, long long int seq1, long long int seq2) {
	CBound::iterator citer;
	DBound::iterator diter;

	citer = (LeftInfluence->c).lower_bound(seq2);

	if ((LeftInfluence->c).end() == (LeftInfluence->c).begin() || (citer == (LeftInfluence->c).begin())) {
		return (LeftInfluence->o).end();
	}

	citer--;

	diter = (LeftInfluence->d).upper_bound(seq2 - seq1);

	if (diter == (LeftInfluence->d).begin()) {
		return citer->second;
	}

	diter--;

	if ((citer->first - diter->first) > seq1) {
		return citer->second;
	} else {
		return diter->second;
	}
}


Fragment * LILookUpOwnerEnd(LI * LeftInfluence,Fragment * current) {
	Owner::iterator own = LILookUpOwnerIterator(LeftInfluence, current->seq1End, current->getSeq2End(LeftInfluence->reflectFlag));

	if (own == (LeftInfluence->o).end()) {
		return &LI_dummy;
	} else {
		return *own;
	}
}


Fragment * LILookUpOwnerStart(LI * LeftInfluence, Fragment * current) {
	Owner::iterator own = LILookUpOwnerIterator(LeftInfluence, current->seq1Start, current->getSeq2Start(LeftInfluence->reflectFlag));

	if (own == (LeftInfluence->o).end()) {
		return &LI_dummy;
	} else {
		return *own;
	}
}


// Returns the column boundary before the current point, if there is none it returns end
CBound::iterator LICColumn(LI * LeftInfluence, long long int /* seq1 */, long long int seq2) {
	CBound::iterator citer;

	citer = (LeftInfluence->c).lower_bound(seq2);

	//should not decrement, also means that the point is before all the column boundaries.
	//FIX #2 if(citer == (LeftInfluence->c).begin())

	if ((LeftInfluence->c).end() == (LeftInfluence->c).begin() || (citer == (LeftInfluence->c).begin())) {
		return (LeftInfluence->c).end();
	} else {
        citer--;
        return citer;
    }
}


Fragment * LICOwner(LI * LeftInfluence, long long int seq1, long long int seq2) {
	CBound::iterator citer;
	citer = LICColumn(LeftInfluence, seq1, seq2);

	if (citer == (LeftInfluence->c).end()) {
		return &LI_dummy;
	} else {
		return *(citer->second);
	}
}


Fragment * LIDOwner(LI * LeftInfluence, long long int seq1, long long int seq2) {
	DBound::iterator diter;
	diter = LIDDiagonal(LeftInfluence, seq1, seq2);

	if (diter == (LeftInfluence->d).end()) {
		return &LI_dummy;
	} else {
		return *(diter->second);
	}
}


//returns the diagonal boundary,  or end if all the point is before all the diagonal boundaries
DBound::iterator LIDDiagonal(LI * LeftInfluence, long long int seq1, long long int seq2) {
	DBound::iterator diter;

	diter = (LeftInfluence->d).upper_bound(seq2-seq1);

	if ((LeftInfluence->d).end() == (LeftInfluence->d).begin() || diter == (LeftInfluence->d).begin()) {
		return (LeftInfluence->d).end();
	} else {
        diter--;
        return diter;
    }
}


// this function should never get called with the LI dummy
// can the scores become negative and how do we handle this?
float LILookUpScore(LI * LeftInfluence, Fragment * current) {
	Fragment * owner = LILookUpOwnerStart(LeftInfluence, current);

	if (owner==NULL) {
		fprintf(stderr,"Owner NULL in call LILookUpScore");
		exit(0);
	}

	if (owner->score == -1) {
		//MUKCHECK
		return -1;
	} else {
		return scoreAll(owner,current,LeftInfluence->scoreIndex);
	}
}


void InitLI(LI * LeftInfluence, long long int scoreIndex) {
	LeftInfluence->scoreIndex = scoreIndex;

	if (((scoreIndex >> RELPOSSHIFT) & 1) == LEFT) {
		LeftInfluence->reflectFlag = TRUE;
	} else {
		LeftInfluence->reflectFlag = FALSE;
	}

	LI_dummy.score = -1;
	LI_dummy.totalScore = 0;
	LI_dummy.back = NULL;

	//there will be a list of structures to insert this into
	(LeftInfluence->o).insert((LeftInfluence->o).begin(), &LI_dummy);
}


long long int LI_Winner(LI * LeftInfluence, Fragment * first, Fragment * second) {
	Fragment dummy;

	if (first->score == -1) { return FALSE; }

	if (second->score == -1) { return TRUE; }

	dummy.seq1Start = max(first->seq1End, second->seq1End) + 2;
	dummy.seq2Start = max(first->getSeq2End(LeftInfluence->reflectFlag), second->getSeq2End(LeftInfluence->reflectFlag)) + 1;

	if (first->getSeq2End(LeftInfluence->reflectFlag) > second->getSeq2End(LeftInfluence->reflectFlag)) {
		dummy.nameIter = first->nameIter;
	} else {
		dummy.nameIter = second->nameIter;
	}

	if (scoreAll(first, &dummy, LeftInfluence->scoreIndex) >= scoreAll(second, &dummy, LeftInfluence->scoreIndex)) {
		return TRUE;
	} else {
		return FALSE;
	}
}


long long int LICommitPoint(LI * LeftInfluence, Fragment * current) {
	Owner::iterator cowner, ownerIter;
	Fragment * owner;
	CBound::iterator citer;
	DBound::iterator diter;
	long long int colFlag;

	ownerIter = LILookUpOwnerIterator(LeftInfluence, current->seq1End, current->getSeq2End(LeftInfluence->reflectFlag));

	citer = LICColumn(LeftInfluence, current->seq1End, current->getSeq2End(LeftInfluence->reflectFlag));
	diter = LIDDiagonal(LeftInfluence, current->seq1End, current->getSeq2End(LeftInfluence->reflectFlag));
	owner = LILookUpOwnerEnd(LeftInfluence, current);

	if (citer == (LeftInfluence->c).end()) {
		colFlag = TRUE;
	} else if (diter == (LeftInfluence->d).end()) {
		colFlag = TRUE;
	} else {
		cowner = citer->second;
		if (cowner == ownerIter) {
			colFlag = TRUE;
		} else {
			colFlag = FALSE;
		}
	}

	if (LI_Winner(LeftInfluence, owner, current)) {
		return FALSE;
	}

	if (colFlag) {
		return LI_CommitColumnOwner(LeftInfluence, current, owner);
	} else {
		return LI_CommitDiagonalOwner(LeftInfluence, current, owner);
	}
}


Owner::iterator LI_OwnerInsertAfter(LI * LeftInfluence, Owner::iterator current, Fragment * curfrag) {
	current++;
	return (LeftInfluence->o).insert(current, curfrag);
}


long long int LI_CommitDiagonalOwner(LI * LeftInfluence, Fragment * current, Fragment * owner) {
	CBound::iterator current_column, next_column;
	DBound::iterator current_diagonal, prevDiag;
	DInter::iterator current_diag_inter, my_diag_inter, prevDiagInter;
	CInter::iterator my_col_inter, next_column_inter, colInter;

	Owner::iterator own, tempowner;

	//searching for the next column to switch on
	current_column = LICColumn(LeftInfluence, current->seq1End, current->getSeq2End(LeftInfluence->reflectFlag));
	current_diagonal = LIDDiagonal(LeftInfluence, current->seq1End, current->getSeq2End(LeftInfluence->reflectFlag));
	current_diag_inter = (LeftInfluence->di).find(current_diagonal->first);
	own = LILookUpOwnerIterator(LeftInfluence, current->seq1End, current->getSeq2End(LeftInfluence->reflectFlag));

	//this implies that the point is before all the cbounds:: THIS CANT HAPPEN!!

	if (current_column == (LeftInfluence->c).end()) {
		//FIX#7
		fprintf(stderr, "\n diagonal owner, but no column before it");
		exit(0);
	} else {
		next_column = current_column;
		next_column++;
	}

	//2cases
	if (next_column == (LeftInfluence->c).end() || next_column->first > current->getSeq2End(LeftInfluence->reflectFlag)) {
		if (current_diagonal->first < current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End) {
			if (DEBUG) { fprintf(stderr, "In Diagonal Commit::FIRSTCASE"); }

			tempowner = LI_OwnerInsertAfter(LeftInfluence, current_diagonal->second, current);
			(LeftInfluence->c)[current->getSeq2End(LeftInfluence->reflectFlag)] = tempowner;
			(LeftInfluence->ci)[current->getSeq2End(LeftInfluence->reflectFlag)] = inter.end();
			my_col_inter = (LeftInfluence->ci).find(current->getSeq2End(LeftInfluence->reflectFlag));

			tempowner = LI_OwnerInsertAfter(LeftInfluence, tempowner, owner);

			(LeftInfluence->d)[current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End] = tempowner;
			(LeftInfluence->di)[current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End] = inter.end();
			my_diag_inter = (LeftInfluence->di).find(current->getSeq2End(LeftInfluence->reflectFlag)-current->seq1End);

			if (next_column!= (LeftInfluence->c).end()) {
				next_column_inter = (LeftInfluence->ci).find(next_column->first);

				if (next_column_inter->second == current_diag_inter->second && current_diag_inter->second!=inter.end()) {
					DeleteIntersectionPoint(next_column_inter->second, next_column_inter, current_diag_inter);
					CreateIntersectionPoint(LeftInfluence, next_column->first,
                                            current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End,
                                            next_column_inter, my_diag_inter);
				} else if (next_column_inter->second == inter.end()) {
					CreateIntersectionPoint(LeftInfluence, next_column->first,
                                            current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End,
                                            next_column_inter, my_diag_inter);
				}
			}

			CreateIntersectionPoint(LeftInfluence, current->getSeq2End(LeftInfluence->reflectFlag),
                                    current_diagonal->first, my_col_inter, current_diag_inter);
		} else {
			if (DEBUG) { fprintf(stderr, "\n In Diagonal Commit:SECONDCASE"); }

			//There will be a previous owner as this is a diagonal case
			own = LILookUpOwnerIterator(LeftInfluence, current->seq1End, current->getSeq2End(LeftInfluence->reflectFlag));
			own--;

			if (LI_Winner(LeftInfluence, *own, current)) {
				return FALSE;
			}

			own++;
			tempowner = (LeftInfluence->o).insert(own, current);
			(LeftInfluence->c)[current->getSeq2End(LeftInfluence->reflectFlag)] = tempowner;
			(LeftInfluence->ci)[current->getSeq2End(LeftInfluence->reflectFlag)] = inter.end();
			colInter = (LeftInfluence->ci).find(current->getSeq2End(LeftInfluence->reflectFlag));

			//There is no diagonal here

			//intersection Point Handling
			// check is the previous intersection Point exists, if it does check if the flag is off in which
			//case insert an intersection Point into Intersect and Handle flags appropriately

			//There is a problem here
			//FIX #7 #4 major fix
			if (current_diagonal != (LeftInfluence->d).begin()) {
				prevDiag = current_diagonal;
				prevDiag--;

				prevDiagInter = (LeftInfluence->di).find(prevDiag->first);
				if (prevDiagInter->second == inter.end()) {
					CreateIntersectionPoint(LeftInfluence, current->getSeq2End(LeftInfluence->reflectFlag),
                                            prevDiag->first, colInter, prevDiagInter);
				}
			}
		}
	} else {
		if (DEBUG) { fprintf(stderr, "\n In Diagonal Commit:THIRDCASE"); }
		if (LI_Winner(LeftInfluence, *(next_column->second), current)) { return false; }

		tempowner = (LeftInfluence->o).insert(next_column->second, current);

		//He does the intersection point processing with lower priority!!?
		//This might mean that the diagonal entry already exists, also this might mean that
		//The intersection point processing removes the entry?!

		(LeftInfluence->d)[current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End] = next_column->second;
		(LeftInfluence->di)[current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End] = inter.end();
		my_diag_inter = (LeftInfluence->di).find(current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End);

		next_column->second = tempowner;

		//checking if the next column exists
		next_column++;

		if (next_column!= (LeftInfluence->c).end()) {
			next_column_inter =(LeftInfluence->ci).find(next_column->first);

			if (next_column_inter->second == inter.end()) {
				CreateIntersectionPoint(LeftInfluence, next_column->first,
                                        current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End,
                                        next_column_inter, my_diag_inter);
			}
		}
	}
	return TRUE;
}


long long int LI_CommitColumnOwner(LI * LeftInfluence, Fragment * current, Fragment * owner) {
	CBound::iterator current_column, next_column;
	CInter::iterator nextColInter, colInter;
	DInter::iterator diagInter;
	Owner::iterator tempowner;

	current_column= LICColumn(LeftInfluence, current->seq1End, current->getSeq2End(LeftInfluence->reflectFlag));

	if ((LeftInfluence->c).end() == (LeftInfluence->c).begin()) {
		//Init has already put in one fragment
		tempowner = LI_OwnerInsertAfter(LeftInfluence, (LeftInfluence->o).begin(), current);
		(LeftInfluence->c)[current->getSeq2End(LeftInfluence->reflectFlag)] = tempowner;
		(LeftInfluence->ci)[current->getSeq2End(LeftInfluence->reflectFlag)] = inter.end();

		//FIX #5 FIRST MAJOR FIX
		tempowner = LI_OwnerInsertAfter(LeftInfluence, tempowner, &LI_dummy);
		(LeftInfluence->d)[current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End] = tempowner;
		(LeftInfluence->di)[current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End] = inter.end();
		return TRUE;
	}

	// If the current_column is the end , that means that we are before all the column boundaries
	//as the other case has been taken care of above

	if (current_column == (LeftInfluence->c).end()) {
		next_column = (LeftInfluence->c).begin();
	} else {
		next_column = current_column;
		next_column++;
	}

	// Either the case that the column boundary is that last column boundary  or that the next column is after the current point

	if (next_column == (LeftInfluence->c).end() || next_column->first > current->getSeq2End(LeftInfluence->reflectFlag)) {
		if (DEBUG) { fprintf(stderr, "\nColCommit::FIRSTCASE"); }
		// this means that the next column is not the first column
		if (current_column != (LeftInfluence->c).end()) {
			tempowner = LI_OwnerInsertAfter(LeftInfluence, current_column->second, current);
		} else {
			// this means that the next column is the first column
			tempowner = LI_OwnerInsertAfter(LeftInfluence, (LeftInfluence->o).begin(), current);
		}

		(LeftInfluence->c)[current->getSeq2End(LeftInfluence->reflectFlag)] = tempowner;
		(LeftInfluence->ci)[current->getSeq2End(LeftInfluence->reflectFlag)] = inter.end();
		//This is inefficient
		colInter = (LeftInfluence->ci).find(current->getSeq2End(LeftInfluence->reflectFlag));
		tempowner = LI_OwnerInsertAfter(LeftInfluence, tempowner, owner);
		(LeftInfluence->d)[current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End] = tempowner;
		(LeftInfluence->di)[current->getSeq2End(LeftInfluence->reflectFlag)-current->seq1End] = inter.end();

		//This is inefficient
		diagInter = (LeftInfluence->di).find(current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End);

		//if there is a next column then there is an issue of an intersection point
		if (next_column != (LeftInfluence->c).end()) {
			nextColInter = (LeftInfluence->ci).find(next_column->first);

			if (nextColInter->second == inter.end()) {
				CreateIntersectionPoint(LeftInfluence, next_column->first,
                                        current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End, nextColInter, diagInter);
			}
		}
	} else {
		if (DEBUG) { fprintf(stderr, "\nColCommit::SECONDCASE"); }

		if (LI_Winner(LeftInfluence, *(next_column->second), current)) {
			return FALSE;
		}

		tempowner = (LeftInfluence->o).insert(next_column->second, current);
		(LeftInfluence->d)[current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End] = next_column->second;
		//FIX #6 SECOND MAJOR FIX
		(LeftInfluence->di)[current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End] = inter.end();

		//I dont think that i need this
		diagInter = (LeftInfluence->di).find(current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End);
		colInter = (LeftInfluence->ci).find(current->getSeq2End(LeftInfluence->reflectFlag));
		next_column->second = tempowner;

		//intersection Point handling
		next_column++;
		if (next_column != (LeftInfluence->c).end()) {
			nextColInter = (LeftInfluence->ci).find(next_column->first);

			if (nextColInter->second == inter.end()) {
				CreateIntersectionPoint(LeftInfluence, next_column->first,
                                        current->getSeq2End(LeftInfluence->reflectFlag) - current->seq1End, nextColInter, diagInter);
			}
		}
	}
	return TRUE;
}


void CreateIntersectionPoint(LI * LeftInfluence, long long int col, long long int diag, CInter::iterator colInter, DInter::iterator diagInter) {
	Point temp;

	InterPoint::iterator tempinter;
	temp.seq1 = col - diag;
	temp.seq2 = col;

	pair<Point,LI*> pairp(temp, LeftInfluence);
	tempinter = inter.insert(pairp);

	colInter->second = tempinter;
	diagInter->second = tempinter;
}


void DeleteIntersectionPoint(InterPoint::iterator tobeerased, CInter::iterator colInter, DInter::iterator diagInter) {
	inter.erase(tobeerased);
	colInter->second = inter.end();
	diagInter->second = inter.end();
}


// handles one intersection point that is at the head of inter
void HandleOneIntersectionPoint() {
	InterPoint::iterator head;
	Owner::iterator delOwner, leftOwner, rightOwner;

	CBound::iterator col, nextCol;
	CInter::iterator nextColInter, colInter;
	DInter::iterator prevDiagInter, diagInter;
	DBound::iterator diag, prevDiag;

	head = inter.begin();

	LI * LeftInfluence;

	//find the three owners that are invloved.
	LeftInfluence = head->second;

	col = (LeftInfluence->c).find((head->first).seq2);

	if (col == (LeftInfluence->c).end()) {
		fprintf(stderr, "\nIn HandleOneIntersectionPoint::The column does not exist. Point is %lld %lld", (head->first).seq1, (head->first).seq2);
		exit(0);
	}

	colInter = (LeftInfluence->ci).find(col->first);
	diag = (LeftInfluence->d).find((head->first).seq2 - (head->first).seq1);

	if (DEBUG) { fprintf(stderr, "\nIn HandleOneIntersectionPoint::The intersection point that is being handled: %lld %lld", (head->first).seq1, (head->first).seq2); }

	if (diag == (LeftInfluence->d).end()) {
		fprintf(stderr, "\nIn HandleOneIntersectionPoint::The diagonal does not exist Point is %lld %lld", (head->first).seq1, (head->first).seq2);
		exit(0);
	}

	diagInter = (LeftInfluence->di).find(diag->first);
	delOwner = diag->second;

	leftOwner = delOwner;
	leftOwner--;
	rightOwner = delOwner;
	rightOwner++;

	if (*leftOwner == *rightOwner) {
		fprintf(stderr, "\nIn HandleOneIter:: The leftOwner is the same as the right owner");
		exit(0);
	}

	if (LI_Winner(LeftInfluence, *leftOwner, *rightOwner)) {
		//the diagonal continues
		if (DEBUG) { fprintf(stderr, "\nIn HandleOneIter:: Diagonal continues"); }
		diag->second = col->second;
		nextCol = col;
		nextCol++;
		nextColInter = (LeftInfluence->ci).find(nextCol->first);
		(LeftInfluence->c).erase(col);
		//FIX #8 MAJOR FIX
		(LeftInfluence->ci).erase(colInter);

		if (nextCol != (LeftInfluence->c).end()) {
			// the column exists
			if (nextColInter->second == inter.end()) {
				// the column is not involved in an intersection
				diagInter->second = inter.end();
				CreateIntersectionPoint(LeftInfluence, nextCol->first, diag->first, nextColInter, diagInter);
			} else {
				//should unset the diagonal
				diagInter->second = inter.end();
			}
		} else {
			diagInter->second = inter.end();
		}
	} else {
		if (DEBUG) { fprintf(stderr, "\nIn HandleOneIter Column continues %f %f %f", (*delOwner)->score, (*leftOwner)->score, (*rightOwner)->score); }

		prevDiag = diag;
		prevDiag--;
		prevDiagInter = (LeftInfluence->di).find(prevDiag->first);

		(LeftInfluence->d).erase(diag);
		(LeftInfluence->di).erase(diagInter);

		if (prevDiag != (LeftInfluence->d).end()) {
			if (prevDiagInter == (LeftInfluence->di).end()) {
				fprintf(stderr, "\nIn HandleOneIter:No diag inter corresponding to  PrevDiag: %lld", prevDiag->first);
				exit(0);
			}

			if (prevDiagInter->second == inter.end()) {
				// the diagonal is not involved in an intersection
				colInter->second = inter.end();
				CreateIntersectionPoint(LeftInfluence, col->first,prevDiag->first, colInter, prevDiagInter);
			} else {
				//should unset the column flag
				colInter->second = inter.end();
			}
		} else {
			colInter->second = inter.end();
		}
	}

	//delete the owner
	(LeftInfluence->o).erase(delOwner);

	inter.erase(inter.begin());
}


long long int printDBound(LI * LeftInfluence) {
	if (DEBUG) { return 0; }
	DBound::iterator i;
	long long int diagCount = 0;
	fprintf(stderr, "\nThe DBound is ::");

	for (i = (LeftInfluence->d).begin(); i != (LeftInfluence->d).end(); i++) {
		fprintf(stderr, "%lld ", i->first);
		diagCount++;
	}

	fprintf(stderr, "Dbound Done/n");
	return diagCount;
}


long long int printCBound(LI * LeftInfluence) {
	if (DEBUG) { return 0; }
	CBound::iterator i;
	long long int colCount = 0;
	fprintf(stderr, "\nThe CBound is ::");

	for (i = (LeftInfluence->c).begin(); i != (LeftInfluence->c).end(); i++) {
		fprintf(stderr, "%lld ", i->first);
		colCount++;
	}

	fprintf(stderr, "Cbound Done/n");
	return colCount;
}


long long int printOwners(LI * LeftInfluence) {
	if (DEBUG) { return 0; }
	Owner::iterator i;
	long long int ownerCount = 0;
	fprintf(stderr, "\nThe Owner is ::");

	for (i = (LeftInfluence->o).begin(); i != (LeftInfluence->o).end(); i++) {
		ownerCount++;
		fprintf(stderr, "%f ", (*i)->score);
	}

	fprintf(stderr, "Owners Done/n");
	return ownerCount;
}


void printState(LI * LeftInfluence) {
	if (DEBUG) { return; }
	long long int colCount, diagCount, ownerCount;

	fprintf(stderr, "\nCurrent State:\n");
	ownerCount = printOwners(LeftInfluence);
	colCount = printCBound(LeftInfluence);
	diagCount = printDBound(LeftInfluence);
	interPointPrint();
}


void interPointPrint() {
	if (DEBUG) { return; }
	InterPoint::iterator i;
	fprintf(stderr, "\nThe Inter is ::");
	for (i = inter.begin(); i != inter.end(); i++) {
		fprintf(stderr, "%lld %lld  ", (i->first).seq1, (i->first).seq2);
	}
	fprintf(stderr, "Inter Done/n");
}
