#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define MAX_CELLS ((long long int) 100000000)
#define MAX_TIME ((long long int) 100000 * (long long int) 100000)

int failed = 0;

void getFileInfo (char *filename, int *numContigs, int *seqLen, int *numHits){
  FILE *file;
  int dummy, i;
  
  if (!(file = fopen (filename, "r"))){
    fprintf (stderr, "contigorder: Error opening file: %s\n");
    exit (1);
  }
  
  fscanf (file, "numContigs = %d\n", numContigs);
  fscanf (file, "seqLen = %d\n", seqLen);
  
  *numHits = 0;
  while (!feof (file)){
    if (fscanf (file, "(%d %d)", &dummy, &dummy) == 2){
      for (i = 0; i < *numContigs; i++){
	fscanf (file, "%&d", &dummy);
      }
      while (fgetc (file) != '\n');
      (*numHits)++;
    }
  }

  fclose (file);
}

void getScores (char *filename, int numContigs, int seqLen, int numHits, int ***score, int ***ranges){
  FILE *file;
  int i, j;

  *score = (int **) malloc (sizeof (int *) * numHits); 
  assert (*score);
  *ranges = (int **) malloc (sizeof (int *) * numHits);
  assert (*ranges);
  for (i = 0; i < numHits; i++){
    (*score)[i] = (int *) calloc (numContigs, sizeof (int));
    assert ((*score)[i]);
    (*ranges)[i] = (int *) calloc (2, sizeof (int));
    assert ((*ranges)[i]);
  }

  if (!(file = fopen (filename, "r"))){
    fprintf (stderr, "contigorder: Error opening file: %s\n");
    exit (1);
  }
  
  fscanf (file, "numContigs = %*d\n");
  fscanf (file, "seqLen = %*d\n");
  
  i = 0;
  while (!feof (file) && i < numHits){    
    if (fscanf (file, "(%d %d)", &((*ranges)[i][0]), &((*ranges)[i][1])) == 2){
      for (j = 0; j < numContigs; j++){
	fscanf (file, "%d", &((*score)[i][j]));
      }
      while (fgetc (file) != '\n');
      i++;
    }
  }

  fclose (file);
}

void floodfill (int *labels, int *first, int *last, int numContigs, int here, int groupNum){
  int i;

  labels[here] = groupNum;
  for (i = 0; i < numContigs; i++){
    if (i != here && labels[i] == -1 && first[i] != -1){
      if (!(first[here] > last[i] || last[here] < first[i])){
	floodfill (labels, first, last, numContigs, i, groupNum);
      }
    }
  }
}

int *getLabels (int **score, int numContigs, int numHits){
  int *labels, *first, *last, i, j;
  
  labels = (int *) calloc (numContigs, sizeof (int)); assert (labels);
  first = (int *) calloc (numContigs, sizeof (int)); assert (first);
  last = (int *) calloc (numContigs, sizeof (int)); assert (last);

  for (j = 0; j < numContigs; j++){
    first[j] = -1;
    for (i = 0; i < numHits; i++){
      if (score[i][j] > 0){
	if (first[j] == -1) first[j] = i;
	last[j] = i;
      }
    }
  }

  j = 0;
  for (i = 0; i < numContigs; i++) labels[i] = -1;
  for (i = 0; i < numContigs; i++){
    if (labels[i] == -1 && first[i] != -1){
      floodfill (labels, first, last, numContigs, i, j++);
    }
  }
  
  free (first);
  free (last);
  return labels;
}

int makeRanges (int **score, int numHits, int *cols, int numCols, int **first, int **last){
  int i, j, k, found, numRanges = 1;
  
  for (i = 0; i < numHits; i++){
    for (j = 0; j <= i; j++){
      for (k = found = 0; !found && k < numCols; k++){
	found = (score[i][cols[k]] > 0) && (score[j][cols[k]] > 0);
      }
      if (found) numRanges++;
    }
  }

  *first = (int *) calloc (numRanges, sizeof (int)); assert (*first);
  *last = (int *) calloc (numRanges, sizeof (int)); assert (*last);
  
  (*first)[0] = -1; // initial range
  (*last)[0] = -1; // initial range
  numRanges = 1;

  for (i = 0; i < numHits; i++){
    for (j = 0; j <= i; j++){
      for (k = found = 0; !found && k < numCols; k++){
	found = (score[i][cols[k]] > 0) && (score[j][cols[k]] > 0);
      }
      if (found){
	(*first)[numRanges] = j;
	(*last)[numRanges] = i;
	numRanges++;
      }
    }
  }

  return numRanges;
}

int **calcRangeScores (int **score, int *cols, int numCols, int *first, int *last, int numRanges){
  int i, j, k, **scoreOf;
  
  scoreOf = (int **) malloc (sizeof (int *) * numCols); assert (scoreOf);
  for (i = 0; i < numCols; i++){
    scoreOf[i] = (int *) malloc (sizeof (int) * numRanges); assert (scoreOf[i]);
    for (j = 0; j < numRanges; j++){
      scoreOf[i][j] = 0;
      
      if (j > 0){
	for (k = first[j]; k <= last[j]; k++){
	  scoreOf[i][j] += score[k][cols[i]];
	}
      }
    }
  }

  
  return scoreOf;
}

void solveOrder (int **score, int numContigs, int numHits, int *cols, int numCols, int **ranges,
		 int **results, int *resultCtr){
  int i, j, k, l, m;
  int numStates = (1 << numCols), numRanges;
  int **best, *first, *last, ptr, newScore, **scoreOf;
  int bestScore = 0, bestState, bestRange, newBest, addedScore;
  int *stateList, *rangeList, *scoreList;
  int work, totwork;

  numRanges = makeRanges (score, numHits, cols, numCols, &first, &last);

  if ((long long int) numRanges * (long long int) numStates > MAX_CELLS ||
      (long long int) numRanges * (long long int) numStates * (long long int) numCols * (long long int) numRanges > MAX_TIME){
    fprintf (stderr, "ordering failed, retrying... (numRanges = %d, numStates = %d)\n", numRanges, numStates);
    printf ("ordering failed\n");
    failed = 1;
    return;
  }

  best = (int **) malloc (sizeof (int *) * numStates); assert (best);
  for (i = 0; i < numStates; i++){    
    best[i] = (int *) calloc (numRanges, sizeof (int)); assert (best[i]);
  }
  for (i = 0; i < numStates; i++) best[i][0] = 0;
  for (j = 1; j < numRanges; j++) best[0][j] = 0;

  scoreOf = calcRangeScores (score, cols, numCols, first, last, numRanges);

  // -- DP solution ---------------

  work = 0;
  totwork = (numRanges - 1) * (numStates - 1);

  // search over all state transitions
  for (i = 1; i < numRanges; i++){
    for (j = 1; j < numStates; j++){
      newBest = -1;

      // compute best previous state
      for (k = 0; k < numCols; k++) if (j & (1 << k)){
	m = j - (1 << k);
	addedScore = scoreOf[k][i];
	for (l = 0; l < numRanges; l++) if (last[l] < first[i]){
	  newScore = best[m][l] + addedScore;
	  if (newScore > newBest){
	    newBest = newScore;
	  }
	}	
      }

      best[j][i] = newBest;
    
      if (best[j][i] > bestScore){
	bestScore = best[j][i];

	bestState = j;
	bestRange = i;
      }
      work++;
      if ((work % 100000) == 0){
	fprintf (stderr, "WORKING %d/%d\n", work, totwork);
      }
    }
  }

  // -- Compute traceback ---------

  l = 0;
  stateList = (int *) calloc (numCols, sizeof (int)); assert (stateList);
  rangeList = (int *) calloc (numCols, sizeof (int)); assert (rangeList);
  scoreList = (int *) calloc (numCols, sizeof (int)); assert (scoreList);

  while (bestState != 0){

    k = 1;
    for (i = 0; k && i < numCols; i++) if (bestState & (1 << i)){
      m = bestState - (1 << i);
      for (j = 0; k && j < numRanges; j++) if (last[j] < first[bestRange]){
	newScore = best[m][j] + scoreOf[i][bestRange];
	if (newScore == best[bestState][bestRange]){
	  stateList[l] = cols[i];
	  rangeList[l] = bestRange;
	  scoreList[l] = scoreOf[i][bestRange];
	  l++;
	  bestState = m;
	  bestRange = j;
	  k = 0;
	}
      }
    }
  }

  // -- Report traceback ----------

  for (i = l - 1; i >= 0; i--){
    results[*resultCtr][0] = stateList[i];
    results[*resultCtr][1] = ranges[first[rangeList[i]]][0];
    results[*resultCtr][2] = ranges[last[rangeList[i]]][1];
    results[*resultCtr][3] = scoreList[i];
    (*resultCtr)++;
  }

  for (i = 0; i < numCols; i++) free (scoreOf[i]);
  free (scoreOf);
  for (i = 0; i < numStates; i++) free (best[i]);
  free (best);
  free (first);
  free (last);
  free (stateList);
  free (rangeList);
  free (scoreList);
}

int compFn (const void *a, const void *b){
  return (*(int **) a)[1] - (*(int **) b)[1];
}

void findGroups (int numContigs, int seqLen, int numHits, int **score, int **ranges){
  int *labels, group, pos, i;
  int *columns, **results, resultCtr = 0;

  labels = getLabels (score, numContigs, numHits);
  columns = (int *) malloc (sizeof (int) * numContigs); assert (columns);
  results = (int **) malloc (sizeof (int *) * numContigs); assert (results);
  for (i = 0; i < numContigs; i++){
    results[i] = (int *) calloc (4, sizeof (int)); assert (results[i]);
  }
  
  group = pos = 0;
  while (!failed){
    for (i = 0; i < numContigs; i++){
      if (labels[i] == group)
	columns[pos++] = i;
    }
    if (pos == 0) break;
    solveOrder (score, numContigs, numHits, columns, pos, ranges, results, &resultCtr);
    pos = 0;
    group++;
  }

  if (!failed){
    qsort (results, resultCtr, sizeof (int *), compFn);
    for (i = 0; i < resultCtr; i++){
      printf ("%d --> (%d %d) %d\n", results[i][0], results[i][1], results[i][2], results[i][3]);
    }
  }

  for (i = 0; i < numContigs; i++) free (results[i]);
  free (results);
  free (labels);
  free (columns);
}

int main (int argc, char **argv){
  int numContigs, seqLen, numHits, i;
  int **score, **ranges;

  if (argc != 2){
    fprintf (stderr, "Usage:\ncontigorder rangefile\n");
    exit (1);
  }
  
  getFileInfo (argv[1], &numContigs, &seqLen, &numHits);
  
  //fprintf (stderr, "numContigs = %d, seqLen = %d, numHits = %d\n", numContigs, seqLen, numHits);
  
  getScores (argv[1], numContigs, seqLen, numHits, &score, &ranges);
  findGroups (numContigs, seqLen, numHits, score, ranges);
  
  for (i = 0; i < numHits; i++){
    free (score[i]);
    free (ranges[i]);
  }
  free (score);
  free (ranges);  
  
  return 0;
}
 
