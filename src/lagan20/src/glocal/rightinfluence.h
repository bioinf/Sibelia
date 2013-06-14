#ifndef RIGHTINFLUENCE
#define RIGHTINFLUENCE

#include<structs.h>
#include<io.h>
#include<score.h>


struct longlongCompare {
	bool operator()(long long int p1,long long int p2) const {
		if (p1 < p2) {
			return 1;
		} else {
			return 0;
		}
	}
};


typedef  map<const long long int , Fragment*,longlongCompare> Active;

typedef struct RI {
  //List of active regions
  Active  act;
  long long int scoreIndex;
  long long int reflectFlag;   
} RI;


void initRI(RI *RightInfluence,long long int scoreIndex);
float lookUpScore(RI * RightInfluence,Fragment *current);
Fragment* lookUpOwnerEnd(RI * RightInfluence,Fragment *current);
Fragment* lookUpOwnerStart(RI * RightInfluence,Fragment *current);
long long int RIWinner(RI *RightInfluence,Fragment *first,Fragment * second);
//long long int processRowofEndPoints(RI *RightInfluence,long long int firstIndex);
long long int diagonal(Fragment * current,RI * RightInfluence);
Fragment * nextOnActive(RI* RightInfluence,Fragment * current);
long long int printActive(RI * RightInfluence);
long long int RICommitEndPoint(RI *RightInfluence,Fragment *current);


#endif
