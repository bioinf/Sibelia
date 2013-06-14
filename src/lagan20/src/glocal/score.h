#ifndef SCORE
#define SCORE

#include <structs.h>
#include <glocal.h>

#define MAXCASES 20
#define MAXOBJECTS 10

struct LI;
struct RI;

class ScoreInterface {
	protected:
	float openConstant,minConstant,maxConstant,diagConstant;

	ScoreInterface (float iopenConstant , float iminConstant ,float imaxConstant,float idiagConstant);
	float getScore(Fragment *up, Fragment * down){return -1;};
};


class Score :public ScoreInterface {
	public:
	Score(float iopenConstant , float iminConstant ,float imaxConstant,float idiagConstant);

	float getScore(Fragment *up, Fragment * down);
};


void initScoreFunctionPointers(char *scoreFileName);
void  createScoreFunctionObjects(char * line);
long long int charToCase(char in);
float scoreAll(Fragment *up,Fragment *down, long long int ret_case);
long long int Myabs(long long int a);
long long int Mymin(long long int a,long long int b);
long long int Mymax(long long int a,long long int b);
float fragmentSetScore(Fragment * current,Fragment *owner,LI *LeftInfluence, RI * RightInfluence,long long int rightInfluenceFlag);

#endif
