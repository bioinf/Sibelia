#ifndef LEFTINFLUENCE
#define LEFTINFLUENCE

#include<structs.h>
#include<score.h>

struct LI;


struct longlongCompare2
{
 
  bool operator()(long long int p1,long long int p2) const
  {
    if(p1< p2)
      return 1;
    else 
      return 0;
          
  }
};


struct paircomp
{
 
  bool operator()(const Point p1,const Point p2) const
  {
    if(p1.seq1< p2.seq1)
      return 1;
    else if((p1.seq1 == p2.seq1) && (p1.seq2 < p2.seq2))
      return 1;
    else 
      return 0;
          
  }
};




typedef list<Fragment*> Owner;
typedef map <long long int ,Owner::iterator,longlongCompare2> CBound;

typedef multimap <Point ,struct LI *,paircomp> InterPoint;

typedef map <long long int ,InterPoint::iterator,longlongCompare2> CInter;
typedef map <long long int,Owner::iterator,longlongCompare2> DBound;

typedef map <long long int,InterPoint::iterator,longlongCompare2> DInter;



typedef struct LI
{
  Owner o;
  CBound c;
  DBound d;
  CInter ci;
  DInter di;
  long long int scoreIndex;
  long long int reflectFlag;

  
}LI;


extern InterPoint inter;
 




Owner::iterator LILookUpOwnerIterator(LI* LeftInfluence,long long int seq1,long long int seq2) ;
Fragment * LILookUpOwnerStart(LI* LeftInfluence,Fragment *current);
Fragment * LILookUpOwnerEnd(LI* LeftInfluence,Fragment *current);
CBound::iterator LICColumn(LI* LeftInfluence,long long int seq1, long long int seq2);
Fragment *LICOwner(LI* LeftInfluence,long long int seq1, long long int seq2);
Fragment *LIDOwner(LI* LeftInfluence,long long int seq1, long long int seq2);
DBound::iterator LIDDiagonal(LI* LeftInfluence,long long int seq1, long long int seq2);
float LILookUpScore(LI *LeftInfluence,Fragment *current);
void InitLI(LI* LeftInfluence, long long int scoreIndex);
long long int LI_Winner(LI* LeftInfluence,Fragment * first,Fragment * second);
long long int LICommitPoint(LI *LeftInfluence,Fragment *current);
Owner::iterator LI_OwnerInsertAfter(LI* LeftInfluence,Owner::iterator current,Fragment * curfrag);
long long int  LI_CommitDiagonalOwner(LI* LeftInfluence,Fragment *current,Fragment *owner);
long long int  LI_CommitColumnOwner(LI* LeftInfluence,Fragment *current,Fragment *owner);
void CreateIntersectionPoint(LI* LeftInfluence,long long int col,long long int diag,CInter::iterator colInter,DInter::iterator diagInter);
void DeleteIntersectionPoint(InterPoint::iterator tobeerased,CInter::iterator colInter,DInter::iterator diagInter);
void HandleOneIntersectionPoint();

long long int printDBound(LI * LeftInfluence);
long long int printOwners(LI * LeftInfluence);
long long int printCBound(LI * LeftInfluence);
void printState(LI* LeftInfluence);
void interPointPrint();



#endif
