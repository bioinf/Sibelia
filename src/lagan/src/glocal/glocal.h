#ifndef GLOCAL
#define GLOCAL

#define DEBUG 1

#ifndef LLONG_MAX
// limits.h entries from ISO C99
#define LLONG_MAX 9223372036854775807LL
#define LLONG_MIN (-LLONG_MAX - 1LL)
#endif

#include<structs.h>
#include<io.h>
#include<rightinfluence.h>
#include<leftinfluence.h>
#include<score.h>

long long int startPointHandler();
long long int endPointHandler();
float fragmentSetScore(Fragment * current,Fragment *owner);
void intersectionPointHandler();

#endif
