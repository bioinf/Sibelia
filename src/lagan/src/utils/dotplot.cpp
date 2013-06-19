#include <stdlib.h>
#include <stdio.h>

int main (int argc, char **argv){
  FILE *file;
  int s1b, s1e, s2b, s2e, pa, pb, maxa = 0, maxb = 0;
  float score;
  char buffer[105];
  char* name1 = NULL;
  char* name2 = NULL;
  char dummy[] = "unknown";
  int PAD, PAD2;

  if (argc < 2){
    fprintf (stderr, "Usage: dotplot anchfile [name1 [name2]] \n");
    exit(1);
  }
  
  if (argc > 2) name1 = argv[2];
  if (argc > 3) name2 = argv[3];
  if (name1 == NULL) name1 = dummy;
  if (name2 == NULL) name2 = dummy;

  pa = -1;
  pb = -1;
  
  file = fopen (argv[1], "r");
  while (!feof (file)){
    if (fscanf (file,
        "(%d %d)=(%d %d) %f", &s1b, &s1e, &s2b, &s2e, &score) == 5 &&
        s2b > 0){
      if (s1b > maxa) maxa = s1b;
      if (s1e > maxa) maxa = s1e;
      if (s2b > maxb) maxb = s2b;
      if (s2e > maxb) maxb = s2e;
    }
    fgets (buffer, 105, file);
  }
  fclose (file);
//  PAD = maxa / 1000;
//  PAD2 = maxb / 1000;

  file = fopen (argv[1], "r");
  printf ("set nokey\n");
  printf ("set xlabel \"%s\"\n", name1);
  printf ("set ylabel \"%s\"\n", name2);
  printf ("set title \"Dotplot: %s vs. %s\"\n", name1, name2);
  printf ("set style line 1 linetype 3 linewidth 3\n");
  printf ("set style line 2 linetype 1 linewidth 4\n");


  while (!feof (file)){
    if (fscanf (file,
        "(%d %d)=(%d %d) %f", &s1b, &s1e, &s2b, &s2e, &score) == 5 && s2b > 0){
      if (s1b > maxa) maxa = s1b;
      if (s1e > maxa) maxa = s1e;
      if (s2b > maxb) maxb = s2b;
      if (s2e > maxb) maxb = s2e;

      if (s2b < s2e){
	// draw forward aligns
	PAD = (s1e-s1b)* 2/10;
	PAD2 = (s2e-s2b)* 2/10;
	printf ("set arrow from %d,%d to %d,%d nohead ls 1\n",
	  s1b-PAD, s2b-PAD2, s1e+PAD, s2e+PAD2);

	// draw connections
	// if (pa != -1 && pb != -1)
	//      printf ("set arrow from %d,%d to %d,%d nohead lt -1 lw 0.01\n", pa, pb, s1b, s2b);
	pa = s1e;
	pb = s2e;
      }
    }
    fgets (buffer, 105, file);
  }
  fclose (file);

  file = fopen (argv[1], "r");
  while (!feof (file)){
    if (fscanf (file,
        "(%d %d)=(%d %d) %f", &s1b, &s1e, &s2b, &s2e, &score) == 5 && s2b > 0){
      if (s2b > s2e){
	// draw rev aligns
	PAD = (s1e-s1b)* 2/10;
	PAD2 = (s2b-s2e)* 2/10;
	printf ("set arrow from %d,%d to %d,%d nohead ls 2\n",
	  s1b-PAD2, s2b+PAD2, 
		s1e+PAD2, s2e-PAD2);

	// draw connections
	// if (pa != -1 && pb != -1)
	//      printf ("set arrow from %d,%d to %d,%d nohead lt -1 lw 0.01\n", pa, pb, s1b, s2b);
	pa = s1e;
	pb = s2b;
      }
    }
    fgets (buffer, 105, file);
  }

  printf ("plot [1:%d][1:%d] -1\n", maxa * 11/10, maxb*11/10);
  printf ("set terminal postscript enhanced color\n");
  printf ("set output \"sin.ps\"\n");
  printf ("replot\n");


  fclose (file);
}
