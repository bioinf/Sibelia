#include <stdio.h>
#include <string.h>

void Add_Tick(char *line, int count, int length);
void Print_Lines(char *line1, char *line2, char *ticks1, char *ticks2,
  char *match);
int Usage(void);

char MyName[1024];

int main(int argc, char **argv) {
  FILE *infile = NULL;
  FILE *snp_file = NULL;
  char *slash;
  int fields, start = -1, end = -1, bp, base1, base2;
  int base1_count = 0;
  int base2_count = 0;
  int start2 = 0;
  int end2 = 0;
  int tick1_done = 0;
  int tick2_done = 0;
  int width = 60;
  int length = 0;
  int html_length = 0;
  int snp_pos = -1;
  int param1 = 1;
  char bases[] = {'-', 'A', 'C', 'T', 'G', 'N'};
  char line1[1024];
  char line2[80];
  char match[80];
  char ticks1[80] = "";
  char ticks2[80] = "";
  char snp_fname[1024] = "";
  char font_start[80] = "<b><font color=red ";
  char font_end[] = "</font></b>";
  char status_start[] = "onmouseover=\"window.status='SNP: ";
  char status_end[] = "'\" onmouseout=\"window.status=''\">";
  char dash[] = " - ";
  char snp_bases[2];

// remove the directory name from the program pathname

  if (((slash = strrchr(argv[0], '/')) != NULL) ||
      ((slash = strrchr(argv[0], '\\')) != NULL))
    strcpy(MyName, slash + 1);
  else
    strcpy(MyName, argv[0]);

// parse my command line and open input file(s)

  if (argc < 2) return Usage();
  if (argv[1][0] == '-')
    if (strcasecmp(argv[1], "-pga") == 0)
      ++param1;
    else if (strcmp(argv[1], "-") != 0)
      return Usage();
  if ((argc <= param1) ||
      ((strcmp(argv[param1], "-") != 0) &&
      ((infile = fopen(argv[param1], "r")) == NULL)) ||
      ((argc > (param1 + 1)) &&
      (((fields = sscanf(argv[param1 + 1], "%d", &start)) != 1) ||
      (start <= 0))) ||
      ((argc > (param1 + 2)) &&
      (((fields = sscanf(argv[param1 + 2], "%d", &end)) != 1) ||
      (start > end))))
    return Usage();
  if (infile == NULL)
    infile = stdin;
  else if (param1 > 1) {
    if (((slash = strrchr(argv[param1], '/')) != NULL) ||
        ((slash = strrchr(argv[param1], '\\')) != NULL)) {
      strncpy(snp_fname, argv[param1], slash - argv[param1] + 1);
      snp_fname[slash - argv[param1] + 1] = '\0';
    }
    strcat(snp_fname, "SNP.txt");
    snp_file = fopen(snp_fname, "r");
  }
  while (!feof(infile)) {
    if ((bp = getc(infile)) == EOF) {  // get next char
      if (!ferror(infile)) {
        end2 = base2_count;
        continue;
      }
      perror("Error reading file");  // stop if an error is found
      return 1;
    }
    // decode bp char
    base1 = bp >> 4;
    base2 = bp & 0xf;
    if (base1 != 0) {
      ++base1_count;
      tick1_done = 0;
    }
    if (base2 != 0) {
      ++base2_count;
      tick2_done = 0;
    }
    if (base1_count < start) continue;
    if (snp_file != NULL) {
      while (base1_count > snp_pos) {
        if ((fields = fscanf(snp_file, "%d %2c", &snp_pos, snp_bases)) == 2)
	  continue;
	fclose(snp_file);
	snp_file = NULL;
	break;
      }
    }
    if (start2 == 0) {
      start2 = base2_count;
      if (base2 == 0) ++start2;
    }
    if (base1_count != snp_pos) {
      line1[html_length] = bases[base1];
      line1[html_length + 1] = 0;
      ++html_length;
    } else {
      strcpy(line1 + html_length, font_start);
      strcat(line1, status_start);
      html_length = strlen(line1);
      line1[html_length] = snp_bases[0];
      strcpy(line1 + html_length + 1, dash);
      line1[html_length + strlen(dash) + 1] = snp_bases[1];
      strcpy(line1 + html_length + strlen(dash) + 2, status_end);
      html_length = strlen(line1);
      line1[html_length] = bases[base1];
      strcpy(line1 + html_length + 1, font_end);
      html_length = strlen(line1);
    }
    line2[length] = bases[base2];
    line2[length + 1] = 0;
    match[length] = ((base1 == base2) && (base1 != 5)) ? '|' : ' ';
    match[length + 1] = 0;
    ++length;
    if ((tick1_done == 0) && ((base1_count % 10) == 0) && (base1_count > 0)) {
      Add_Tick(ticks1, base1_count, length);
      tick1_done = 1;
    }
    if ((tick2_done == 0) && ((base2_count % 10) == 0) && (base2_count > 0)) {
      Add_Tick(ticks2, base2_count, length);
      tick2_done = 1;
    }
    if (length == 60) {
      Print_Lines(line1, line2, ticks1, ticks2, match);
      length = 0;
      html_length = 0;
    }
    if (base1_count == end) {
      end2 = base2_count;
      break;
    }
  }
  if (length != 0)
    Print_Lines(line1, line2, ticks1, ticks2, match);
  fclose(infile);
  if (param1 > 1)
    printf("start2=%d\nend2=%d\n", start2, end2);
  return 0;
}

void Add_Tick(char *line, int count, int length) {
  int space;
  char tick[20];
  
  sprintf(tick, "%d", count);
  space = length + 9 - strlen(line) - strlen(tick);
  if (space > 0) {
    while (space > 0) {
      strcat(line, " ");
      --space;
    }
    strcat(line, tick);
  }
}

void Print_Lines(char *line1, char *line2, char *ticks1, char *ticks2,
    char *match) {
  printf("\n%s\nseq1     %s\n         %s\nseq2     %s\n%s\n",
    ticks1, line1, match, line2, ticks2);
  line1[0] = line2[0] = ticks1[0] = ticks2[0] = match[0] = 0;
}

int Usage() {
  fprintf(stderr, " \
Usage: %s [-pga] { - | alignment_file } [start [end]]\n",
    MyName);
  return 1;
}
