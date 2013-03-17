Supplementary materials for paper "Sibelia: A fast synteny block generation tool for closely related microbial genomes".

Here you can find datasets and command-lines/scripts, which were used to obtain results,
described in paper.


0 Tools preparation
===================

0.0 Sibelia
-----------
Download Sibelia (http://github.com/bioinf/Sibelia) and compile it (as described in INSTALL.md).
To run some of scripts, you need to copy Sibelia executable in script directory, or change 
script variables (see scripts comments).

0.1 Mugsy
---------
Download mugsy (http://mugsy.sourceforge.net/). See installtion guide for futher instructions.
Note, that mugsy is need to be run in default working directory (so you should copy scripts to it)

0.2 TBA
-------
Download/compile multiz-tba (http://www.bx.psu.edu/miller_lab/) and  muscle 
(http://www.drive5.com/muscle/downloads.htm). Our run scripts placed in "mugsy-tba" folder.
Edit installation paths in tbaenv.sh first. To run mugsy in tba mode use --tba option.

0.3 Mauve
---------
Download and compile mauve (http://gel.ahabs.wisc.edu/mauve/).


1 Small synthetic test
======================

Scripts to to run small syntenic test are located in "synthetic_test" folder. See comments inside.
Note, that if you want to compare results or make a circos diagram, you should use some scripts, 
described below.


2 Comparison on 3 e.coli genomes
================================

Scripts located in "3_genomes_test" folder. Usage: "./toolRunner.sh parameters_file genomes_file". 
Files with parameters for corresponding tools and file with genomes are located in same folder.
It`s better to copy these scripts in installation directories of corresponding tools. As a result,
folders with tools output will be produced (named as "tool_name_params"). To compare obtained results,
use scripts described below.

4 Performance comparison
========================

Scripts are located in "performance_test" folder. First, put any number of genomes in "testgenomes" folder
(for example, you can get 59 e.coli genomes from http://www.ncbi.nlm.nih.gov/genome/?term=txid562[orgn]).
Then, run script for corresponding tool. See instructions inside scripts. It`s better to run scripts
from tools installation directories. All tools are run with their default parameters.

3 Comparison scripts
====================

Scripts are located in "scripts" folder. Some of them needs AlignIO branch of Biopython library
(http://biopython.org/wiki/Multiple_Alignment_Format). Folder with cloned git repo is expected to be
named "alignio-maf". Most of them accept files with synteny blocks description, which can be obtained
by conversion scripts (see below). All scripts outputs in stdout.

* sibConvet.py sib_file
Converts sibelia output (blocks_coords.txt) to plain format, recognised by other scripts.

* mafConvert.py maf_file min_block_size strains_file
Converts maf file (from mugsy ot mugsy-tba) to plain format. Needs minimum block size of run
(which can bo found in "blocksize.txt" file in output folder) and file with ordered input genomes names
("strains.txt" file in "3_genomes_test" folder).

* mauveConvert.py xmfa_file min_block_size strains_file
Converts mauve output to plain format. Parameters are described ahead.

* highlight.py plain_file
Converts blocks in plain format to circos highlight file (for making circos diagrams).

* coverage.py file
Calculates synteny blocks coverage rate. Script requires fasta files with genomes. Edit path to them,
if you need.

* f1stat.py file1 file2
Calculates F1-statistic for two pain blocks files.

* boundaryDiscord.py file1 file2
Calculates boundary discordane for two plain blocks files.