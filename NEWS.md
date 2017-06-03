Sibelia 3.0.7
============
* Fixed a compilation error for newer versions of GCC

Sibelia 3.0.6
=============
* Fixed a bug with crash when 0 simplificaiont stages specified
* When both --graphfile" and --visualize are set, graph is output at each stage
* Added commandline option --noblocks

Sibelia 3.0.5
=============
* Removed commandline options: -a, --alignment
* Added commandline options: --maf
* Added the third parameters set "far"
* Changed permissions of LAGAN scripts to run on some platforms
* Changed structure of the directory "share" (Linux platforms only)
* Fixed a bug with misorder of aligned repeats 
  (affected -a option of C-Sibelia only)
* Fixed some compiler warnings yielded by LAGAN sources
* Fixed a bug in hierarchy circos output that resulted in uncompilable
  configs for some datasets

Sibelia 3.0.4
=============
* Released postprocessing scripts for circos pictures (see CIRCOS_HELPER.md)
* Fixed a few bugs caused by large amount of 'N' characters in the input
* Fixed crash of C-Sibelia with large block size threshold

Sibelia 3.0.3
=============
* Updated snpEff to version 3.3f
* C-Sibelia now produces alignments of repeats within genomes

Sibelia 3.0.2
=============
* Fixed a bug in FASTA parser caused by Windows line endings on Unix-like
systems

Sibelia 3.0.1
=============
* Fixed a bug in C-Sibelia: incorrect reporting of variants within synteny
  blocks having opposite directions
* Added output of synteny blocks in GFF format
* Removed sequences without synteny blocks from Circos diagrams
* Blocks orientations are depicted on Circos diagram by different colors
* Added Sibelia cmd options: --nopostprocess, --allstages, --correctboundaries,
  --gff, --lastk
* Added C-Sibelia cmd options: -o, -a
* C-Sibelia outputs whole genome alignments in XMFA format (-o flag)
* Slightly changed "fine" parameters simplification set

Sibelia 3.0.0
=============
* C-Sibelia is functional now
* Variant annotation with snpEff
* Directory structure on Unix system is more conventional

Sibelia 2.1.2
=============
* Fixed crashes on some input genomes

Sibelia 2.1.1
=============
* Got rid of overlap between synteny blocks, overlapping regions are granted to
  the blocks with higher multiplicity
* Added generation of circos picture that shows blocks at multiple stages
  simultaneously
* Added cmd option: -v

Sibelia 2.1.0
=============
* Improved performance in the case of many genomes. Now it is possible to
  compute synteny blocks for 59 E.coli genomes within 70 minutes and 8 GB of
  RAM on a standard desktop
* Synteny block boundaries now are more precise

Sibelia 2.0.1
=============
* Improved diagram appearance
* Fixed bug with multiple input files
* Fixed circos generation in case of many input sequencs
* Changed format of coordinates of synteny blocks, see USAGE.md, section
  "Blocks coordinates"

Sibelia 2.0.0
=============
* Improved performance
* Completed manual
* Added d3 visualization
* Removed cmd options: -c -p -d
* Changed behaviour of the options: -q -g -r
* Added cmd options: -t -o

Sibelia 1.0.0
=============
* Initial release
