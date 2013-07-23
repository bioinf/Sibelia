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