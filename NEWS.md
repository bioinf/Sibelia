Sibelia 2.20
============
* Got rid of overlap between synteny blocks, overlapping regions are granted to
the blocks with higher multiplicity
* Added generation of circos picture that shows blocks at multiple stages
simultaneously
* Added cmd option: -v

Sibelia 2.1.0
=============
* Improved performance in the case of many genomes. Now it is possible to
compute synteny blocks for 59 E.coli genomes within 70 minutes and 8 GB of RAM
on a standard desktop
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