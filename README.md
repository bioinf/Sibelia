Sibelia 3.0.7

Release date: 2th June 2017

Authors
=======

* Ilya Minkin (St. Petersburg Academic University)
* Nikolay Vyahhi (St. Petersburg Academic University)
* Mikhail Kolmogorov (St. Petersburg Academic University)
* Ekaterina Starostina (St. Petersburg Academic University)
* Son Pham (University of California, San Diego)

Introduction
============
This package contains two programs:

* Sibelia -- "Sibelia" is a tool for finding synteny blocks in closely related
genomes, like different strains of the same bacterial species. It takes a set
of FASTA files with genomes and locates coordinates of the synteny blocks in
these sequences. It also represents genomes as permutations of the blocks.

* C-Sibelia -- This tool is designed for comparison between two genomes
represented either in finished form or as sets of contigs. It is able to detect
SNPs/SNVs and indels of different scales. "C-Sibelia" works by locating synteny
blocks between the input genomes and aligning different copies of a block.
It considers only unique blocks, i.e. blocks having one copy in reference and
one copy in another genome. It is also possible to get the alignments itself,
in this case you will get the alignments of repeats within genomes as well.

It also contains a script for annotation of variants found by "C-Sibelia" using
the "snpEff" tool.

Installation
============
See INSTALL.md file.

Usage
=====
See SIBELIA.md for "Sibelia", C-SIBELIA.md for "C-Sibelia" and "ANNOTATION.md"
for the annotation script.

System requirements
===================
This version of "C-Sibelia" supports only "Unix"-like operating systems, but
"Sibelia" runs fine on "Windows". To use "C-Sibelia", "Windows" users may use a
virtual machine or try our web server:

	http://etool.me/software/csibelia

"Windows" support will be retained soon in the future releases. Binary releases
for "Windows" contain only "Sibelia" binaries.

"C-Sibelia" requires "Python" and "Perl" to be istalled in your system. Please
note that "C-Sibelia" runs only under "Python2" of version at least 2.7. The
annotation script requires "Java" runtime. 

Citation
========
If you use "Sibelia" in your research, please cite:

Ilia Minkin, Anand Patel, Mikhail Kolmogorov, Nikolay Vyahhi, and Son Pham. "Sibelia: a scalable and comprehensive synteny block generation tool for closely related microbial genomes." In Algorithms in Bioinformatics, pp. 215-229. Springer Berlin Heidelberg, 2013.

License
=======
"Sibelia" is distributed under GNU GPL v2 license, see LICENSE.

It also uses third-party librarires and programs:
* libdivsufsort (MIT License), author Yuta Mori
https://code.google.com/p/libdivsufsort
* TCLAP (MIT License), authors Michael E. Smoot and Daniel Aarno 
http://tclap.sourceforge.net
* Boost (Boost Software License)
http://www.boost.org
* D3.js (BSD License)
http://d3js.org
* Seqan (BSD/3-clause)
http://www.seqan.de

Contacts
========
E-mail your feedback at ivminkin@gmail.com.

You also can report bugs or suggest features using issue tracker at GitHub
https://github.com/bioinf/Sibelia/issues
