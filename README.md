Sibelia 3.0.0 Beta version

Release date: TBD

Authors
=======

* Ilya Minkin (St. Petersburg Academic University)
* Nikolay Vyahhi (St. Petersburg Academic University)
* Mikhail Kolmogorov (St. Petersburg Academic University)
* Son Pham (University of California, San Diego)

Introduction
============
This package contains two programs:

* Sibelia -- "Sibelia" is a tool for finding synteny blocks in closely related
genomes, like different strains of the same bacterial species. It takes a set
of FASTA files with genomes and locates coordinates of the synteny blocks in
these sequences. It also represents genomes as permutations of the blocks.

* C-Sibelia -- This tool is designed for comparison between a reference and
either a finished genome, or a genome represented as a set of contigs. It is
able to detect SNPs/SNVs and indels of different scales. It also finds synteny
blocks between these two genomes.

Installation
============
See INSTALL.md file.

Usage
=====
See SIBELIA.md for "Sibelia" and C-SIBELIA.md for "C-Sibelia".

Citation
========
If you use "Sibelia" in your research, please cite:

Ilya Minkin, Nikolay Vyahhi, Son Pham. "SyntenyFinder: A Synteny Blocks 
Generation and Genome Comparison Tool" (poster), WABI 2012
http://bioinf.spbau.ru/sites/default/files/SyntenyFinder.pdf

License
=======
"Sibelia" is distributed under GNU GPL v2 license, see LICENSE.

It also uses third-party librarires:
* kseq.h (MIT License), author Heng Li
http://lh3lh3.users.sourceforge.net/kseq.shtml
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
E-mail your feedback at ilya.minkinen@gmail.com.

You also can report bugs or suggest features using issue tracker at GitHub
https://github.com/bioinf/Sibelia/issues
