Sibelia 2.1.1

Release date: 21th March 2013

Authors
=======

* Ilya Minkin (St. Petersburg Academic University)
* Nikolay Vyahhi (St. Petersburg Academic University)
* Mikhail Kolmogorov (St. Petersburg Academic University)
* Son Pham (University of California, San Diego)

Introduction
============
"Sibelia" is a tool for finding synteny blocks (regions of conserved DNA)
in closely related genomes, like different strains of the same bacterial
specie. It takes a set of FASTA files with genomes and locates coordinates of
the synteny blocks in these sequences.

System Requirements
===================
This version is designed for small-scale genomes. For example, it can find
synteny blocks in 8 bacterial genomes (51 MB of genomic data) in ~10 minutes
and using 900 MB of RAM and 800 MB HDD space on an Intel i5 laptop. 

Installation
============
See INSTALL.md file.

Usage
=====
See USAGE.md file.

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
