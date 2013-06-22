Basic usage
===========
The easiest way to run "C-Sibelia" is to type:

	C-Sibelia.py <FASTA file with reference> <FASTA file with a genome>

For example, directory "examples/C-Sibelia/Staphylococcus_aureus" contains a
set of two bacterial genomes. To run "C-Sibelia" on this dataset, type:

	C-Sibelia.py NCTC8325.fasta RN4220.fasta

The command above will run "C-Sibelia" on two bacterial genomes with the "fine"
simplification parameters. After this run you can find a file "variant.vcf" 
with detected differences between these two genomes in VCF format.

There is another simplification parameters set, called "loose". To understand
the difference between them, please read manual for "Sibelia" (SIBELIA.md),
section "Basic usage".

"C-Sibelia" requires some free disk space to write it's temporary files. By
default, they're stored within the current working directory. If you want to
store temporary files in some other location, use the '-t' option.

Output description
==================
Output file "variants.vcf" contains variants found between two genomes. The 
variants are presented in VCF format, version 4.1. If you want "C-Sibelia" to
write variants to some other file, set the "variant" option":

	-v <file name> or --variant <file name>

Technical parameters
====================
"C-Sibelia" works by aligning different copies of a synteny block. It relies on
"Sibelia" in finding the blocks. Thus, it accepts some arguments of "Sibelia"
that control synteny block generation.

Directory for temporary files
-----------------------------
Default directory is the current working directory. You can change this by
setting a cmd parameter:

	-t <dir name> or --tempdir <dir name>

"Sibelia" creates some temporary files while running. By default these files
are placed in the output directory. If you want to place temporary files in
another folder due to some reasons, use this parameter. 

Number of processes for alignment
---------------------------------
"C-Sibelia" can perform alignment of synteny blocks using multiple processes
to speed-up. By default, it uses only one process. To change the number of 
concurrent alignment processes, use the "-p" option:

	-p <processes number> or --processcount <processes number>

Minimum block size
------------------
See the corresponding section in SIBELIA.md. For C-Sibelia default value of
this parameter is 500.

Parameters set
--------------
See the corresponding section in SIBELIA.md.

Maximum number of iterations
----------------------------
See the corresponding section in SIBELIA.md.
