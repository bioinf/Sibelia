Basic usage
===========
Directory "examples/C-Sibelia" contains a set of two bacterial genomes. The easiest
way to run "C-Sibelia" is to type:

	C-Sibelia.exe -s fine -r <a FASTA file with reference> -a <a FASTA file with a genome>

For example (for two genomes in "examples/C-Sibelia" dir):

	C-Sibelia.exe -s fine -r nctc8325.fasta -a RN4220S.fasta

The command above will run "C-Sibelia" on two bacterial genomes with the "fine"
simplification parameters. After this run you will find a file "variant.vcf" 
with found differences between these two genomes in VCF format. You'll also
find some other files that contain information about synteny blocks between
these two genomes.

Please note that the file with reference must contain only one sequence, while
another file may contain multiple sequences (contigs).

There is another simplification parameters set, called "loose". To understand
the difference between them, please read manual for "Sibelia" (SIBELIA.md),
section "Basic usage".

Technical parameters
====================

Directory for output files
--------------------------
See the corresponding section in SIBELIA.md.

Output description
==================
By default, "Sibelia" produces 5 files: 

1. Blocks coordinates
2. Genomes represented as permutations of the synteny blocks
3. Coverage report
4. Files for generating a "Circos" picture
5. Interactive html-diagram of synteny blocks
6. Sequences file
7. Variant report

All these files are described below in details.

Variant report
--------------
File "variants.vcf" contains variants found between two genomes. The variants
are presented in VCF format.

Blocks coordinates
------------------
See the corresponding section in SIBELIA.md.

Genomes permutations
--------------------
See the corresponding section in SIBELIA.md.

Coverage report
---------------
See the corresponding section in SIBELIA.md.

Sequences file
--------------
See the corresponding section in SIBELIA.md.

"d3" visualization
------------------
See the corresponding section in SIBELIA.md.

"Circos" visualization
----------------------
See the corresponding section in SIBELIA.md.

Fine tuning
===========
Here we will describe parameters that can affect computation results.

Minimum block size
------------------
See the corresponding section in SIBELIA.md. For C-Sibelia default value of
this parameter is 1000.

Parameters set
--------------
See the corresponding section in SIBELIA.md.

Maximum number of iterations
----------------------------
See the corresponding section in SIBELIA.md.
