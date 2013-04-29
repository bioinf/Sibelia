Basic usage
===========
Directory "examples/C-Sibelia/Staphylococcus_aureus" contains a set of two
bacterial genomes. The easiest way to run "C-Sibelia" is to type:

	C-Sibelia -s fine -r <a FASTA file with reference> -a <a FASTA file with a genome>

For example (for two genomes in "examples/C-Sibelia" dir):

	C-Sibelia -s fine -r NCTC8325.fasta -a RN4220.fasta

The command above will run "C-Sibelia" on two bacterial genomes with the "fine"
simplification parameters. After this run you will find a file "variant.vcf" 
with found differences between these two genomes in VCF format. You also can
get additional information about synteny blocks between these two genomes, to
do this run "C-Sibelia" with the "-o" option set:

	C-Sibelia -s fine -r NCTC8325.fasta -a RN4220.fasta -o out

To understand the information in these files, see "Directory for output files"
section. Please note that the file with reference must contain only one
sequence, while another file may contain multiple sequences (contigs).

There is another simplification parameters set, called "loose". To understand
the difference between them, please read manual for "Sibelia" (SIBELIA.md),
section "Basic usage".

Technical parameters
====================

Directory for output files
--------------------------
If you want to see files, describing synteny blocks between input genomes, set
the "outdir" option:

	-o <dir name> or --outdir <dir name>

"C-Sibelia" will write these files there. If this option is set, VCF file with
variants will also be written in output directory.

Output description
==================
By default, "Sibelia" writes only VCF file with found variants ("variant.vcf").
Additionally, it can produce output similar to Sibelia (see section "Directory
for output files"):

1. Blocks coordinates
2. Genomes represented as permutations of the synteny blocks
3. Coverage report
4. Files for generating a "Circos" picture
5. Interactive html-diagram of synteny blocks
6. Sequences file

All these files are described below in details.

Variant report
--------------
Output file "variants.vcf" contains variants found between two genomes. The 
variants are presented in VCF format, version 4.1. If you want "C-Sibelia" to
write variants to some other file, set the "variant" option":

	-v <vile name> or --variant <file name>

If "-o" option is set, then VCF file be written in the output directory.

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
this parameter is 500.

Parameters set
--------------
See the corresponding section in SIBELIA.md.

Maximum number of iterations
----------------------------
See the corresponding section in SIBELIA.md.
