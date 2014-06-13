Basic usage
===========
In this manual it is assumed that "C-Sibelia" is properly installed and the
directory with "C-Sibelia.py" is in your "PATH" environment variable or input
files are in the same folder with "C-Sibelia.py".

Please note that "C-Sibelia" is available for "Unix"-like OS only, use 
virtual machine or our web-server to run it:

	http://etool.me/software/sibelia

"C-Sibelia" finds a set of differences between a reference and an assembly. All
of them are reported relative to the given reference. The easiest way to run
"C-Sibelia" is to type:

	C-Sibelia.py <FASTA file with reference> <FASTA file with an assembly>

For example, directory "examples/C-Sibelia/Staphylococcus_aureus" contains a
set of two bacterial genomes. To run "C-Sibelia" on this dataset, type:

	C-Sibelia.py NCTC8325.fasta RN4220.fasta

The command above will run "C-Sibelia" on two bacterial genomes with the "fine"
simplification parameters. After this run you can find a file "variant.vcf" 
with detected differences between these two genomes in VCF format. You can also
obtain alignments between these two genomes (see section "Alignment output").

There is another simplification parameters set, called "loose". To understand
the difference between them, please read manual for "Sibelia" (SIBELIA.md),
section "Basic usage".

"C-Sibelia" requires some free disk space to write it's temporary files. By
default, they're stored within the current working directory. If you want to
store temporary files in some other location, use the "-t" option.

You can also obtain synteny blocks which were used by "C-Sibelia" in variant
calling. To do this, run "C-Sibelia" with "-o" option set.

Output description
==================

Variant output
--------------
Output file "variants.vcf" contains variants found between two genomes. The 
variants are presented in VCF format, version 4.1. If you want "C-Sibelia" to
write variants to some other file, set the "variant" option":

	-v <file name> or --variant <file name>

If you specify the "-o" option, then variant file will be put in that directory.

Unmapped insertions
-------------------
Sometimes it is impossible to determine exact position of variation relative to
a given reference. Particularly, suppose that you compare an assembly in form
of contigs against a fully assembled reference. And of the contigs doesn't 
contain synteny blocks at all. It seem to be a novel insertion, but how can
it's position relative to the reference be determined? Such insertions are
shown in VCF file on the first lines using breakends (see specification of the
VCF format for reference) which are located arbitrarily within the first 
sequence from the reference genome. However, this way of reporting such
variants may be unconvenient, so if you want to see unmapped insertions in the
FASTA format, use the key:

	-u <file name> or --unmapped <file name>

If you specify the "-o" option, then insertion file will be put in that
directory.

Alignment output
----------------
You can also obtain alignment of obtained synteny blocks:
	
	--maf <file_name>

If you specify the "-o" option, then alignment file will be put in that
directory. Please note that these alignments may also indicate intergenomic
repeats. "C-Sibelia" output was verified by the validator from "mafTools"
(https://github.com/dentearl/mafTools/).

Directory for synteny blocks output
-----------------------------------
By default, "C-Sibelia" prodcues only a single VCF file. However, you can also
store the synteny block information. To do this, use the "-o" cmd parameter:

	-o <dir name> or --outdir <dir name>

The directory <dir name> will contain output produced by "Sibelia" that was
used for variant calling.

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
