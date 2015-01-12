Basic usage
===========
In this manual it is assumed that "Sibelia" is properly installed and the
directory with "Sibelia.py" is in your "PATH" environment variable or input
files are in the same folder with "Sibelia" executable.

Directory "examples/Sibelia" contains two sets of bacterial genomes. The easiest
way to run "Sibelia" is to type:

	Sibelia -s loose <input FASTA file(s)>

For example:

	Sibelia -s loose Helicobacter_pylori.fasta

Important note -- "Sibelia" requires some free space on HDD to run. If you
experience any problems and see error messages that mention temporary files,
try to change directory used for temporary files (see section "Directory for
temporary files").

Above commands will run "Sibelia" on the file "Helicobacter_pylori.fasta" with
the "loose" simplification parameters. There are two another simplification
parameters sets, called "fine" and "far". To run "Sibelia" on "Helicobacter_pylori.fasta"
with "fine" parameters set, type:

	Sibelia -s fine Helicobacter_pylori.fasta

The difference between "loose" and "fine" set is that "loose" usually produces
fewer blocks, but longer. And it may lose some small synteny blocks 
(shorter than 15 000 BP), while "fine" option produces more blocks (but shorter)
and their coverage is worse. Usually "loose" is the best choice, but if you do
not want to lose information about small-scale rearrangements,  use "fine". See
also section "Output description" for detailed depiction of the output format.
The "far" parameters set can be used for analysis of distantly-related genomes.
However, usage of this set requires more computational time and space and may
be not appropriate in case of many genomes.

If you are not satisfied by the results (poor coverage, for example), try to set
simplification parameters manually (see section "Fine tuning"). 

By default, Sibelia filters out synteny blocks shorter than 5 000 BP. You can
change this behaviour, see section "Minimum block size".

You also may be interested in blocks that occur exactly once in each input
sequence, such blocks are used as input for MGRA algorithm, for example. To get
such output use option "-a":

	Sibelia -s loose -a Helicobacter_pylori.fasta

Synteny blocks are visualized with an interactive diagram (see section "d3"
visualization"). The blocks are also can be visualized with "Circos" (see 
section "Circos" visualization"). While "Circos" is better for publications,
"d3" diagram is more suitable for analysis. Please note that sequences that do
not contain any synteny blocks instances are not shown on these diagrams.

Genomes from the "examples/Sibelia" dir were taken from [5, 6]. Note that you
can specify multiple FASTA files, just separate them with spaces.

Technical parameters
====================

Directory for output files
--------------------------
Default directory = "." You can change this by setting cmd parameter:

	-o <dir name> or --outdir <dir name>

By default, "Sibelia" places output files in current working directory. Setting
this parameter will change output directory.

Directory for temporary files
-----------------------------
Default directory is output directory. You can change this by setting cmd parameter:

	-t <dir name> or --tempdir <dir name>

"Sibelia" creates some temporary files while running. By default these files
are placed in the output directory. If you want to place temporary files in
another folder due to some reasons, use this parameter. Although the files
exist only for a very short period of time, they can be quite big -- ~20*N
bytes, where N is the total size of all input genomes. You can also use switch

	-r or --inram
	
that will force "Sibelia" to not create any temporary files and store all it's
data in RAM.

Output description
==================
By default, "Sibelia" produces following files: 

1. Blocks coordinates
2. Genomes represented as permutations of the synteny blocks
3. Coverage report
4. Files for generating a "Circos" picture
5. Interactive html-diagram of synteny blocks

There are also optional output files:

1. Sequences file
2. Dot file with resulting de Bruijn graph

All these files are described below in details.

Blocks coordinates
------------------
File name = "blocks_coords.txt". First part of this file lists input sequences,
their IDs, sizes and descriptions. IDs are just index numbers of sequences (in
the same order as they apper in the input files).

Second part of the file describes synteny blocks in sections separated by 
dashes. Each section starts with the synteny block ID. Then all instances
of this block are listed in tabular format. Each row of a table depicts
one instance of this synteny block. Columns of the table designate following:

1. Seq_id -- ID of the sequence, that the instance belongs to
2. Strand -- strand of the synteny block instance, either '+' or '-'. Input
sequences are treated as positive strands
3. Start -- one-based index of the first base pair of the instance
4. End -- one-based index of the last base pair of the instance
5. Length -- length of the instance of the synteny block

Note that all coordinates are given relatively to POSITIVE strand of the
sequence. If an instance of a synteny block is located on the positive strand,
then start < end, otherwise start > end. 

If you wish, you can obtain coordinates in GFF format. If you use the flag:

	-gff

Then coordinates of synteny blocks are listed in file "blocks_coords.gff". For
description of this format, see:

	https://cgwb.nci.nih.gov/FAQ/FAQformat.html#format3

Each record represents different copies of a synteny block. Copies having the
same number in the "tag" field (last column) are instances of the same synteny
block.

Genomes permutations
--------------------
File name = "genomes_permutations.txt".

This file contains input sequences represented as permutations of the synteny
block. It has two lines for each input sequence:

1. Header line -- FASTA header of the sequence (starting with '>')
2. Genome line -- sequence of synteny blocks instances as they appear on the 
positive strand of the sequence. Each instance is represented as a signed
integer, where '+' sign depicts direct version of the block, and '-' depicts
reversed block

Coverage report
---------------
File name = "coverage_report.txt".

The file describes portion of the genomes, that found synteny block cover.
First part of the file describes input sequencess (see "Blocks coordinates"
section). Second part of the file is a table with the following columns:

1. Degree -- multiplicity of a synteny block. For example, if a synteny block
has degree = 3, then the are three instances of this block in the input
sequence
2. Count -- number of synteny blocks with a given degree
3. Total -- portion of all the input sequences that are covered by the blocks
with a given degree
4. Seq <n> -- portion of the sequence with id <n> that is covered by blocks
with a given degree

Here is an example of such table from a report file.

| Degree | Count | Total   | Seq 1   | Seq 2   | Seq 3   | Seq 4   |
| :----- | :---: | :-----: | :-----: | :-----: | :-----: | :-----: |
| 2	 | 11	 | 3.82%   | 6.59%   | 2.41%   | 2.96%   | 3.30%   |
| 3	 | 4	 | 1.68%   | 2.24%   | 2.19%   | 1.44%   | 0.85%   |
| 4	 | 21	 | 91.93%  | 91.34%  | 94.71%  | 87.54%  | 94.53%  |
| All	 | 36	 | 95.66%  | 97.44%  | 97.98%  | 90.67%  | 96.89%  |

This table contains one row for each degree (2, 3, 4) and one ("All") row for
the overall coverage. It means that there are 11 blocks with degree = 2, i.e.
11 * 2 instances, and they cover 3.82% of all four genomes, 6.59% of Seq 1 and
2.41% of Seq 2. And also there are 21 synteny blocks with degree = 4, i.e.
4 * 21 instances and they cover 91.93% of all genomes. All the blocks cover
95.66% of all the input sequences.

Sequences file
--------------
File name = "blocks_sequences.fasta". By default this file is not written. To
output this file, set cmd parameter:

	-q or --sequencesfile

This FASTA file contains sequences of instances of the synteny block. Each
sequence has header in following format:

	>Seq="<header>",Strand=<sign>,Block_id=<block id>,Start=<x>,End=<y>

Where "<header>" is a header of the FASTA sequence where the block instance
is located. Other fields are described in section "Coordinates file".
Sequences of synteny blocks are also written in SAM format in the file
"blocks_sequences.sam".

"d3" visualization
------------------
File name = "d3_blocks_diagram.html".

"Sibelia" generates an interactive html diagram that shows found synteny blocks.
Coordinates follow the same convention as described in section "Coordinates 
file".

"Circos" visualization
----------------------
You can visualize synteny blocks with a colorful circular diagram by using
the "Circos" software [3]. Files for generating such diagram are written in
directory "circos" inside the output directory. To generate "Circos" diagram
do following:

1. Download and install "Circos" software
2. Run "Sibelia"
3. Run "Circos" in "circos" directory

For example of such diagrams (generated from "Helicobacter_pylori.fasta),
see "examples/Sibelia/Helicobacter_pylori/circos/circos.png". Also note that
such diagrams can become very piled with larger number of genomes. To overcome
this, plot only big blocks, see section "Minimum block size". Blocks located
on the positive strand are colored green, while blocks from negative strand
are red.

By default, "Sibelia" plots only blocks obtained after the last stage.
You can also view blocks at the intermediate stages by using switch:

	-v or --visualize

On the resulting diagram the outermost circle shows blocks obtained at the
first stage, then the second stage and so on. Please note that this option
slows down the computation.

Resulting de Bruijn graph
-------------------------
File name = "de_bruijn_graph.dot". By default this file is not written. To
output this file, set cmd parameter:

	-g or --graphfile
	
If you are a curious person, you can also view condensed de Bruijn graph that
is used for generating synteny blocks. To understand the graph, see [1].
Condensed means that only bifurcations in the graph are plotted and 
non-branching paths are collapsed into a single edge. Blue edges are generated
from positive strand and red edges are from negative strand respectively. Note
that this graph is generated for K = min(Kn, MinimumBlockSize) or for value of
"--lastk" cmd parameter if it is set, where Kn is the value of K used for the
last stage (see section "Fine tuning").

If one is interested in graph output only, he or she can use the "--noblocks"
option:

	--noblocks

In this case Sibelia doesn't compute the synteny blocks and doesn't ouput them,
but can output the graph. For example, to get only non-modified compressed de
Bruijn grahp for k=25, one can use the following command line:

	Sibelia -k run.stage --noblocks -g -m 25 <input_file>

Where "run.stage" contains single number "0".

Fine tuning
===========
Here we will describe parameters that can affect computation results.

Minimum block size
------------------
Default value = 5000. To change this value, set cmd parameter:

	-m <integer> or --minblocksize <integer>

If you are interested only in big synteny blocks, like > 100 000 BP, set
this parameter to an appropriate value.

Outputting only shared blocks
-----------------------------
Default = not set. Add flag to cmd parameters to set:

	-a or --sharedonly

Output only blocks that occur exactly once in each input sequence. This option
assumes that all input genomes contain single chromosome.

Postprocessing
--------------
By default, synteny blocks are postprocessed after computation by gluing
"stripes" consisting of the same synteny blocks. For example, if each 
occurence of synteny block 1 is followed by synteny block 2 and vice a versa,
their directions are consistent, then they are "glued" together to form a
single synteny block. Postprocessing could be turned off by specifying flag:

	--nopostprocess

Output blocks from all stages
-----------------------------
"Sibelia" performs computations in multiple stages. Every stage produces it's
own synteny blocks. You can get blocks from all stages by specifying flag:

	--allstages

Files "blocks_coordsN.(txt|gff)" will contain their coordinates, where N is the
number of stage. Zero corresponds to blocks obtained without any simplification.

Boundaries correction
---------------------
Algorithm of "Sibelia" depends on presence of solid k-mers within syntenic
regions. If such regions contain variations close to their borders, they will
be truncated. In case of two genomes, some synteny blocks (of multiplicity two)
could be corrected using local alinment algorithms. Use flag:

	--correctboundaries

This flag is used by "C-Sibelia".

Parameters set
--------------
Default value = not set. To select the parameters set, use cmd parameter:

	-s <loose|fine|far> or --parameters <loose|fine>

This option is incompatible with "-k", you must specify one of these, not both.
Approach used in "Sibelia" is parameter dependent. To understand the details,
please see the next section and [1]. The "loose" option produces longer blocks
and better coverage, while "fine" can capture small-scale rearrangements, for
example, inversions of size < 15 000 BP. The "far" set is for distant genomes.

Custom parameters set
---------------------
Default value = not set. To specify the file that contains custom parameters
set, use cmd parameter:

	-k <file name> or --stagefile <file name>

The algorithm consists of several stages of computations. Each stage has two 
parameters, K and D. Let's call K-mer a substring of length K. At each stage
"Sibelia" constructs so called de Bruijn graph, graph of K-mers that occur
in the genome, and simplifies it by removing special type of undirected cycles
called "bulges", see [1] for more details.

Graph is a good model for describing the algorithm, but to understand "physical
meaning" of the parameters it is useful to consider operations that are 
actually performed with the genome behind the graph model. Suppose that 
somewhere in the genome exist two pairs of K-mers K1 and K2:

1st pair: ... K1 ABCD K2 ...  
2nd pair: ... K1 FGHE K2 ...  

If the distance between K1 and K2 within each pair is less than D, then "Sibelia"
replaces FGHE with ABCD to obtain longer "synteny block":

1st pair: ... K1 ABCD K2 ...  
2nd pair: ... K1 ABCD K2 ...  

More concrete example. Suppose that K = 3, D = 5 and somewhere in the genome we
find:

1st pair: ... act gaga ggc ...  
2nd pair: ... act gatg ggc ...  

As we see, distance between "act" and "ggc" is less than 5 nucleotides so we
replace "gatg" by "gaga":

1st pair: ... act gaga ggc ...  
2nd pair: ... act gaga ggc ...  

"Sibelia" keeps track of all changes so it is able to locate original locations
of the synteny blocks obtained by the simplification. This process continues 
step by step, we start with small values of K to obtain longer K-mers shared
between synteny regions and then increase K and D. The "loose" parameters set
has 4 stages:

| K        | D         |
| :------- | --------: |
| 30       | 150       |
| 100      | 1000      |
| 1000     | 5000      |
| 5000     | 15000     |

The "fine" set consists of 3 stages and it's final values are less:

| K        | D         |
| :------- | --------: |
| 30       | 150       |
| 100      | 500       |
| 500      | 1500      |

The "far" set is for distant genomes:

| K        | D         |
| :------- | --------: |
| 15       | 120       |
| 100      | 500       |
| 500      | 1500      |

As you can see, "loose" set is more aggressive -- at it's final stage it glues
together 5000-mers that are separated from each other by at most 15000 symbols.
Although this description is very simplified and lacks many important technical
details, it is enough to infer your own parameter set. Stage file that you may
use to specify your own parameters has following simple format:

M  
K1 D1  
K2 D2  
...  
KM KM  

Where M is the number of stages. So, running with the stage file:

4  
30 150  
100 1000  
1000 5000  
5000 15000  

Is equivalent to running with the -s "loose" cmd option. As you may notice, the
algorithm relies on exact K-mers shared between the genomes. If input genomes
doesn't have such shared substrings, then "Sibelia" won't be able to locate the
synteny blocks.

If you cannot find synteny blocks with the default parameters, try to start
with smaller values of K (~20), increase D values or vary number of stages.

Last value of K
---------------
In "Sibelia" de Bruijn graph is constructed (N + 1) times, where N is the
number of simplification stages. De Bruijn graph which is constructed last is
used for inferring synteny blocks. Value of K for this grahp is determined by
min(Kn, MinimumBlockSize)  where Kn is the value of K used for the last stage.
You can override this value by setting the "--lastk" parameter:

   --lastk <integer > 1>


Maximum number of iterations
----------------------------
Default value = 4. Tho change this value, set cmd parameter:

	-i <integer> or --maxiterations <integer>

Maximum number of iterations during a stage of simplification. Increasing
this parameter may slightly improve coverage.

References
==========
1. Ilya Minkin, Nikolay Vyahhi, Son Pham. "SyntenyFinder: A Synteny Blocks 
Generation and Genome Comparison Tool" (poster), WABI 2012
http://bioinf.spbau.ru/sites/default/files/SyntenyFinder.pdf
2. Max A. Alekseyev and Pavel A. Pevzner. "Breakpoint graphs and ancestral
genome reconstructions", Genome Res. 2009. 19: 943-957.
3. Circos. http://circos.ca
4. D3. http://d3js.org/
5. Helicobacter pylori. http://www.ncbi.nlm.nih.gov/genome/169
6. Staphylococcus aureus. http://www.ncbi.nlm.nih.gov/genome/154
