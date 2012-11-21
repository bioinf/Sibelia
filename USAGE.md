Basic usage
===========
Directory "examples" contains two sets of bacterial genomes. The easiest way 
to run "Sibelia" is to type:

	Sibelia -s loose <input FASTA file(s)>

For example:

	Sibelia -s loose Helicobacter_pylori.fasta

Important note -- "Sibelia" requires some free space on HDD to run. If you
experience any problem and see error messages that mention temporary files,
try to change directory used for temporary files (see section "Directory for
temporary files").

Above commands will run "Sibelia" on the file "Helicobacter_pylori.fasta" with
the "loose" simplification parameters. There is another simplification
parameters set, called "fine". To run "Sibelia" on "Helicobacter_pylori.fasta"
with "fine" parameters set, type:

	Sibelia -s fine Helicobacter_pylori.fasta

The differenece between "loose" and "fine" set is that "loose" usually produces
fewer blocks, but longer. And it may lose some small synteny blocks 
(shorter than 15000 BP), while "fine" option produces more blocks (but shorter)
and their coverage is worse. Usually "loose" is the best choice, but if you do
not want to lose information about small-scale rearrangaments, use "fine". 

If you are not satisfied by the results (poor coverage for example), set
simplification parameters manually (see section "Fine tuning"). 

By default, Sibelia filters out synteny blocks shorter than 5000 BP. You can
change this behaviour, see section "Minimum block size".

You also may be interested in blocks that occur exactly once in each input
sequence, such blocks are used as input for MGRA algorithm, for example. To get
such output use option "-a":

	Sibelia -s loose -a Helicobacter_pylori.fasta

Synteny blocks can be visualized with interactive diagrams. To generate such
diagram use option "--d3" (see section "d3" visualization"):

	Sibelia -s loose --d3 Helicobacter_pylori.fasta

Synteny blocks are also can be visualized with "Circos" (see section "Circos" 
visualization"). 

Genomes from the "examples" dir were taken from [5, 6].

Technical parametes
===================

Directory for temporary files
-----------------------------
Default directory = "." You can change this by setting cmd paramter:

	-t <dir name> or --tempdir <dir name>

"Sibelia" creates some temporary files while running. By default these files
are placed in the current working directory. Altghough the files exist only for
a very short period of time, they can be quite big -- ~ 20 * N bytes, where N
is the total size of all input genomes.

Output description
==================
By default, "Sibelia" produces 3 files: 

1. Blocks coordinates
2. Genomes represented as permutations of the synteny blocks
3. Coverage report

There are also optional output files:

1. Sequences file
2. Files for generating a "Circos" picture
3. Dot file for resulting de Bruijn graph
4. Interactive html-diagram of synteny blocks

All these files are described below in details.

Blocks coordinates
------------------
Default file name = "block_coords.txt". You can change it by setting cmd
parameter:

	-c <file name> or --coordsfile <file name>

First part of this file lists input chromosomes, their IDs, sizes and
descriptions. IDs are just sequence numbers of chromosomes (in the same
order as they apper in input files).

Second part of the file describes synteny blocks in sections separated by 
dashes. Each section starts with the synteny block ID. Then all instances
of this block are listed in tabular format. Each row of a table depicts
one instance of this synteny block. Columns of the table designate
following:

1. Chr_id -- ID of the chromosome, that the instance belongs to
2. Strand -- strand of the synteny block instance, either '+' or '-'. Input
sequences are treated as positive strands of the chromosomes
3. Start -- zero based index of the starting base pair of the instance. All
indices are given relative to POSITIVE strand
4. End -- zero based index of the base pair following last base pair of the
instance
5. Length -- length of the instance of the synteny block

Genomes permutations
--------------------
Default file name = "genomes_permutations.txt". You can change it by setting
cmd parameter:

	-p <file name> or --permfile <file name>

This file contains input chromosomes represented as permutations of the synteny
block. It has two lines for each input chromosome:

1. Header line -- FASTA header of the sequence (starting with '>')
2. Genome line -- sequence of synteny blocks instances as they appear on the 
positive strand of the chromosome. Each instance is represented as a signed
integer, where '+' sign depicts direct version of the block, and '-' depicts
reversed block

Coverage report
---------------
Default file name = "coverage_report.txt". You can change it by setting
cmd parameter:
	
	-r <file name> or  --reportfile <file name>

The file describes portion of the genomes, that found synteny block cover.
First part of the file describes input chromosomes (see "Blocks coordinates"
section). Second part of the file is a table with following collumns:

1. Degree -- multiplicity of the synteny block. For example, if synteny block
has degree = 3, then the are three instances of this block in the input
chromosomes
2. Count -- number of synteny blocks with a given degree
3. Total -- portion of all the input chromosomes that cover blocks with a given
degree
4. Chr <n> -- portion of the chromosome with id <n> that cover blocks with a 
given example

This table contains one row for each degree and one ("All") row for overall
coverage. For example (output from "Helicobacter_pylori.fasta"):

	Degree	Count	Total	Chr 1	Chr 2
	2	23	94.68%	96.40%	93.09%	
	All	23	94.68%	96.40%	93.09%	

It means that there are 23 blocks with degree = 2, i.e. 23 * 2 instances, and
they cover 94.68% of both genomes, 96.40% of Chr 1 and 93.09% of Chr2. Note 
that synteny blocks can overlap (by at most 5000 BP for loose parameter set),
so sum in each column may not equal to value at last row.

Sequences file
--------------
Default file name = not set. To output this file, set cmd parameter 

	-q <file name> or --sequncesfile <file name>

This FASTA file contains sequences of instances of the synteny block. Each
sequence has header in following format:

	>Seq="<description>",Strand=<sign>,Block_id=<block id>,Start=<x>,End=<y>

Where <description> is a header of the FASTA sequence where the block instance
is located. Other fields are described in section "Coordinates file".

"d3" visualization
------------------
Default file name = not set. To output thes file, set cmd parameter

	--d3 <file name>

With this option set "Sibelia" will generate an interactive html diagram that
show found synteny blocks. 

"Circos" visualization
----------------------
Default directory name = not set. To output these files, set cmd parameter

	-d <dir name> or --circosdir <dir name>

You can visualize synteny blocks with a colorful circular diagram by using
the "Circos" software [3]. Do following:

1. Download and install "Circos" software
2. Run "Sibelia" with "-d" options set
3. Run "Circos" on files generated by "Sibelia"

For example, to generate Circos diagram for example "Helicobacter_pylori.fasta"
perform following:

1. Run "Sibelia" with following parameters:

	Sibelia -s loose -d ./circos Helicobacter_pylori.fasta

2. Run circos in the "circos" directory

For example of such diagrams (generated from "Helicobacter_pylori.fasta),
see "examples/Helicobacter_pylori/circos/circos.png". Also note that such
diagrams can become very piled with larger genomes. To overcome this, plot only
big blocks, see section "Minimum block size".

Resulting de Bruijn graph
-------------------------
Default file name = not set. To output this file, set cmd parameter

	-g <file name> or  --graphfile <file name>

If you are a curious person, you can also view condensed de Bruijn graph that
is used for generating synteny blocks. To understand the graph, see [1].
Condensed means that only bifurcations in the graph are plotted and 
non-branching paths are collapsed into a single edge. Blue edges are generated
from positive strand and red edges are from negative strand respectively.

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

Output only blocks that occur exactly once in each input sequence.

Parameters set
--------------
Default value = not set. To select the parameters set, use cmd parameter:

	-s <loose|fine> or --parameters <loose|fine>

This option is incompatible with "-k", you must specify one of these, not both.
Approach used in "Sibelia" is parameter dependent. To understand the details,
please see the next section and [1]. The "loose" option produces longer blocks
and better coverage, while "fine" can capture small-scale rearrangements, for
example, inversions of size < 15000 BP. 

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
step by step, we start with small values of K to obatin longer K-mers shared
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
| 100      | 1000      |
| 1000     | 2500      |

As you can see, "loose" set is more agressive -- at it's final stage it glues
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
