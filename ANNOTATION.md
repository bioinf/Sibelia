Variants annotation
===================
The script "snpEffAnnotate.py" is designed for variant annotation and effect
prediction for variants found by C-Sibelia. It can be launched without 
arguments from directory with C-Sibelia results:
	
	snpEffAnnotate.py

If you want to specify your own vcf file as input for this tool, type:

	snpEffAnnotate.py -i variants.vcf

By default, all outup files are located in './annotation' directory, name of
annotated vcf is variant_ann.vcf. You can specify other directory with the "-o"
option, if specified directory is not exists it will be created. For example:
	
	snpEffAnnotate.py -i variants.vcf -o annotation

Description of output file format can be found at

	http://snpeff.sourceforge.net/SnpEff_manual.html#output
	
Take into account, that by default, script uses ##assembly field from vcf file
to get NCBI database name. For example, a possible value of this field can be
"##assembly=gi|57116681|ref|NC_000962.2|". 

One more important note about vcf is that #CHROM column must contain cromosome
name exactly the same as in snpEff database. For assembly name above it will be
"NC_000962".

If you wish, you can specify snpEff database name with the "--db" argument.
