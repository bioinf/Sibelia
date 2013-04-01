#!/usr/bin/python

from argparse import ArgumentParser
import os.path
import re
import sys
import subprocess

script_dir = os.path.dirname(sys.argv[0])

parser = ArgumentParser(description='Script for variants annotation')
parser.add_argument("-i", action="store", dest="source", help="source vcf file with variants")
parser.add_argument("-o", action="store", dest="dest", help="destination path for annotated variants file")
parser.add_argument("--db", action="store", dest="db", help="snpEff database name")
parser.add_argument("-c", action="store", dest="config", help="config file for snpEff")

args = parser.parse_args()

# define default names if needed
if (not args.source):
	args.source = "./variant.vcf"

if (not os.path.exists(args.source)):
	print "Please, specify source variants file"
	sys.exit(-1)

if (not args.dest):
	args.dest = "./variant_ann.vcf"
if (not args.config):
	args.config = script_dir + "/snpEff.config"

if (not args.db):
# get genome and cromosome name from vcf file
	chrom_name = ""
	source_file = open(args.source)
	for line in source_file:
		if (line[:10] == "##assembly"):
			assembly_name = line.strip().split("=")[1]
			chrom_match = re.search("\|(\w+)(\.\d+)?\|$", assembly_name)
			chrom_name = assembly_name if (not chrom_match) else chrom_match.group(1)
	source_file.close()

# find snpEff database name
	genomes_file = open(script_dir + "/genomes.txt")
	for line in genomes_file:
		fields = line.strip().split("/")
		if (fields[-1] == chrom_name + ".val"):
			args.db = fields[1]

print args

if (not os.path.exists("snpEff_v3_1_" + args.db + ".zip")):
	if (subprocess.call("java -jar " + script_dir + "/snpEff.jar download -c " 
			+ args.config + " " + args.db, shell=True) != 0):
		print "Database was not loaded"
		sys.exit(-1)

dest_file = open(args.dest, "w")
subprocess.call("java -jar " + script_dir + "/snpEff.jar eff -c " + args.config
	+ " " + args.db + " " + args.source, stdout=dest_file, shell=True)
