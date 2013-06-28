#!/usr/bin/python

from argparse import ArgumentParser
import os.path
import os
import re
import sys
import subprocess

FILENAME = "variant_ann.vcf"
script_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib', 'Sibelia', 'snpEff')

parser = ArgumentParser(description='Script for variants annotation')
parser.add_argument("-i", action="store", dest="source", help="source vcf file with variants")
parser.add_argument("-o", action="store", dest="dest", help="destination directory for annotated variants file")
parser.add_argument("--db", action="store", dest="db", help="snpEff database name")
parser.add_argument("-c", action="store", dest="config", help="config file for snpEff")

args = parser.parse_args()

# define default names if needed
if (not args.source):
	args.source = "./variant.vcf"
args.source = os.path.abspath(args.source)

if (not os.path.exists(args.source)):
	print "Please specify source variants file"
	sys.exit(-1)

if (not args.dest):
	args.dest = "annotation"
if (not args.config):
	args.config = script_dir + "/snpEff.config"
args.config = os.path.abspath(args.config)

if (not os.path.exists(args.dest)):
	os.mkdir(args.dest)

if (not args.db):
# get genome and cromosome name from vcf file
	chrom_name = ""
	source_file = open(args.source)
	for line in source_file:
		if (line[:11] == "##reference"):
			assembly_name = line.strip().split("=")[1]
			chrom_match = re.search("\|(\w+)(\.\d+)?\|$", assembly_name)
			chrom_name = assembly_name if (not chrom_match) else chrom_match.group(1)
			break
	source_file.close()

# find snpEff database name
	genomes_file = open(script_dir + "/genomes.txt")
	for line in genomes_file:
		fields = line.strip().split("/")
		if (fields[-1] == chrom_name + ".val"):
			args.db = fields[1]
	if (not args.db):
		print "Couldn't get database name from vcf, please provide it manually"
		sys.exit(-1)

print args

os.chdir(args.dest)
if (not os.path.exists("snpEff_v3_1_" + args.db + ".zip")):
	if (subprocess.call("java -jar " + script_dir + "/snpEff.jar download -c " 
			+ args.config + " " + args.db, shell=True) != 0):
		print "Database was not loaded"
		sys.exit(-1)

dest_file = open(FILENAME, "w")
subprocess.call("java -jar " + script_dir + "/snpEff.jar eff -c " + args.config
	+ " " + args.db + " " + args.source, stdout=dest_file, shell=True)
