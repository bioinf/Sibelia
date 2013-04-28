#!/usr/bin/python

from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import os.path
import re
import sys

ME_THRESHOLD = 2
MAX_INDENT_THRESHOLD = 85

class Block:
	def __init__(self, id, strand, start = 0, end = 0):
		self.id = id
		self.strand = strand
		self.start = start
		self.end = end
	
	def __str__(self):
		return ("+" if self.strand else "-") + str(self.id) + ", start=" + str(self.start) + ", end=" + str(self.end)
	
	def __eq__(self, other):
		return self.id == other.id and self.strand == other.strand

	def reverse(self):
		return Block(self.id, not self.strand, self.start, self.end)

class Permutation:
	def __init__(self, name, id, permutation = []):
		self.name = name
		self.id = id
		self.permutation = permutation
	
	def __str__(self):
		res = "Name: " + self.name + " " + str(self.id) + "\n" + "Permutation:"
		for block in self.permutation:
			res += str(block) + " "
		return res
	
	def reverse(self):
		return [block.reverse() for block in self.permutation[::-1]]
	
	@staticmethod
	def reverse_seq(seq):
		return [block.reverse() for block in seq[::-1]]

	@staticmethod
	def compare_permutations(seq1, seq2):
		if seq1 == seq2:
			return 0
		elif seq1 == Permutation.reverse_seq(seq2):
			return 1
		else:
			return -1

	def find_mobile_elements(self, mobile_elements):
		# very simple version, find single repeating blocks
		repeats = {} # block_id -> block
		for i, block in enumerate(self.permutation):
			if block.id in repeats.keys():
				repeats[block.id].append(block)
			else:
				repeats[block.id] = [block]
		for block_id, entry in repeats.iteritems():
			if len(entry) >= ME_THRESHOLD:
				for block in entry:
					insert_to_mobile_elements([Block(block.id, block.strand)], mobile_elements, self.id, block.start, block.end)

class MobileElement:
	def __init__(self, locations, blocks_seq=[], genome_seq="", descriptions=[]):
		self.blocks_seq = blocks_seq
		self.locations = locations
		self.genome_seq = genome_seq
		self.descriptions = descriptions

	def __str__(self):
		res = "Blocks:\n"
		for block in self.blocks_seq:
			res += str(block) + " "
		res += "Locations:" + str(self.locations) + "\n"
		res += "Description:" + str(self.descriptions) + "\n"
		res += "Sequence:" + self.genome_seq
		return res

def blast_sequence(sequence):
	blast_handle = NCBIWWW.qblast("blastn", "nr", sequence, megablast = True)
	blast_records = NCBIXML.parse(blast_handle)
	insertion_hits = []
	for record in blast_records:
		for aln in record.alignments:
			if aln.title.find("insert") > 0 or aln.title.find("transpos") > 0:
				insertion_hits.append(aln)
	return insertion_hits

def insert_blast_res_to_mobile_elements(blast_results, mobile_elements, permutations, seq_id):
	permutation = next(p for p in permutations if p.name == seq_id)
	for hit in blast_results:
		for hsp in hit.hsps:
			start_block = next((i for i, block in enumerate(permutation.permutation) if block.start >= hsp.query_start), None)
			end_block = next((len(permutation.permutation) - 1 - i for i, block in enumerate(permutation.permutation[::-1]) 
				if block.end <= hsp.query_start + hit.length), None)
			if start_block and end_block and start_block >= end_block:
				insert_to_mobile_elements(permutation.permutation[start_block:end_block + 1], mobile_elements, permutation.id, 
					permutation.permutation[start_block].start, permutation.permutation[end_block].end, hit.title)
			else:
				mobile_elements.append(MobileElement(locations = {permutation.id: [(hsp.query_start, hsp.query_start + hit.length, True)]}, descriptions = [hit.title]))

def insert_to_mobile_elements(blocks_seq, mobile_elements, permutation_id, start = 0, end = 0, description = ""):
	if not start and not end:
		start = blocks_seq[0].start
		end = blocks_seq[-1].end
	for me in mobile_elements:
		compare_res = Permutation.compare_permutations(me.blocks_seq, blocks_seq)
		if compare_res == -1:
			continue
		else:
			if permutation_id not in me.locations.keys():
				me.locations[permutation_id] = []
			if compare_res == 1:
				me.locations[permutation_id].append((start, end, False))
			elif compare_res == 0:
				me.locations[permutation_id].append((start, end, True))
			if description:
				me.descriptions.append(description)
		return
	mobile_elements.append(MobileElement(blocks_seq = blocks_seq, locations = {permutation_id: [(start, end, True)]}, descriptions = [description]))

def get_descriptions(blocks_seq, mobile_elements):
	for me in mobile_elements:
		if any((blocks_seq == me.blocks_seq[i : i + len(blocks_seq)]) for i in xrange(0, len(me.blocks_seq) - len(blocks_seq) + 1)):
			return me.descriptions
	return None

def fill_mobile_element_seq(sequences, mobile_elements):
	for me in mobile_elements:
		if me.descriptions:
			continue
		descriptions = get_descriptions(me.blocks_seq, mobile_elements)
		if (descriptions):
			me.descriptions = descriptions
			continue
		for block in me.blocks_seq:
			me.genome_seq += sequences[block.id - 1] if block.strand else sequences[block.id - 1].reverse_complement()

def add_block_start_end(permutations, block_id, seq_id, start, end, prev = 0):
	permutation = next(p for p in permutations if p.name == seq_id)
	(move, block) = next(((i, b) for i, b in enumerate(permutation.permutation[prev:]) if b.id == block_id), (0, None))
	if not block:
		print "It should never happen"
		sys.exit(-1)
	block.start = start if start <= end else end
	block.end = end if start <= end else start
	prev += move + 1
	return prev

def parse_permutation_string(perm_str):
	res = []
	for block in perm_str.split(' ')[:-1]:
		res.append(Block(int(block[1:]), block[0] == '+'))
	return res

def parse_permutations_file(filename, permutations):
	perm_file = open(filename)
	for l in perm_file:
		if l[0] == ">":
			permutations.append(Permutation(l.strip()[1:], len(permutations)))
		else:
			permutations[-1].permutation = parse_permutation_string(l.strip())

def parse_sequences_file(filename, sequences, permutations):
	prev = 0
	prev_seq_id = ""
	for block_record in SeqIO.parse(filename, "fasta"):
		block_params = {param.split('=')[0]: param.split('=')[1] for param in block_record.id.split(',')}
		block_id, strand, seq_id, start, end = (int(block_params['Block_id']), (block_params['Strand'] == '\'+\''), block_params['Seq'][1:-1], int(block_params['Start']), int(block_params['End']))

		if block_id == len(sequences) + 1:
			sequences.append(block_record.seq if strand else block_record.seq.reverse_complement())
			prev = 0
		if prev_seq_id != seq_id:
			prev_seq_id = seq_id
			prev = 0

		prev = add_block_start_end(permutations, block_id, seq_id, start, end, prev)

def parse_assembly_file(filename, sequences):
	for record in SeqIO.parse(filename, "fasta"):
		sequences.append(record)

if (__name__ == "__main__"):

	script_dir = os.path.dirname(sys.argv[0])

	parser = ArgumentParser(description='Script for mobile elements annotation')
	parser.add_argument("-p", action="store", dest="permutations", help="genome permutations file")
	parser.add_argument("-s", action="store", dest="sequences", help="block sequences file")
	parser.add_argument("-a", action="store", dest="assembly", help="FASTA file with assembly")
	parser.add_argument("-o", action="store", dest="output", help="output file")

	args = parser.parse_args()

	# define default names if needed
	if (not args.permutations):
		args.permutations = "./genomes_permutations.txt"

	if (not os.path.exists(args.permutations)):
		print "Please, specify permutations file"
		sys.exit(-1)

	if (not args.sequences):
		args.sequences = "./blocks_sequences.fasta"

	if (not os.path.exists(args.sequences)):
		print "Please, specify sequences file"
		sys.exit(-1)

	if (not args.output):
		args.output = "./mob_el.txt"
	
	permutations = []
	mobile_elements = []
	blocks_sequences = []
	contigs = []

	parse_permutations_file(args.permutations, permutations)
	parse_sequences_file(args.sequences, blocks_sequences, permutations)
	parse_assembly_file(args.assembly, contigs)

	for i, seq_record in enumerate(contigs):
		if (len(seq_record.seq) < 70000):
			print seq_record.id + " start"
			blast_results = blast_sequence(seq_record.seq)
			insert_blast_res_to_mobile_elements(blast_results, mobile_elements, permutations, seq_record.id)			
			print seq_record.id + " done"
	
	for me in mobile_elements:
		print me

	for permutation in permutations:
		permutation.find_mobile_elements(mobile_elements)

	for me in mobile_elements:
		print me

	fill_mobile_element_seq(blocks_sequences, mobile_elements)

	print "After fill"

	for me in mobile_elements:
		print me
