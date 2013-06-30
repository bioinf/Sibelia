#!/usr/bin/python

import re
import os
import sys
import time
import shutil
import tempfile
import argparse
import itertools
import functools
import subprocess
import collections
import multiprocessing

COVER = 1
UNCOVER = 0
MINIMUM_CONTEXT_SIZE = 30
BLOCKS_FILE = 'blocks_sequences.fasta'
INSTALL_DIR = os.path.dirname(os.path.abspath(__file__))
LAGAN_DIR = os.path.join(INSTALL_DIR, '..', 'lib', 'Sibelia', 'lagan')
os.environ['LAGAN_DIR'] = LAGAN_DIR

FastaRecord = collections.namedtuple('FastaRecord', ['seq', 'description', 'id'])

def strip_chr_id(chr_id):
	part = chr_id.split('|')
	if len(part) == 5:
		return part[-2].split('.')[0]
	return chr_id

def parse_fasta_file(file_name):	
	handle = open(file_name)
	line = [line.strip() for line in handle if line.strip() != '']
	record = []
	i = 0
	while i < len(line):
		if line[i][0] == '>':
			j = i + 1
			while j < len(line) and line[j][0] != '>':
				j += 1
			seq = ''.join(line[i + 1:j])
			description = line[i][1:].strip()
			seq_id = description.split()[0]
			record.append(FastaRecord(seq=seq, description=description, id=seq_id))
			i = j
		else:
			i += 1		
	handle.close()
	return record

def write_fasta_records(fasta_record, file_name):
	LINE_LENGTH = 60
	handle = open(file_name, 'w')
	for record in fasta_record:
		print >> handle, '>' + record.description		
		pos = 0
		while pos < len(record.seq):
			end = min(pos + LINE_LENGTH, len(record.seq))
			print >> handle, record.seq[pos:end]
			pos = end
	handle.close()

class Variant(object):
	def __init__(self, reference_chr_id, reference_pos, contig_id, assembly_pos, reference_allele,
				 assembly_allele, reference_context, assembly_context, synteny_block_id):
		self._reference_chr_id = reference_chr_id
		self._reference_pos = reference_pos		
		self._contig_id = str(contig_id)
		self._assembly_pos = assembly_pos
		self._reference_allele = '' if reference_allele is None else reference_allele.upper()  
		self._assembly_allele = '' if assembly_allele is None else assembly_allele.upper()
		self._reference_context = str('' if reference_context is None else reference_context.upper())
		self._assembly_context = str('' if assembly_context is None else assembly_context.upper())
		self._synteny_block_id = synteny_block_id

	def __str__(self):
		return "\t".join([str(self._reference_pos), self._reference_allele, self._assembly_allele,
						str(self._synteny_block_id), self._contig_id,
						self._reference_context, self._assembly_context])
				
	def get_synteny_block_id(self):
		return self._synteny_block_id
	
	def get_assembly_pos(self):
		return self._assembly_pos
	
	def get_reference_context(self):
		return self._reference_context
	
	def get_assembly_context(self):
		return self._assembly_context
			
	def get_reference_pos(self):
		return self._reference_pos
	
	def get_reference_chr_id(self):
		return self._reference_chr_id
	
	def get_contig_id(self):
		return self._contig_id
	
	def get_reference_allele(self):
		return self._reference_allele
	
	def get_assembly_allele(self):
		return self._assembly_allele
	
	def get_vcf_record(self):
		data = [strip_chr_id(self.get_reference_chr_id()), str(self.get_reference_pos() + 1),
			'.', self.get_reference_allele(), self.get_assembly_allele(), '.', '.', '.']
		return '\t'.join(data)

def no_gaps(sequence):
	return ''.join([ch for ch in sequence if ch != '-'])

def parse_blocks_coords(blocks_file):
	line = [l.strip() for l in open(blocks_file) if l.strip()]
	pos = 1
	num_seq_id = dict()
	while pos < len(line) and line[pos][0] != '-':
		line[pos] = line[pos].split()
		num_seq_id[line[pos][0]] = line[pos][2]
		pos += 1
	pos += 1
	block_list = []	
	while pos < len(line):				
		pos += 2
		block = []
		while line[pos][0] != '-':
			instance = line[pos].split()
			chr_id = instance[0]
			strand = instance[1]
			start = int(instance[2])
			end = int(instance[3])
			if strand == '+':
				start -= 1
			else:
				temp = start
				start = end - 1
				end = temp			
			block.append((num_seq_id[chr_id], start, end))			
			pos += 1
		pos += 1
		block_list.append(block)
	return block_list

def get_context(alignment, alignment_segment, segment_index):
	context = []	
	if segment_index > 0:
		segment = alignment_segment[segment_index - 1]
		start = segment[1] - min(segment[1] - segment[0], MINIMUM_CONTEXT_SIZE)
		context.append(str(alignment[0][start:segment[1]]))
	else:
		context.append('')

	if segment_index + 1 < len(alignment_segment):
		segment = alignment_segment[segment_index + 1]
		end = segment[0] + min(segment[1] - segment[0], MINIMUM_CONTEXT_SIZE)
		context.append(str(alignment[0][segment[0]:end]))
	else:
		context.append('')
		
	segment = alignment_segment[segment_index]
	reference_context = context[0] + no_gaps(alignment[0][segment[0]:segment[1]]) + context[1]
	assembly_context = context[0] + no_gaps(alignment[1][segment[0]:segment[1]]) + context[1]
	return reference_context, assembly_context

def parse_alignment(alingment_file_name, reference_chr_id, synteny_block_id, contig_id, reference_start):	
	last_match = None
	start_position = None
	alignment_segment = []
	alignment = [record.seq for record in parse_fasta_file(alingment_file_name)]		
	for now_position, symbol in enumerate(zip(alignment[0], alignment[1])):
		now_match = symbol[0] == symbol[1]
		if last_match is None:
			last_match = now_match
			start_position = 0
		elif last_match != now_match:
			if last_match == False or now_position - start_position >= MINIMUM_CONTEXT_SIZE or start_position == 0:
				alignment_segment.append([start_position, now_position, last_match])
				start_position = now_position
			elif alignment_segment:
				start_position = alignment_segment[-1][0]
				del alignment_segment[-1]
			last_match = now_match
			
	alignment_segment.append([start_position, len(alignment[0]), last_match])
	position = reference_start
	reference_position_map = []
	for symbol in alignment[0]:
		reference_position_map.append(position)
		position += 1 if symbol != '-' else 0
	
	variant = []	
	for segment_index, segment in enumerate(alignment_segment):
		start, end, match = segment
		if match == False:
			shift = 1
			variant_reference_start = reference_position_map[start]
			reference_context, assembly_context = get_context(alignment, alignment_segment, segment_index)
			SNP = end - start == 1 and alignment[0][start] != '-' and alignment[1][start] != '-'
			if start == 0 or SNP:
				shift = 0
			reference_allele = no_gaps(alignment[0][start - shift:end])
			assembly_allele = no_gaps(alignment[1][start - shift:end])
			variant.append(Variant(reference_chr_id, variant_reference_start - shift, contig_id, None,
								reference_allele, assembly_allele, reference_context, assembly_context,
								synteny_block_id))				
	return variant

def get_seq(file_name):
	all_seq = [record for record in parse_fasta_file(file_name)]
	seq = [(record.id, record.seq) for record in all_seq ]	
	return (all_seq[0], dict(seq))

def parse_header(header):
	ret = dict()
	header = header.split(',')
	for item in header:
		item = item.split('=')
		key = item[0]
		value = ''.join([ch for ch in item[1] if ch != "'" and ch != '"'])
		ret[key] = value
	return ret

def find_instance(instance_list, reference_seq_id, in_reference):
	for instance in instance_list:
		header = parse_header(instance.description)
		if (header['Seq'] in reference_seq_id) == in_reference:
			return instance
	return None

def process_unique_block(unique_block, block_index):	
	pid = str(os.getpid()) + '_'	
	alignment_file = pid + 'align.fasta'
	reference_block_file = pid + 'blockr.fasta'
	assembly_block_file = pid + 'blocka.fasta'
	lagan_cmd = ['perl', LAGAN_DIR + "/lagan.pl", reference_block_file, assembly_block_file, '-mfa']
	alignment_handle = open(alignment_file, 'w')	
	synteny_block_id, reference_instance, assembly_instance = unique_block[block_index]	
	reference_start = int(parse_header(reference_instance.description)['Start'])
	reference_chr_id = parse_header(reference_instance.description)['Seq']	
	contig_id = parse_header(assembly_instance.description)['Seq']
	instance = [reference_instance, assembly_instance]
	file_name = [reference_block_file, assembly_block_file]	
	for index, record in enumerate(instance):
		write_fasta_records([record], file_name[index])		
	with open(os.devnull, "w") as fnull:
		subprocess.call(lagan_cmd, stdout=alignment_handle, stderr=fnull)
	alignment_handle.close()
	return parse_alignment(alignment_file, reference_chr_id, synteny_block_id, contig_id, reference_start)
	
def call_variants(directory, reference_seq, assembly_seq, min_block_size, proc_num):	
	os.chdir(directory)
	block_seq = dict()	
	for record in parse_fasta_file(BLOCKS_FILE):				
		block_id = parse_header(record.description)['Block_id']
		if block_id not in block_seq:
			block_seq[block_id] = []			
		block_seq[block_id].append(record)
	
	pool = multiprocessing.Pool(proc_num)
	unique_block = []	
	for synteny_block_id, instance_list in block_seq.items():		
		if len(instance_list) == 2:		
			reference_instance = find_instance(instance_list, reference_seq.keys(), True)
			assembly_instance = find_instance(instance_list, reference_seq.keys(), False)
			if (not reference_instance is None) and (not assembly_instance is None):
				unique_block.append((synteny_block_id, reference_instance, assembly_instance))											
	
	process_block = functools.partial(process_unique_block, unique_block)
	variant = pool.map(process_block, range(len(unique_block)))	
	variant = [obj for obj in itertools.chain.from_iterable(variant)]	
	coords_file_re = re.compile('blocks_coords[0-9]*.txt')	
	coords_file_list = [coords_file for coords_file in os.listdir('.') if coords_file_re.match(coords_file)]
	blocks_coords = itertools.chain.from_iterable([parse_blocks_coords(coords_file) for coords_file in coords_file_list])
	base_cover = dict()
	for seq_group in (reference_seq, assembly_seq):		
		for seq_id, seq in seq_group.items():
			base_cover[seq_id] = [UNCOVER for _ in seq]
			
	for block in blocks_coords:
		reference = [obj for obj in block if obj[0] in reference_seq]
		if reference and len(reference) < len(block):
			for instance in block:
				size = instance[2] - instance[1]
				base_cover[instance[0]][instance[1]:instance[2]] = [COVER] * size
	
	insertion = []
	os.chdir('..')	
	for seq_id, cover in base_cover.items():
		i = 0
		while i < len(cover):
			if cover[i] == UNCOVER:
				start = i
				while i < len(cover) and cover[i] == UNCOVER:
					i += 1
				end = i
				if end - start > min_block_size:
					if seq_id in reference_seq:					
						common_char = str(reference_seq[seq_id][start - 1]) if start > 0 else ''
						assembly_allele = common_char if common_char else None						
						reference_allele = common_char + str(reference_seq[seq_id][start:end])
						variant.append(Variant(seq_id, start, None, None, reference_allele, assembly_allele,
											reference_allele, assembly_allele, None))
					else:												
						assembly_allele = str(assembly_seq[seq_id][start:end])					
						insertion.append(Variant(None, None, seq_id, start, None, assembly_allele,
												 None, assembly_allele, None))
			else:
				i += 1
	
	return (variant, insertion)

def generate_conventional_output(variant_list, handle):
	for variant in variant_list:
		print >> handle, variant

def write_vcf_header(reference, handle):
	vcf_header = ['##fileformat=VCFv4.1', '##source=Sibelia', '##reference=' + strip_chr_id(reference.id)]
	table_header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
	print >> handle, '\n'.join(vcf_header)
	print >> handle, '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
	print >> handle, '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">'
	print >> handle, '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">'
	print >> handle, '\t'.join(table_header)

def write_variants_vcf(variant_list, handle):	
	for variant in variant_list:
		print >> handle, variant.get_vcf_record()
		
def write_insertions_vcf(variant_list, reference_organism, handle):
	ref_len = str(len(reference_organism.seq))	
	reference_chr = strip_chr_id(reference_organism.id)
	for index, variant in enumerate(variant_list):
		ref_pos = '1'
		ref_allele = reference_organism.seq[0]
		contig = variant.get_contig_id()
		assembly_start = variant.get_assembly_pos() + 1
		assembly_end = assembly_start + len(variant.get_assembly_allele())
		start_alt_allele = ref_allele + '[' + contig + ':' + str(assembly_start) + '['
		end_alt_allele = ']' + contig + ':' + str(assembly_end) + ']' + ref_allele 
		start_bnd = 'bnd_' + str(index * 2)
		end_bnd = 'bnd_' + str(index * 2 + 1)
		info = ';'.join(('IMPRECISE', 'SVTYPE=BND', 'CIPOS=0,' + ref_len))
		start_record = [reference_chr, ref_pos, start_bnd, ref_allele, start_alt_allele, '.', '.', info]
		end_record = [reference_chr, ref_pos, end_bnd, ref_allele, end_alt_allele, '.', '.', info]
		for record in (start_record, end_record):
			print >> handle, '\t'.join(record)

def write_insertions_text(variant_list, handle):
	header = ['SEQ_ID', 'POS', 'FRAGMENT']
	print >> handle, '\t'.join(header)
	for variant in variant_list:
		record = [variant.get_contig_id(), str(variant.get_assembly_pos() + 1), variant.get_assembly_allele()]
		print >> handle, '\t'.join(record)
		
def write_insertions_fasta(variant_list, file_name):
	record = []	
	for variant in variant_list:
		description = 'Seq="' + variant.get_contig_id() + '",Start=' + str(variant.get_assembly_pos() + 1)
		record.append(FastaRecord(seq=variant.get_assembly_allele(), id=description, description=description))
	write_fasta_records(record, file_name)

start = time.time()
parser = argparse.ArgumentParser(description='A tool for comparing two microbial genomes.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('reference', help='A multi-FASTA file with the reference genome')
parser.add_argument('assembly', help='A multi-FASTA file with the assembly genome')
parser.add_argument('-s', '--parameters', help='Parameters set, used for the simplification. \
					Option \"loose\" produces fewer blocks, but they are larger (\"fine\" is opposite).', 
					default='fine')
parser.add_argument('-m', '--minblocksize', help='Minimum size of a synteny block', type=int, default=500)
parser.add_argument('-t', '--tempdir', help='Directory for temporary  files', default='.')
parser.add_argument('-p', '--processcount', help='Number of running processes', type=int, default=1)
parser.add_argument('-i', '--maxiterations', help='Maximum number of iterations during a stage of simplification',
					default=4)
parser.add_argument('-v', '--variant', help='Output file with detected variants', default='variant.vcf')
parser.add_argument('-u', '--unmapped', help='Name of the file to store unmapped insertions in text format', type=str)
args = parser.parse_args()

temp_dir = tempfile.mkdtemp(dir=args.tempdir)
sibelia_cmd = ' '.join([os.path.join(INSTALL_DIR, 'Sibelia'), 					
					args.reference, args.assembly,
					'-q', '--matchrepeats', '--allstages',
					'-m', str(args.minblocksize),
					'-o', temp_dir,
					'-s', args.parameters,
					'-i', str(args.maxiterations)],
					'-r')
print >> sys.stderr, "Calculating synteny blocks..."
os.system(sibelia_cmd)
_, assembly_seq = get_seq(args.assembly)
reference_organism, reference_seq = get_seq(args.reference)
print >> sys.stderr, "Calling variants..."
variant_list, insertion_list = call_variants(temp_dir, reference_seq, assembly_seq,
											 args.minblocksize, args.processcount)
variant_list.sort(key=Variant.get_reference_pos)
vcf_output = open(args.variant, 'w') 
write_vcf_header(reference_organism, vcf_output)
if args.unmapped is None:
	write_insertions_vcf(insertion_list, reference_organism, vcf_output)
else:
	write_insertions_fasta(insertion_list, args.unmapped)
#conventional = open('variant.txt', 'w')
#generate_conventional_output(variant_list, conventional)
#generate_conventional_output(insertion_list, conventional)
#conventional.close()
write_variants_vcf(variant_list, vcf_output)
shutil.rmtree(temp_dir)

