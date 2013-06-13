import re
import os
import sys
import time
import argparse
import itertools
import functools
import subprocess
import multiprocessing
from Bio import SeqIO
from Bio import AlignIO

MINIMUM_CONTEXT_SIZE = 30
BLOCKS_FILE = 'blocks_sequences.fasta'
LAGAN_DIR = 'C:/Temp/lagan20'
COVER = 1
UNCOVER = 0
os.environ['LAGAN_DIR'] = LAGAN_DIR

def strip_chr_id(chr_id):
	part = chr_id.split('|')
	if len(part) == 5:
		return part[-2].split('.')[0]
	return chr_id

class Variant(object):
	def __init__(self, reference_chr_id, reference_pos, contig_id, reference_allele,
				 assembly_allele, reference_context, assembly_context, synteny_block_id):
		self._reference_chr_id = reference_chr_id
		self._reference_pos = reference_pos		
		self._contig_id = str(contig_id)
		self._reference_allele = reference_allele.upper()  
		self._assembly_allele = assembly_allele.upper()
		self._reference_context = str(None if reference_context is None else reference_context.upper())
		self._assembly_context = str(None if assembly_context is None else assembly_context.upper())
		self._synteny_block_id = synteny_block_id

	def __str__(self):
		return "\t".join([str(self._reference_pos), self._reference_allele, self._assembly_allele,
						str(self._synteny_block_id), self._contig_id,
						self._reference_context, self._assembly_context])
				
	def get_synteny_block_id(self):
		return self._synteny_block_id	
	
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
		context.append(str(alignment[0][start:segment[1]].seq))
	else:
		context.append('')

	if segment_index + 1 < len(alignment_segment):
		segment = alignment_segment[segment_index + 1]
		end = segment[0] + min(segment[1] - segment[0], MINIMUM_CONTEXT_SIZE)
		context.append(str(alignment[0][segment[0]:end].seq))
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
	alignment = AlignIO.read(open(alingment_file_name), "fasta")
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
	alignment_segment.append([start_position, now_position, last_match])
	alignment_segment[-1][1] = alignment.get_alignment_length()

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
			variant.append(Variant(reference_chr_id, variant_reference_start - shift, contig_id,
								reference_allele, assembly_allele, reference_context, assembly_context,
								synteny_block_id))				
	return variant		   

def get_reference_seq(file_name):
	seq = [(record.id, record.seq) for record in SeqIO.parse(file_name, "fasta")]	
	return (seq[0][0], dict(seq))

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
	lagan_cmd = [LAGAN_DIR + "/lagan.pl", reference_block_file, assembly_block_file, '-mfa']
	alignment_handle = open(alignment_file, 'w')	
	synteny_block_id, reference_instance, assembly_instance = unique_block[block_index]
	reference_start = int(parse_header(reference_instance.description)['Start'])
	reference_chr_id = parse_header(reference_instance.description)['Seq']
	contig_id = parse_header(assembly_instance.description)['Seq']
	instance = [reference_instance, assembly_instance]
	handle = [open(reference_block_file, 'w'), open(assembly_block_file, 'w')]
	for index, record in enumerate(instance):
		SeqIO.write(record, handle[index], 'fasta')
		handle[index].close()
	with open(os.devnull, "w") as fnull:
		subprocess.call(lagan_cmd, stdout=alignment_handle, stderr=fnull)
	alignment_handle.close()
	return parse_alignment(alignment_file, reference_chr_id, synteny_block_id, contig_id, reference_start)
	
def call_variants(directory, reference_seq, min_block_size):	
	os.chdir(directory)
	'''
	block_seq = dict()	
	for record in SeqIO.parse(BLOCKS_FILE, 'fasta'):				
		block_id = parse_header(record.description)['Block_id']
		if block_id not in block_seq:
			block_seq[block_id] = []			
		block_seq[block_id].append(record)		
			
	cpu_count = multiprocessing.cpu_count()
	pool = multiprocessing.Pool(cpu_count)
	unique_block = []	
	for synteny_block_id, instance_list in block_seq.items():		
		if len(instance_list) == 2:		
			reference_instance = find_instance(instance_list, reference_seq_id, True)
			assembly_instance = find_instance(instance_list, reference_seq_id, False)
			if (not reference_instance is None) and (not assembly_instance is None):
				unique_block.append((synteny_block_id, reference_instance, assembly_instance))											
		
	process_block = functools.partial(process_unique_block, unique_block)	
	variant = pool.map(process_block, range(len(unique_block)))	
	variant = [obj for obj in itertools.chain.from_iterable(variant)]'''
	variant = []
	coords_file_re = re.compile('blocks_coords[0-9]*.txt')	
	coords_file_list = [file for file in os.listdir('.') if coords_file_re.match(file)]
	blocks_coords = itertools.chain.from_iterable([parse_blocks_coords(coords_file) for coords_file in coords_file_list])
	base_cover = dict()
	for seq_id, seq in reference_seq.items():
		base_cover[seq_id] = [UNCOVER for _ in seq]
	for block in blocks_coords:
		reference = [obj for obj in block if obj[0] in reference_seq]
		if reference and len(reference) < len(block):
			for instance in reference:
				size = instance[2] - instance[1]
				base_cover[instance[0]][instance[1]:instance[2]] = [COVER] * size
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
					common_char = str(reference_seq[seq_id][start - 1]) if start > 0 else ''
					reference_allele = common_char + str(reference_seq[seq_id][start:end])
					variant.append(Variant(seq_id, start, None, reference_allele, 
										common_char, None, None, None))
			else:
				i += 1
	
	return variant

def generate_conventional_output(variant_list, handle):
	for variant in variant_list:
		print >> handle, variant

def generate_vcf_output(variant_list, reference, handle):
	vcf_header = ['##fileformat=VCFv4.1', '##source=Sibelia', '##reference=' + strip_chr_id(reference)]
	table_header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
	print >> handle, '\n'.join(vcf_header)
	print >> handle, '\t'.join(table_header)
	for variant in variant_list:
		print >> handle, variant.get_vcf_record()

start = time.time()
parser = argparse.ArgumentParser(description='A tool for comparing two microbial genomes.')
parser.add_argument('reference', help='A multi-FASTA file with the reference genome')
parser.add_argument('assembly', help='A multi-FASTA file with the assembly genome')
parser.add_argument('-t', '--tempdir', help='Directory for temporary  files', default='out')
parser.add_argument('-m', '--minblocksize', help='Minimum size of a synteny block', type=int, default=500)


					
'''
args = parser.parse_args()
sibelia_cmd = ' '.join(['Sibelia', '-s fine', '-m', str(args.minblocksize), args.reference, 
						args.assembly, '--comparative', '-o', args.tempdir])
'''

reference_organism, reference_seq = get_reference_seq('H37Rv.fasta')
variant_list = call_variants('.tuberout', reference_seq, 500)
variant_list.sort(key=Variant.get_reference_pos)
generate_conventional_output(variant_list, open('variant.txt', 'w'))
generate_vcf_output(variant_list, reference_organism, open('variant.vcf', 'w'))
print >> sys.stderr, (time.time() - start), "seconds elapsed"