import argparse
import os, sys
lib_path = os.path.abspath('./lib')
sys.path.append(lib_path)
import vcf
from Bio import AlignIO

MINIMUM_CONTEXT_SIZE = 3
LAGAN_PATH = 'C:/Temp/lagan20'

class Variant(object):
	def __init__(self, reference_pos, contig_id, reference_allele, assembly_allele,
				 reference_context, assembly_context, synteny_block_id):
		self._reference_pos = reference_pos
		self._contig_id = contig_id
		self._reference_allele = reference_allele.upper()  
		self._assembly_allele = assembly_allele.upper()
		self._reference_context = None if reference_context is None else reference_context.upper()
		self._assembly_context = None if assembly_context is None else assembly_context.upper()
		self._synteny_block_id = synteny_block_id

	def __str__(self):
		return "\t".join([str(self._reference_pos), self._reference_allele, self._assembly_allele, 
						  self._contig_id, str(self._synteny_block_id),
						  self._reference_context, self._assembly_context])
				
	def get_synteny_block_id(self):
		return self._synteny_block_id	
	
	def get_reference_context(self):
		return self._reference_context
	
	def get_assembly_context(self):
		return self.assembly_context
	
	def get_reference_pos(self):
		return self.reference_pos
	
	def get_contig_id(self):
		return self.contig_id
	
	def get_reference_allele(self):
		return self.reference_allele
	
	def get_assembly_allele(self):
		return self.assembly_allele

def no_gaps(sequence):
	return ''.join([ch for ch in sequence if ch != '-'])

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
		context.append('')
		
	segment = alignment_segment[segment_index]
	reference_context = context[0] + no_gaps(alignment[0][segment[0]:segment[1]]) + context[1]
	assembly_context = context[0] + no_gaps(alignment[1][segment[0]:segment[1]]) + context[1]
	return reference_context, assembly_context

def parse_alignment(alingment_file_name, synteny_block_id, contig_id, reference_start):	
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
			reference_allele = no_gaps(alignment[0][start - shift: end])
			assembly_allele = no_gaps(alignment[1][start - shift: end])
			variant.append(Variant(variant_reference_start - shift, contig_id, reference_allele, assembly_allele,
								   reference_context, assembly_context, synteny_block_id))
				
	return variant		   


parser = argparse.ArgumentParser(description='A tool for comparing two microbial genomes.')
parser.add_argument('reference', help='A multi-FASTA file with the reference genome')
parser.add_argument('assembly', help='A multi-FASTA file with the assembly genome')
parser.add_argument('-t', '--tempdir', help='Directory for temporary  files', default='out')
parser.add_argument('-m', '--minblocksize', help='Minimum size of a synteny block', type=int, default=500)

					
'''
args = parser.parse_args()
cmd = ' '.join(['Sibelia', '-s fine', '-m', str(args.minblocksize), args.reference,
				 args.assembly, '--comparative', '-o', args.tempdir])'''
os.environ['LAGAN_PATH'] = LAGAN_PATH
lagan_cmd = ' '.join([LAGAN_PATH + "/lagan.pl", 'block10_1.fasta block10_2.fasta', '-mfa > out'])
os.system(lagan_cmd)


#os.system(cmd)