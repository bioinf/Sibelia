import argparse
import os, sys
lib_path = os.path.abspath('./lib')
sys.path.append(lib_path)
import vcf
from Bio import AlignIO

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
                          self._reference_context, self._assembly_context, 
                          self._contig_id, str(self._synteny_block_id)])
                
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

def get_context(alignment_segment, segment_index):
    return ('', '')

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
            alignment_segment.append((start_position, now_position, last_match))
            last_match = now_match
            start_position = now_position

    variant = []
    alignment_segment.append((start_position, alignment.get_alignment_length(), last_match))
    for segment_index, segment in enumerate(alignment_segment):
        start, end, match = segment
        if match == False:
            variant_reference_start = reference_start + start                 
            reference_context, assembly_context = get_context(alignment_segment, segment_index)
            if end - start == 1 and alignment[0][start] != '-' and alignment[1][start] != '-':
                variant.append(Variant(variant_reference_start, contig_id, alignment[0][start], alignment[1][start],
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

variant = parse_alignment("block115.fasta", 0, "contig_1", 0)
for v in variant:
    print v



#os.system(cmd)