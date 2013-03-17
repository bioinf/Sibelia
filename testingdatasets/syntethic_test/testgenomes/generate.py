#This file generates small synthetic test (see Fig. 4 in the paper).
#Generated genomes contain different kinds of repeats.

import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

alphabet = "ACGT"
def randomString(length):
	return Seq(''.join([alphabet[random.randint(0, len(alphabet) - 1)] for _ in xrange(length)]), generic_dna)

def R(seq):
	return seq.reverse_complement()

def S():
	return randomString(500)

def M(seq):
	mutatedSeq = [' '] * len(seq)
	for i in xrange(len(seq)):
		mutatedSeq[i] = seq[i] if random.random() <= 0.97 else alphabet[random.randint(0, len(alphabet) - 1)]
	return Seq(''.join(mutatedSeq), generic_dna)

A = randomString(20000)		#0
I = randomString(20000)		#1
I0 = randomString(20000)	#2
I1 = randomString(20000)	#3
C = randomString(20000)		#4

G1list = [M(I0), S(), M(A), S(), M(I), S(), M(C), S(), M(I), S(), R(M(I0))]
G2list = [M(I1), S(), M(A), S(), R(M(I)), S(), M(C), S(), M(I), S(), M(I1)]
G1 = Seq(''.join([str(x) for x in G1list]), generic_dna)
G2 = Seq(''.join([str(x) for x in G2list]), generic_dna)

r1 = SeqRecord(seq = G1, id = 'genome1')
r2 = SeqRecord(seq = G2, id = 'genome2')

genomeOutHandle = open('genome.fasta', 'w')
for g in [r1, r2]:
	SeqIO.write(g, genomeOutHandle, 'fasta')
genomeOutHandle.close()
