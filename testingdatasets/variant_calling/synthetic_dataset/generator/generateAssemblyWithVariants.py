import sys
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

element = ['A', 'C', 'G', 'T']

def cutGenomeIntoContigs(sequence, pattern):
	ret = []
	begin = 0
	contigCount = 0
	for record in pattern:
		end = min(begin + len(record.seq), len(sequence))
		cseq = sequence[begin:end]
		if random.randint(0, 1) == 0:
			cseq.reverse_complement()
		ret.append(SeqRecord(seq = cseq, id = 'contig_' + str(contigCount)))
		contigCount += 1
		begin = end
		if begin >= len(sequence):
			break

	if begin < len(sequence):
		ret.append(SeqRecord(seq = sequence[begin:], id = 'contig_' + str(contigCount)))
	return ret
	
def generateMismatch(sequence, pos, length, refPos, diffList, _):
	refAllele = sequence[pos:pos + length]
	for i in xrange(length):
		j = pos + i
		if j < len(sequence):
			toReplace = [ch for ch in element if ch != sequence[j]]
			sequence[j] = random.choice(toReplace)
	altAllele = sequence[pos:pos + length]
	diffList.append((refPos[pos], refAllele, altAllele))

def generateInsertion(sequence, pos, length, refPos, diffList, _):
	refAllele = '.'
	varPos = refPos[pos]
	for i in xrange(length):
		j = pos + i
		sequence.insert(pos, random.choice(element))
		refPos.insert(pos, varPos)
	altAllele = sequence[pos:pos + length]
	diffList.append((varPos, refAllele, altAllele))

def generateDeletion(sequence, pos, length, refPos, diffList, _):
	varPos = refPos[pos]
	refAllele = sequence[pos:pos + length]	
	del sequence[pos:pos + length]
	del refPos[pos:pos + length]
	altAllele = '.'
	diffList.append((varPos, refAllele, altAllele))	

def generateReversal(sequence, pos, length, refPos, _, rearrList):
	begin = pos
	end = min(begin + length, len(sequence))
	rearrList.append(('Reversal', refPos[begin], refPos[end]))
	revSeq = sequence[begin:end]
	revSeq.reverse_complement()
	sequence[begin:end] = revSeq
	revPos = refPos[begin:end]
	revPos.reverse()
	refPos[begin:end] = revPos

def generateTranslocation(sequence, pos, length, refPos, _, rearrList):
	posFrag = refPos[pos:pos + length]
	seqFrag = sequence[pos:pos + length]
	del refPos[pos:pos + length]
	del sequence[pos:pos + length]
	nextPos = random.randint(0, len(sequence))
	rearrList.append(('Translocation', posFrag[0], posFrag[-1], refPos[nextPos]))
	for i in xrange(length):
		sequence.insert(nextPos + i, seqFrag[i])
		refPos.insert(nextPos + i, refPos[i])
		
def generateVariants(genome, pattern, nSNV, nSmallIndel, nLargeIndel, nRearr, genomeOutFile, diffOutFile, rearrOutFile):
	diffList = []
	rearrList = []
	sequence = MutableSeq(str(genome.seq), generic_dna)
	refPos = [i + 1 for i in xrange(len(sequence))]
	variationSize = [1, 100, 2000]
	variationCount = [nSNV, nSmallIndel, nLargeIndel]
	for i in xrange(len(variationSize)):
		variationType = [generateInsertion, generateDeletion]
		if variationSize[i] == 1:
			variationType.append(generateMismatch)
		for _ in xrange(variationCount[i]):
			position = random.randint(0, len(sequence) - variationSize[i])
			variation = random.choice(variationType)
			variation(sequence, position, variationSize[i], refPos, diffList, rearrList)
	rearr = [generateTranslocation, generateReversal]
	rearrSize = 15000
	for variation in rearr:
		for _ in xrange(nRearr):
			position = random.randint(0, len(sequence) - rearrSize)
			variation(sequence, position, rearrSize, refPos, diffList, rearrList)			
	contig = cutGenomeIntoContigs(sequence, pattern)
	diffOutHandle = open(diffOutFile, 'w')
	genomeOutHandle = open(genomeOutFile, 'w')
	rearrOutHandle = open(rearrOutFile, 'w')
	for record in contig:
		SeqIO.write(record, genomeOutHandle, 'fasta')
	diffList.sort()
	for record in diffList:
		print >> diffOutHandle, '\t'.join([str(x) for x in record])
	for record in rearrList:
		print >> rearrOutHandle, '\t'.join([str(x) for x in record])
	for handle in [genomeOutHandle, rearrOutHandle, diffOutHandle]:
		handle.close()

if len(sys.argv) != 6:
	print 'Usage: generateAssemblyWithVariants.py <genome FASTA file> <pattern FASTA file> <genome output file> <diff output file> <rearrangaments output file>'
	exit(1)

genome = SeqIO.parse(sys.argv[1], "fasta").next()
pattern = [record for record in SeqIO.parse(sys.argv[2], "fasta")]
generateVariants(genome, pattern, 100, 20, 10, 2, sys.argv[3], sys.argv[4], sys.argv[5])