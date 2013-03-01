# only contigs of length >= 1000
./fastapl -ge 'length $seq >= 1000' $1

# sort contigs by length (asc)
# ./fastapl  -se 'length $seq1 <=> length $seq2' $1
