#!/usr/bin/env python

# Parses output from SignalP 3 output (short format)
# Proteins with HMM S prob >= 0.9, NN Ymax score >= 0.5 and NN D-score >= 0.5 are considered secreted.
# These were then submitted to TMHMM to predict TM domains. Proteins with a TM domain after the signal peptide
# are removed and not considered secreted

# Example usage: python parse_signalp_results.py signalp_results.txt input.fasta output.fasta
# Mature sequences of proteins that are potentially secreted are written to output.fasta

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys

sigp_file = open(sys.argv[1], 'r')
original_fasta_file = sys.argv[2]
output_fasta_file = sys.argv[3]

secreted_ids = {}

for line in sigp_file:
	# Check if its a comment line
	if line[0] != '#' and not 'ERROR' in line:
		l = line.split()
		prediction = l[-1]
		hmm_score = float(l[-2]) # S prob
		y_max = float(l[4])
		cleavage_site = int(l[5]) # Cleavage site according to Ymax
		d_score = float(l[12])

		if prediction == 'Y' and y_max >= 0.5 and d_score >= 0.5 and hmm_score >= 0.9:
			secreted_ids[l[0]] = cleavage_site - 1

# Writes the MATURE proteins sequences to a FASTA file
output_fasta_file = open(output_fasta_file, 'w')
for record in SeqIO.parse(original_fasta_file, 'fasta'):
	if str(record.id) in secreted_ids:
		new_record = SeqRecord(Seq(str(record.seq)[secreted_ids[str(record.id)]:]), id = str(record.id), name = str(record.id), description = "")
		SeqIO.write(new_record, output_fasta_file, 'fasta')

output_fasta_file.close()
