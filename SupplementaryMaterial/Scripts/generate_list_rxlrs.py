#!/usr/bin/env python

# Reads in Fasta files, SignalP output, HMMsearch output and BLAST output to generate tabular list of putative RxLRs
# Uses three methods described by Haas et al. (2009): Win method, Regex method and HMM method.
# Also uses HMMsearch to detect presence of WYL domains and uses BLASTp to identify proteins
# with significant sequence similarity to previously definied RxLRs (from Win et al., 2009)

# Input: protein fasta file, signalP output, cropped.hmm HMMsearch tabular, wyl.hmm HMMsearch tabular, BLAST tabular output

from Bio import SeqIO
import re, sys

def main():
	fasta_file = sys.argv[1]
	signalp_file = sys.argv[2]
	hmm_file = sys.argv[3]
	wyl_hmm_file = sys.argv[4]
	blast_file = sys.argv[5]

	regex_win = re.compile("^.{29}.{0,26}R.LR")
	regex_eer = re.compile("^.{10,40}.{1,96}R.LR.{1,40}[ED][ED][KR]")
	regex_rxlr_eer = re.compile("^.{10,40}.{1,96}R.LR")

	win_rlxrs = {}
	regex_rxlrs = {}
	hmm_rxlrs = {}

	# Parse SignalP output to get list of secreted proteins
	secreted_proteins = dict()

	with open(signalp_file, 'r') as f:
		for line in f:
			if line[0] != '#' and not 'ERROR' in line:
				l = line.split()
				protein_id = l[0]
				prediction = l[-1]
				hmm_score = float(l[-2])
				cleavage_site = int(l[5]) - 1

				if prediction == 'Y' and hmm_score >= 0.9:
					secreted_proteins[protein_id] = [hmm_score, cleavage_site]

	for record in SeqIO.parse(fasta_file, 'fasta'):
		protein_sequence = str(record.seq)
		protein_id = str(record.id)

		if protein_id in secreted_proteins:
			cleavage_site = secreted_proteins[protein_id][1]
			sigp_score = secreted_proteins[protein_id][0]

			# Win method, looks for a signal peptide < 30, RxLR between 30 and 60 
			hit = motif_search(protein_sequence, regex_win)
			if hit:
				rxlr_position = hit[1] - 3 # 0 indexed, get index of start of RxLR
				rxlr_seq = protein_sequence[rxlr_position - 1: rxlr_position + 3]
				if rxlr_position > cleavage_site and cleavage_site <= 30:
					win_rlxrs[protein_id] = [rxlr_position, rxlr_seq, len(protein_sequence), protein_sequence]

			# Regex method, looks for RxLR and EER-like motif
			hit = motif_search(protein_sequence, regex_eer)
			if hit:
				eer_position = hit[1] - 2 # 0 indexed, get index of start of EER
				eer_seq = protein_sequence[eer_position - 1: eer_position + 2]

				# Have the location of the EER, also get RxLR pos
				rxlr_position = motif_search(protein_sequence[0: eer_position], regex_rxlr_eer)[1] - 3
				rxlr_eer_seq = protein_sequence[rxlr_position - 1: eer_position + 2]

				if rxlr_position > cleavage_site and cleavage_site >= 10 and cleavage_site <= 40:
					regex_rxlrs[protein_id] = [rxlr_position, rxlr_eer_seq, len(protein_sequence), protein_sequence]

	# HMM method, bit scores > 0 are hits
	with open(hmm_file, 'r') as f:
		for line in f:
			if line[0] != '#':
				l = line.split()
				hit_id = l[0]
				hit_bitscore = float(l[5])

				if hit_bitscore > 0 and hit_id in secreted_proteins:
					hmm_rxlrs[hit_id] = [hit_bitscore]

	# Get protein length and sequence for HMM hits
	for record in SeqIO.parse(fasta_file, 'fasta'):
		if str(record.id) in hmm_rxlrs:
			hmm_rxlrs[str(record.id)] += [len(str(record.seq)), str(record.seq)]


	# Get hits to wyl.hmm, bit scores > 0 are included
	wyl_hits = {}
	with open(wyl_hmm_file, 'r') as f:
		for line in f:
			if line[0] != '#':
				l = line.split()
				hit_id = l[0]
				hit_bitscore = float(l[5])

				if hit_bitscore > 0 and hit_id in secreted_proteins:
					wyl_hits[hit_id] = hit_bitscore

	# Also identify homologs to reference RxLRs (from Haas et al. 2009) via Blast
	secreted_blast_hits = {}
	secreted_blast_hits_seqs = {}
	with open(blast_file, 'r') as f:
		for line in f:
			l = line.split()
			protein = l[1]

			if protein in secreted_proteins:
				hit = l[0]
				if protein in secreted_blast_hits:
					secreted_blast_hits[protein].add(hit)
				else:
					secreted_blast_hits[protein] = set([hit])

	# Get the sequences for proteins that have a hit to a reference effector (inefficient)
	for record in SeqIO.parse(fasta_file, 'fasta'):
		if str(record.id) in secreted_blast_hits:
			secreted_blast_hits_seqs[str(record.id)] = [len(str(record.seq)), str(record.seq)]

	print "Protein ID\tWin\tRegex\tHMM\tWYL domain\tSimilar To\tSignalP HMM Score\t Signal Cleavage Site\t RxLR Position\tRxLR-EER Sequence\tProtein Length\tProtein Sequence"
	for protein in secreted_proteins:
		if protein in win_rlxrs:
			win = True
		else:
			win = False
		if protein in regex_rxlrs:
			regex = True
		else:
			regex = False
		if protein in hmm_rxlrs:
			hmm = True
		else:
			hmm = False
		if protein in wyl_hits:
			wyl = True
		else:
			wyl = False
		if protein in secreted_blast_hits:
			homolog = True
		else:
			homolog = False

		# Format info about protein for output to tabular file
		line = ""

		if win or regex or hmm or homolog:
			line += protein + "\t"
			if win:
				line += "Y\t"
			else:
				line += "N\t"

			if regex:
				line += "Y\t"
			else:
				line += "N\t"

			if hmm:
				line += "Y\t"
			else:
				line += "N\t"

			if wyl:
				line += "Y\t"
			else:
				line += "N\t"
			if homolog:
				line += ";".join(sorted(secreted_blast_hits[protein])) + "\t"
			else:
				line += "N\t"

			line += str(secreted_proteins[protein][0]) + "\t"
			line += str(secreted_proteins[protein][1]) + "\t"

			if regex:
				regex_info = regex_rxlrs[protein]
				line += "\t".join(map(str, regex_info))
			elif win:
				win_info = win_rlxrs[protein]
				line += "\t".join(map(str, win_info))
			elif hmm:
				hmm_info = hmm_rxlrs[protein]
				line += "\t".join(map(str, ['-', '-'] + hmm_info[1:]))
			elif homolog:
				homolog_info = secreted_blast_hits_seqs[protein]
				line += "\t".join(map(str, ['-', '-'] + homolog_info))

			print line

	print ""
	print "Win:", len(win_rlxrs)
	print "Regex:", len(regex_rxlrs)
	print "HMM:", len(hmm_rxlrs)
	
	print "Both HMM and Regex:", len(set([x for x in hmm_rxlrs]).intersection(set([y for y in regex_rxlrs])))
	print "Either HMM or Regex:", len(set([x for x in hmm_rxlrs]).union(set([y for y in regex_rxlrs])))

	print "Similar to reference effectors (BLAST):", len(secreted_blast_hits)

	win_regex = (set([x for x in win_rlxrs]).union(set([y for y in regex_rxlrs])))
	total_rxlrs = (set([x for x in hmm_rxlrs]).union(win_regex))
	total_rxlrs = (set([x for x in secreted_blast_hits]).union(total_rxlrs))

	print "WYL hmm hits:", len(wyl_hits)
	print "All hits: " + str((len(total_rxlrs)))

	print ""

def motif_search(sequence, regex):
	m = re.search(regex, sequence)
	if m:
		return m.start(), m.end()
	else:
		return False

if __name__ == '__main__':
	main()
