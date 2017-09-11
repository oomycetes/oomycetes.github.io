The Python script `generate_list_rxlrs.py` is used to generate lists of candidate RxLR effectors. Four methods were used to identify potential RxLRs:

1. Win method
2. Regex method
3. HMM method
4. BLAST method

Additionally potential RxLRs were screened to detect if they contain a WYL domain using HMMsearch.

Example usage:

	python generate_list_rxlrs.py PHIF.fasta PHIF.signalp hmmsearch_cropped_PHIF.tab hmmsearch_wyl_PHIF.tab blast_PHIF.tab

Input for `generate_list_rxlrs.py` is:

- Fasta file for proteome
- SignalP output for proteome (short format)
- HMMsearch output for cropped.hmm (tabular format)
- HMMsearch output for wyl.hmm search (tabular format)
- BLASTp output for proteome vs reference RxLR effectors (outfmt 6 / tabular)

The output from this script is a tabular list of potential RxLR effectors found in the proteome.
Output columns:

- **Protein ID**
- **Win** (Y/N if the protein had a hit according to Win method)
- **Regex** (Y/N if the protein had a hit according to Regex method)
- **HMM** (Y/N if the protein had a hit according to HMM method)
- **WYL domain** (Y/N if the protein has a W/Y/L domain identified by HMMsearch)
- **Similar to** (List of reference RxLRs that share significant sequence similarity)
- **SignalP HMM Score**
- **SignalP cleavage site**
- **RxLR position**
- **RxLR-EER sequence**
- **Protein length**
- **Protein sequence**

------

The Python script `parse_signal_results.py` is used to parse SignalP results to generate the predicted secretome.

Example usage:

	parse_signalp_results.py signalp_results.txt input.fasta output.fasta

Input for `parse_signal_results.py` is:

- Output file from SignalP (3.0) in short format
- Original fasta sequences corresponding to the SignalP output file (proteomes)
- Ouput fasta filename which stores the mature protein sequences

Proteins with HMM S prob >= 0.9, NN Ymax score >= 0.5 and NN D-score >= 0.5 were considered potentially secreted and were then
submitted to [TMHMM Web Server](http://www.cbs.dtu.dk/services/TMHMM/) to predict transmembrane helices. Proteins that
did not contain transmembrane helicies after the SignalP cleavage site and matched the above criteria were considered
secreted.