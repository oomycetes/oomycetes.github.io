The Python script `generate_list_rxlrs.py` is used to generate lists of candidate RxLR effectors. Four methods were used to identify potential RxLRs:

1. Win method
2. Regex method
3. HMM method
4. BLAST method

Additionally potential RxLRs were screened to detect if they contain a WYL domain using HMMsearch.

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

