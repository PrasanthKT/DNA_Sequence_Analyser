# DNA Sequence Analyzer

### Overview
The DNA Sequence Analyzer is a Python-based tool for analyzing DNA sequences. This project was my first experience with Python programming, developed as part of my transition from wet lab research to computational bioinformatics.

### Features
1. Calculates GC content, AT content, and nucleotide frequencies.
2. Translates DNA sequences into protein sequences.
3. Identifies restriction enzyme cut sites (e.g., EcoRI, BamHI, HindIII).
4. Simulates mutations (substitutions, insertions, deletions).
5. Analyzes codon usage in DNA sequences.
6. Outputs results in the terminal with an option to save to a .txt file.

### Dependencies
1. Python: Version 3.12
2. Libraries: Biopython, tqdm

Installation:
```
pip install biopython tqdm
```
### Usage

1. Clone the repository

```
git clone https://github.com/<your-username>/dna-sequence-analyzer.git
```

2. Run the Script
```
python main.py
```
You can input a DNA sequence directly or provide a FASTA file.

### Project Files
1. main.py: Handles user interaction and integrates all functions.
2. analysis.py: Contains the core functions for DNA sequence analysis.

### Note
For detailed information about the project, refer to the report uploaded in this repository.
