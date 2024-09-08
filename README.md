# Protospacer_Finder

Protospacer_Finder is a Python-based tool designed to identify protospacers for various CRISPR-Cas systems. It serves as a universal protospacer finder for all Cas nucleases with different PAM (Protospacer Adjacent Motif) sequences.

# Overview

The CRISPR-Cas system has revolutionized gene and genome editing. This tool aids researchers by:
Identifying protospacers for specific PAM sequences
Locating protospacers at desired positions
Allowing customization of protospacer lengths
Detecting potential off-targets in the genome
While this tool doesn't rank gRNAs or predict their efficiency, it provides valuable information on off-target presence for each protospacer sequence.

# Features

Universal compatibility with various CRISPR-Cas systems
Customizable PAM sequence and location (3' or 5')
Adjustable protospacer length
Off-target detection in reference genomes

## Installation/Usage: 
Cloning the repository, use the following command:
```
git clone https://github.com/PrasanthKT/Protospacer_Finder.git
cd Protospacer_Finder
```
## Dependencies
Ensure you have the following dependencies installed:
1. Conda 24.7.1
2. Python 3.12.2
3. Biopython 1.84
4. BWA 0.7.18

## Running the Tool 
```bash
conda activate your_environment_name
```
Run the tool using the following command:
```bash
python main.py <fasta_file> <pam_sequence> <pam_location> <protospacer_length> <reference_genome> <output_file>
```
Arguments:
fasta_file: Target gene FASTA file
pam_sequence: PAM sequence of the nuclease
pam_location: Location of the PAM (3prime or 5prime)
protospacer_length: Length of the protospacer
reference_genome: Reference genome for off-target detection
output_file: Name of the output file

## Example Run:
```
python main.py gene.fasta NGG 3prime 20 ref_genome.fasta Output.txt
```
## Output
The script generates a tab-separated output file containing:
Protospacer sequence
PAM sequence
Matching sequence in the genome
Off-target count
```
  Protospacer         PAM Sequence    Matching Sequence       Off-Target Count
GTTCCCGTTACAATCCTTAC    AGG        GTTCCCGTTACAATCCTTACAGG 	      1
ACAATCCTTACAGGATTCTT    AGG        ACAATCCTTACAGGATTCTTAGG 	      1
CTTACAGGATTCTTAGGTTC    AGG        CTTACAGGATTCTTAGGTTCAGG 	      1
```
# Future Developments
We plan to enhance the tool with:

Customizable mismatch allowance for off-target detection
Specification of allowed mismatch positions
Identification of intronic/exonic targets (with GTF/GFF file input)
User-friendly interface for file uploads and genome selection
Customizable protospacer ID assignment 

# Contributing
This is an open-source project, and contributions are welcome. Please feel free to submit pull requests or open issues for bugs and feature requests.

## Contact
In case of any problems, bugs, or feedback, please contact: prasanthkumar.t.22@gmail.com
