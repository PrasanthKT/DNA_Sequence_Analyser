# Proto_Finder

Proto_Finder is a Python-based tool designed to identify protospacers for various CRISPR-Cas systems. It serves as a universal protospacer finder for all Cas nucleases with different PAM (Protospacer Adjacent Motif) sequences.

### Overview
The CRISPR-Cas system has revolutionized gene and genome editing. Proto_Finder enhances CRISPR research by offering:
1. Universal compatibility with various CRISPR-Cas systems
2. Identification of protospacers for specific PAM sequences
3. Customizable PAM sequence and location (3' or 5')
4. Adjustable protospacer length
5. Off-target detection in reference genomes

This comprehensive set of features allows researchers to efficiently identify and analyze potential CRISPR targets while considering off-target effects, making it a valuable tool for precise gene editing experiments.

### Installation/Usage: 
Cloning the repository, use the following command:
```
git clone https://github.com/PrasanthKT/Proto_Finder.git
cd Proto_Finder
```
### Dependencies
Ensure you have the following dependencies installed:
1. Conda 24.7.1
2. Python 3.12.2
3. Biopython 1.84
4. BWA 0.7.18

### Running the Tool 
```bash
conda activate your_environment_name
```
Run the tool using the following command:
```bash
python main.py <fasta_file> <pam_sequence> <pam_location> <protospacer_length> <reference_genome> <output_file>
```
Arguments:
1. fasta_file: Target gene FASTA file
2. pam_sequence: PAM sequence of the nuclease
3. pam_location: Location of the PAM (3prime or 5prime)
4. protospacer_length: Length of the protospacer
5. reference_genome: Reference genome for off-target detection
6. output_file: Name of the output file

### Example Run:
```
python main.py gene.fasta NGG 3prime 20 ref_genome.fasta Output.txt
```
### Output
The script generates the output file containing:
1. Protospacer sequence
2. PAM sequence
3. Matching sequence in the genome
4. Off-target count
```
  Protospacer         PAM Sequence    Matching Sequence       Off-Target Count
GTTCCCGTTACAATCCTTAC    AGG        GTTCCCGTTACAATCCTTACAGG 	      1
ACAATCCTTACAGGATTCTT    AGG        ACAATCCTTACAGGATTCTTAGG 	      1
CTTACAGGATTCTTAGGTTC    AGG        CTTACAGGATTCTTAGGTTCAGG 	      1
```
### Future Developments
1. Customizable mismatch allowance for off-target detection
2. Specification of allowed mismatch positions
3. Identification of intronic/exonic targets (with GTF/GFF file input)
4. User-friendly interface for file uploads and genome selection
5. Customizable protospacer ID assignment 

### Contributations
This is an open-source project, and contributions are welcome. Please feel free to submit pull requests or open issues for bugs and feature requests.

### Contact
In case of any problems, bugs, or feedback, please contact: prasanthkumar.t.22@gmail.com
