# Proto_Finder

Proto_Finder is a Python-based tool designed to identify protospacers for various CRISPR-Cas systems, perform off-target analysis. It serves as a universal protospacer and off-target analyzer, supporting a wide range of PAM sequences.

### Overview
The CRISPR-Cas system has revolutionized gene and genome editing. Proto_Finder enhances CRISPR research by offering:
   1.	PAM Flexibility: Support for diverse PAM sequences (e.g., NGG, TTTV, TYCV) and locations (3’ or 5’).
   2.	Protospacer Identification: Extract protospacers of customizable lengths from any input gene.
   3.	Off-Target Analysis: Detect off-target sites in reference genomes, allowing up to 3 mismatches.
   4.	Genomic Annotation: Integrates GTF files to annotate genomic regions (e.g., exon, intron, intergenic).
   5.	Optimized Workflow: Automatically checks and uses existing indexed genomes, avoiding redundant processing.
   6.	Customizable Outputs: Includes protospacer, off-target sequences, genomic positions, and annotations.
This comprehensive set of features allows researchers to efficiently identify and analyze potential CRISPR targets while considering off-target effects, making it a valuable tool for precise gene editing experiments.

### WorkFlow
1. Input the target gene FASTA file, PAM sequence, and reference genome.
2. Identify protospacers based on the given PAM and location.
3. Map protospacers to the reference genome using BWA.
4. Detect and count off-targets (allowing mismatches).
5. Annotate off-target regions using a GTF file.
6. Output results in a tabular format.

### Dependencies
Ensure you have the following dependencies installed:
1. Conda 24.7.1
2. Python 3.12.2
3. Biopython 1.84
4. BWA 0.7.18
5. Python - Lavenshtein 0.26.1

Create the Conda Environment
```bash
conda create --name proto_finder python=3.12.2 -y
conda activate proto_finder
```
Install the dependencies using the following commands
```bash
conda install -c conda-forge biopython=1.84 -y
conda install -c bioconda bwa=0.7.18 -y
pip install python-Levenshtein==0.26.1
```
### Installation/Usage: 
Cloning the repository, use the following command:
```
git clone https://github.com/PrasanthKT/Proto_Finder.git
cd Proto_Finder
```
### Running the Tool 
```bash
conda activate proto_finder
```
Run the tool using the following command:
```bash
python main.py <fasta_file> <pam_sequence> <pam_location> <protospacer_length> <reference_genome> <gtf_file> <output_file>
```
Arguments:
1. fasta_file: Target gene FASTA file
2. pam_sequence: PAM sequence of the nuclease
3. pam_location: Location of the PAM (3prime or 5prime)
4. protospacer_length: Length of the protospacer
5. reference_genome: Reference genome for off-target detection
6. gtf_file: GTF file for annotating genomic regions.
7. output_file: Name of the output file

### Example Run:
```
python main.py TP53.fasta NGG 3prime 20 Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.110.gtf TP53_results.tsv
```
### Output
The output file (.tsv) includes the following columns:
	1.	Protospacer: Identified protospacer sequences.
	2.	PAM Sequence: Associated PAM sequence.
	3.	Off-Target Count: Total number of off-targets detected.
	4.	Off-Target Position: Genomic coordinates of off-target sites.
	5.	Off-Target Sequence: Matched off-target sequences in the genome.
	6.	Region Annotation: Genomic feature (e.g., exon, intron, intergenic).

### Future Developments
1. Support for more complex mismatch patterns.
2. Multi-threading support for faster analysis.
3. Pre-built CRISPR libraries for popular nucleases.

### Contributations
This is an open-source project, and contributions are welcome. Please feel free to submit pull requests or open issues for bugs and feature requests.

### Contact
In case of any problems, bugs, or feedback, please contact: prasanthkumar.t.22@gmail.com
