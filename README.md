# Protospacer_Finder
<br>
<img width="361" alt="CRISPR" src="https://github.com/user-attachments/assets/2dce1e4d-3965-4eb2-9e30-85801cb3325b">

The CRISPR-Cas system is a groundbreaking advancement in gene and genome editing, guided by guide RNAs (gRNAs) that search for specific sequences known as Protospacer Adjacent Motif (PAM). These sequences enable gRNAs to locate and edit targets precisley within the genome. Each CRISPR-Cas system has its own unique PAM sequence to identify and edit target positions. For instance, PAMs can be present either on the 3' end for Cas9 (NGG) or the 5' end for Cas12a (TTTV).

<br>
<img width="604" alt="Protospacer" src="https://github.com/user-attachments/assets/6b4fe841-e859-49e3-ae88-b1f1351df074">

This project introduces a tool developed using Python that identifies protospacers for various CRISPR-Cas systems, essentially serving as a universal Protospacer finder for all Cas Nucleases which have different PAMs. The primary objective of a CRISPR-Cas system is to accurately edit specific targets while minimizing off-target effects, which are critical in gene/genome editing. This tool can be used to locate protospacers for specific PAM sequences, at desired locations, and with defined protospacer lengths. Additionally, it can identify off-targets present in the genome, enhancing the precision of CRISPR-Cas gene editing.

<br>
<img width="922" alt="Workflow" src="https://github.com/user-attachments/assets/436c5b81-f5c1-4f34-a631-e1c57082c0c5">

The results from this tool can be used to design gRNAs. Although this tool does not have features to rank gRNAs or predict their efficiency, it provides the number of off-targets present in the genome for each protospacer Sequence. These results can also be used to determine whether target regions are locate in intronic or exonic regions if you have the GTF/GFF format file of the genome.

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
python main.py 'fasta_file' 'pam_sequence' 'pam_location' 'protospacer_length' 'reference_genome' 'output_file'
```
Arguments:
1. fasta_file - Target Gene Fasta file
2. pam_sequence - PAM Sequence of the Nuclease.
3. pam_location - Location of the PAM (either 3' or 5')
4. protospacer_length - Length of the Protospacer (ex: 20bp for Cas9 24bp for Cas12a)
5. reference_genome - Reference genome of the organism to find the off targets for each protospacer
6. output_file - Name of the Output file for the script to write.

## Example Run:
```
python main.py gene.fasta NGG 3prime 20 ref_genome.fasta Output.txt
```
## Output
The output of the script consists of all the protospacers for the target gene, the PAM sequence of the protospacer, the matching sequence in the genome, and the off-target count.
```
  Protospacer         PAM Sequence    Matching Sequence       Off-Target Count
GTTCCCGTTACAATCCTTAC    AGG        GTTCCCGTTACAATCCTTACAGG 	      1
ACAATCCTTACAGGATTCTT    AGG        ACAATCCTTACAGGATTCTTAGG 	      1
CTTACAGGATTCTTAGGTTC    AGG        CTTACAGGATTCTTAGGTTCAGG 	      1
```
## Roadmap
Future versions of the tool could include the following features:
1. Allowing input for the number of mismatches for finding off-targets.
2. Specifying positions where mismatches are allowed.
3. Identifying whether targets are present in intronic or exonic regions if a GTF/GFF file of the organism is provided.
4. A clear user interface can be developed in future versions, allowing users to simply upload their gene FASTA file and select the reference genome file that are already indexed using BWA and the associated GTT/GTF file to find the protospacers and their regions.
5. Unique protospacer_id can be assigned according to the user specification and the next ids can be followed accordingly for labelling them for their experiments 

## Open Tool
This is open-source, and contributions are welcome.

## Contact
In case of any problems, bugs, or feedback, please contact: prasanthkumar.t.22@gmail.com
