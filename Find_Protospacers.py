# Usage: python Find_Protospacers.py <FASTA file> <PAM sequence> <PAM location> <Protospacer length>

import sys
from Bio import SeqIO 
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 

"""git 
This script will take a FASTA Sequence as input and return all the protospacers for a particular PAM sequence, PAM location
and length of the Protospacer.
"""

def find_protospacers(fasta_file, pam_sequence, pam_location, protospacer_length):
    # Open the FASTA file and parse the sequences
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        protospacers = []

        # Find protospacers based on PAM sequence and location of the PAM.
        if pam_location == '3prime':
            for i in range(len(sequence) - protospacer_length + 1):
                protospacer = sequence[i:i + protospacer_length]
                if sequence[i + protospacer_length:i + protospacer_length + len(pam_sequence)] == pam_sequence:
                    protospacers.append(protospacer)
        elif pam_location == '5prime':
            for i in range(len(sequence) - protospacer_length + 1):
                protospacer = sequence[i:i + protospacer_length]
                if sequence[i - len(pam_sequence):i] == pam_sequence:
                    protospacers.append(protospacer)

        # Print or process the protospacers
        for protospacer in protospacers:
            print(f"Found protospacer: {protospacer}")
