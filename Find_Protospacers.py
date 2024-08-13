'''
This Script takes the the Target Gene FASTA file as Input
along with the PAM Sequence, PAM LOcation and Length of the Protospacer Sequence. 
Returns the Protospacers in the given FASTA file based on the PAM sequence and location.
'''

from Bio import SeqIO # Import the SeqIO module from the Bio package
import re # Import the re module for regular expressions

def find_protospacers(fasta_file, pam_sequence, pam_location, protospacer_length):
    '''
    The Standard IUPAC Codes for the Nucleotides.
    This makes the users to just give the PAM Sequence of their choice, instead of giving each PAM Sequence Separetely
    to identify the Protospacers.
    '''
    iupac_codes = {
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]',
        'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
    }
    pam_regex = ''.join([iupac_codes.get(base, base) for base in pam_sequence])
    
    # Read the input FASTA file and parse the sequences.
    protospacers = [] #Intialize an empty list to store the Protospacers.
    for record in SeqIO.parse(fasta_file, "fasta"): 
        sequence = str(record.seq)

        '''
        Identify the Protospacers based on the PAM Sequence and the Location of the PAM.
        '''
        if pam_location.lower() == '3prime':
            for i in range(len(sequence) - protospacer_length - len(pam_sequence) + 1):
                protospacer = sequence[i:i + protospacer_length]
                pam_site = sequence[i + protospacer_length:i + protospacer_length + len(pam_sequence)]
                if re.fullmatch(pam_regex, pam_site):
                    protospacers.append((protospacer, pam_site, sequence[i:i + protospacer_length + len(pam_sequence)]))
        elif pam_location.lower() == '5prime':
            for i in range(len(sequence) - protospacer_length + 1):
                protospacer = sequence[i:i + protospacer_length]
                if i >= len(pam_sequence):
                    pam_site = sequence[i - len(pam_sequence):i]
                    if re.fullmatch(pam_regex, pam_site):
                        protospacers.append((protospacer, pam_site, sequence[i - len(pam_sequence):i + protospacer_length]))
        else:
            raise ValueError("Invalid PAM location. Use '3prime' or '5prime'.") # Raise an error if the PAM location is invalid.
        
    return protospacers # Returns the list of Protospacers which is the input for the Off-Target Analysis.