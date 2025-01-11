import os
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from collections import Counter
import random
from tqdm import tqdm

def read_fasta(file_path):
    """
    Reads a FASTA file and extracts the DNA sequence.
    :param file_path: Path to the FASTA file
    :return: DNA sequence (str)
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    
    with open(file_path, "r") as file:
        lines = file.readlines()
        sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
        return sequence.upper()

def validate_sequence(sequence):
    """
    Validates a DNA sequence.
    :param sequence: DNA sequence (str)
    :return: Tuple
    """
    sequence = sequence.upper()
    if not sequence:
        return False, "Sequence is empty!"
    
    if not all(base in "ATCG" for base in sequence):
        return False, "Invalid sequence! Only A, T, C, and G are allowed."
    
    if len(sequence) < 3:
        return False, "Sequence is too short for analysis (minimum 3 bases required)."
    
    return True, ""

def calculate_gc_at_content(sequence):
    """
    Calculates GC and AT content of a DNA sequence.
    :param sequence: DNA sequence (str)
    :return: Tuple (gc_content: float, at_content: float)
    """
    gc_content = gc_fraction(sequence) * 100
    at_content = 100 - gc_content
    return gc_content, at_content

def nucleotide_frequencies(sequence):
    """
    Calculates nucleotide frequencies in a DNA sequence.
    :param sequence: DNA sequence (str)
    :return: Dict of nucleotide frequencies
    """
    count = Counter(sequence)
    return {base: count.get(base, 0) for base in "ATCG"}

def translate_dna(sequence):
    """
    Translates a DNA sequence into a protein sequence.
    :param sequence: DNA sequence (str)
    :return: Protein sequence (str)
    """
    trimmed_sequence = sequence[:len(sequence) - (len(sequence) % 3)]  # Trim to multiple of 3
    dna_seq = Seq(trimmed_sequence)
    protein_seq = dna_seq.translate()
    return str(protein_seq)

def find_restriction_sites(sequence, enzymes):
    """
    Finds restriction enzyme cut sites in a DNA sequence.
    :param sequence: DNA sequence (str)
    :param enzymes: Dictionary of enzyme names and recognition sequences
    :return: Dictionary of enzyme names and their cut positions
    """
    sequence = sequence.upper()
    sites = {}
    for enzyme, site in tqdm(enzymes.items(), desc="Analyzing restriction sites"):
        positions = [i for i in range(len(sequence) - len(site) + 1) if sequence[i:i+len(site)] == site]
        if positions:
            sites[enzyme] = positions
    return sites

def simulate_mutations(sequence, mutation_rate, indel_rate=0.0):
    """
    Simulates random mutations in a DNA sequence, including substitutions, insertions, and deletions.
    :param sequence: DNA sequence (str)
    :param mutation_rate: Probability of substitution (0 to 1)
    :param indel_rate: Probability of insertion or deletion (0 to 1)
    :return: Mutated sequence (str)
    """
    sequence = sequence.upper()
    mutated = []
    for base in sequence:
        if random.random() < indel_rate:
            # Insertion
            mutated.append(random.choice("ATCG"))
        elif random.random() < mutation_rate:
            # Substitution
            mutated.append(random.choice("ATCG"))
        else:
            # No mutation
            mutated.append(base)
    # Random deletions
    mutated = [base for base in mutated if random.random() > indel_rate]
    return ''.join(mutated)

def codon_usage(sequence):
    """
    Calculates codon usage in a DNA sequence.
    :param sequence: DNA sequence (str)
    :return: Dict of codons and their counts
    """
    codons = [sequence[i:i+3] for i in range(0, len(sequence) - len(sequence) % 3, 3)]
    return Counter(codons)