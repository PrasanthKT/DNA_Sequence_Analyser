'''
This Script takes the Protospacers that were Identified in the Target Gene FASTA file as Input
and Aligns them to the Reference Genome to find the same or similar sequences in the Genome (Off-Targets).
This script use BWA Package to Align the Protospacers to the Reference Genome.
'''

import argparse
import subprocess # Import the subprocess module to run command-line programs
import os

def index_reference_genome(reference_genome):
    '''
    Checks if the Index files for the Reference Genome exists.
    If exists, skips the indexing process.
    '''
    index_files = [reference_genome + ext for ext in ['.bwt', '.pac', '.ann', '.amb', '.sa']]
    if all(os.path.exists(file) for file in index_files):
        print("Index files already exist. Skipping indexing.")
        return


    bwa_index_command = f"bwa index {reference_genome}"
    '''
    Runs the BWA Indexing command to Index the Reference Genome.
    '''
    try:
        print(f"Running command: {bwa_index_command}")
        subprocess.run(bwa_index_command, shell=True, check=True)
        print("Indexing completed.")
    except subprocess.CalledProcessError as e:
        print(f"Error during indexing: {e}")
        raise

def run_bwa(reference_genome, protospacer_fasta, output_file):
    # Ensure the output file can be written.
    try:
        with open(output_file, 'w') as f:
            pass
    except IOError as e:
        print(f"Error: Cannot write to output file {output_file}: {e}")
        raise

    
    bwa_command = f"bwa mem {reference_genome} {protospacer_fasta} > {output_file}"
    '''
    Runs the BWA Alignment command to Align the Protospacers to the Reference Genome.
    '''
    try:
        print(f"Running command: {bwa_command}")
        subprocess.run(bwa_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error during BWA alignment: {e}")
        raise

def parse_bwa_output(output_file):
    '''
    Parses the BWA output file to find the off-targets for each protospacer.
    '''
    off_targets = {}
    with open(output_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.split('\t')
            protospacer = fields[0]
            flag = int(fields[1])
            if flag & 4 == 0:  # Check if the read is mapped
                if protospacer not in off_targets:
                    off_targets[protospacer] = 0
                off_targets[protospacer] += 1
    return off_targets

def find_off_targets(reference_genome, protospacer_fasta, output_file):
    '''
    Main function to find off-targets
    '''
    index_reference_genome(reference_genome)
    run_bwa(reference_genome, protospacer_fasta, output_file)
    return parse_bwa_output(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check for off-targets using BWA.")
    parser.add_argument("reference_genome", help="Path to the reference genome FASTA file.")
    parser.add_argument("protospacer_fasta", help="Path to the protospacer sequences FASTA file.")
    parser.add_argument("output_file", help="Path to the output file for BWA alignment results.")
    
    args = parser.parse_args()
    
    off_targets = find_off_targets(args.reference_genome, args.protospacer_fasta, args.output_file)
    for protospacer, count in off_targets.items():
        print(f"Protospacer: {protospacer}, Off-targets: {count}")
