import argparse
from Find_Protospacers import find_protospacers # Import the function from the Find_Protospacers.py file
from Find_Off_Targets import find_off_targets # Import the function from the Find_Off_Targets.py file
import subprocess # Import the subprocess module to run command-line programs

def main():
    parser = argparse.ArgumentParser(description="Find protospacers in a given FASTA file based on PAM sequence and location.")
    parser.add_argument("fasta_file", help="Path to the input FASTA file.")
    parser.add_argument("pam_sequence", help="PAM sequence (e.g., NGG, GTTTR).")
    parser.add_argument("pam_location", help="Location of the PAM sequence relative to the protospacer (3prime or 5prime).")
    parser.add_argument("protospacer_length", type=int, help="Length of the protospacer sequence.")
    parser.add_argument("reference_genome", help="Path to the reference genome FASTA file.")
    parser.add_argument("output_file", help="Path to the output file for off-target analysis.")

    args = parser.parse_args()

    # Find the Protospacers in the input Target_Gene FASTA_file.
    protospacers = find_protospacers(args.fasta_file, args.pam_sequence, args.pam_location, args.protospacer_length)

    # Write all the Protospacers into a temperoray FASTA_file for the off-target analysis/align the Protospacers to the reference genome. 
    with open("protospacers.fasta", "w") as f:
        for i, (protospacer, pam_sequence, matching_sequence) in enumerate(protospacers):
            f.write(f">protospacer_{i}\n{protospacer}\n")

    # Index the reference genome and save the files.
    subprocess.run(["bwa", "index", args.reference_genome])

    # Align all the protospacers to the reference genome.
    result = subprocess.run(["bwa", "mem", args.reference_genome, "protospacers.fasta"], capture_output=True, text=True)
    
    # Find the off-targets for each protospacer Sequence and Count the Off_targets for each Indentified Protospacer.
    off_target_counts = {}
    for line in result.stdout.split('\n'):
        if line.startswith('@'):
            continue
        fields = line.split('\t')
        if len(fields) > 9:
            protospacer_seq = fields[9]
            if fields[2] == '*':  # Unmapped
                off_target_counts[protospacer_seq] = off_target_counts.get(protospacer_seq, 0) + 1

    # Write the Output to a file.
    with open(args.output_file, "w") as f:
        f.write("Protospacer\tPAM Sequence\tMatching Sequence\tOff-Target Count\n") # The Output file Header
        for (protospacer, pam_sequence, matching_sequence) in protospacers:
            count = off_target_counts.get(protospacer, 0)
            f.write(f"{protospacer}\t{pam_sequence}\t{matching_sequence}\t{count}\n")

if __name__ == "__main__":
    main()