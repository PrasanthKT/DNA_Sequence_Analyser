import sys
import argparse
from Find_Protospacers import find_protospacers


def main():
    parser = argparse.ArgumentParser(description="Process a FASTA file to find protospacers.")
    parser.add_argument("fasta_file", type=str, help="Input FASTA file")
    parser.add_argument("pam_sequence", type=str, help="PAM sequence to search for")
    parser.add_argument("pam_location", type=str, help="Location of PAM sequence (e.g., '3prime' or '5prime')")
    parser.add_argument("protospacer_length", type=int, help="Length of the protospacer")

    args = parser.parse_args()

    find_protospacers(args.fasta_file, args.pam_sequence, args.pam_location, args.protospacer_length)

if __name__ == "__main__":
    main()