import logging
import analysis

# Configure logging
logging.basicConfig(filename="dna_analyzer.log", level=logging.INFO, 
                    format="%(asctime)s - %(levelname)s - %(message)s")

def save_results(file_name, content):
    """
    Saves the analysis results to a text file.
    :param file_name: Name of the output file
    :param content: The content to save
    """
    with open(file_name, "w") as f:
        f.write(content)

def main():
    print("Welcome to the DNA Sequence Analyzer!")
    
    # Choose input method
    input_method = input("Do you want to input a FASTA file (f) or enter a sequence manually (m)? ").lower()
    if input_method == "f":
        file_path = input("Enter the path to the FASTA file: ").strip()
        try:
            sequence = analysis.read_fasta(file_path)
            print(f"Sequence loaded from file. Length: {len(sequence)}")
        except FileNotFoundError as e:
            print(e)
            return
    elif input_method == "m":
        sequence = input("Enter your DNA sequence: ").upper()
    else:
        print("Invalid input method. Please choose 'f' or 'm'.")
        return
    
    # Validate sequence
    is_valid, message = analysis.validate_sequence(sequence)
    if not is_valid:
        print(f"Error: {message}")
        return
    
    print(f"Sequence is valid! Length: {len(sequence)}")
    logging.info(f"Sequence length: {len(sequence)}")
    
    # GC and AT content
    gc_content, at_content = analysis.calculate_gc_at_content(sequence)
    print("\n--- Nucleotide Content ---")
    print(f"GC Content: {gc_content:.2f}%")
    print(f"AT Content: {at_content:.2f}%")
    logging.info(f"GC Content: {gc_content:.2f}%, AT Content: {at_content:.2f}%")
    
    # Nucleotide frequencies
    frequencies = analysis.nucleotide_frequencies(sequence)
    print("\n--- Nucleotide Frequencies ---")
    for base, count in frequencies.items():
        print(f"{base}: {count}")
    logging.info(f"Nucleotide Frequencies: {frequencies}")
    
    # DNA to protein translation
    protein_seq = analysis.translate_dna(sequence)
    print("\n--- Translated Protein Sequence ---")
    print(protein_seq)
    logging.info(f"Protein Sequence: {protein_seq}")
    
    # Restriction enzyme prediction
    restriction_enzymes = {
        "EcoRI": "GAATTC",
        "BamHI": "GGATCC",
        "HindIII": "AAGCTT"
    }
    sites = analysis.find_restriction_sites(sequence, restriction_enzymes)
    print("\n--- Restriction Enzyme Sites ---")
    if sites:
        for enzyme, positions in sites.items():
            print(f"{enzyme} ({restriction_enzymes[enzyme]}): Found at positions {positions}")
    else:
        print("No restriction sites found.")
    logging.info(f"Restriction Sites: {sites}")
    
    # Codon Usage Analysis
    codon_counts = analysis.codon_usage(sequence)
    print("\n--- Codon Usage ---")
    for codon, count in codon_counts.items():
        print(f"{codon}: {count}")
    logging.info(f"Codon Usage: {codon_counts}")
    
    # Mutation simulation
    mutation_rate = float(input("\nEnter mutation rate (0 to 1, e.g., 0.05 for 5%): "))
    indel_rate = float(input("Enter insertion/deletion rate (0 to 1, e.g., 0.01 for 1%): "))
    mutated_sequence = analysis.simulate_mutations(sequence, mutation_rate, indel_rate)
    print("\n--- Mutated DNA Sequence ---")
    print(mutated_sequence)
    
    # Translate mutated sequence to protein
    mutated_protein_seq = analysis.translate_dna(mutated_sequence)
    print("\n--- Mutated Protein Sequence ---")
    print(mutated_protein_seq)
    logging.info(f"Mutated Sequence: {mutated_sequence}")
    logging.info(f"Mutated Protein Sequence: {mutated_protein_seq}")
    
    # Save results to file
    save_to_file = input("\nDo you want to save the results to a file? (y/n): ").lower()
    if save_to_file == "y":
        output_content = f"""
        DNA Sequence Analyzer Results
        =============================

        Sequence Length: {len(sequence)}
        GC Content: {gc_content:.2f}%
        AT Content: {at_content:.2f}%
        
        Nucleotide Frequencies:
        {frequencies}
        
        Translated Protein Sequence:
        {protein_seq}
        
        Restriction Enzyme Sites:
        {sites if sites else "No restriction sites found."}
        
        Codon Usage:
        {codon_counts}
        
        Mutated DNA Sequence:
        {mutated_sequence}
        
        Mutated Protein Sequence:
        {mutated_protein_seq}
        """
        save_results("dna_analysis_results.txt", output_content)
        print("\nResults saved to 'dna_analysis_results.txt'.")

if __name__ == "__main__":
    main()