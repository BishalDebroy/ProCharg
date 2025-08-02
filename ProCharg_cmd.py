from Bio.SeqUtils.ProtParam import ProteinAnalysis
import os

def read_protein_sequence(filename):
    """Reads protein sequence from a FASTA file, ignoring headers."""
    try:
        with open(filename, 'r') as file:
            sequence_lines = []
            for line in file:
                if line.startswith('>'):
                    continue
                sequence_lines.append(line.strip())
            sequence = ''.join(sequence_lines)
            
            # Validate sequence contains only amino acid characters
            valid_aa = "ACDEFGHIKLMNPQRSTVWY"
            if any(aa.upper() not in valid_aa for aa in sequence):
                raise ValueError("Sequence contains invalid amino acid characters")
            
            return sequence
    except Exception as e:
        print(f"Error reading file: {str(e)}")
        exit(1)

def main():
    input_file = "protein.txt"
    output_file = "protein_output.txt"
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        exit(1)
    
    sequence = read_protein_sequence(input_file)
    
    # Validate sequence length
    if len(sequence) < 3:
        print("Error: Protein sequence is too short for analysis.")
        exit(1)
    
    # Get pH value from user
    try:
        ph = float(input("Enter the pH value for charge calculation (0-14): "))
        if ph < 0 or ph > 14:
            raise ValueError("pH must be between 0 and 14")
    except ValueError as e:
        print(f"Invalid pH value: {str(e)}")
        exit(1)
    
    # Analyze protein
    protein_analysis = ProteinAnalysis(sequence)
    mw = protein_analysis.molecular_weight()
    pI = protein_analysis.isoelectric_point()
    aromaticity = protein_analysis.aromaticity()
    instability_index = protein_analysis.instability_index()
    gravy = protein_analysis.gravy()
    sec_struct = protein_analysis.secondary_structure_fraction()
    charge_at_ph = protein_analysis.charge_at_pH(ph)
    charge_at_pI = protein_analysis.charge_at_pH(pI)
    
    # Create output
    output = f"Protein Sequence Source: {input_file}\n"
    output += f"Sequence Length: {len(sequence)} amino acids\n"
    output += f"Molecular Weight: {mw:.2f} Da\n"
    output += f"Isoelectric Point (pI): {pI:.4f}\n"
    output += f"Aromaticity: {aromaticity:.4f}\n"
    output += f"Instability Index: {instability_index:.2f}\n"
    output += f"GRAVY (Hydrophobicity): {gravy:.4f}\n"
    output += f"Secondary Structure Fractions:\n"
    output += f"  Helix: {sec_struct[0]:.4f}\n"
    output += f"  Turn: {sec_struct[1]:.4f}\n"
    output += f"  Sheet: {sec_struct[2]:.4f}\n"
    output += f"Charge at pH {ph}: {charge_at_ph:.4f}\n"
    output += f"Charge at pI (pH {pI:.4f}): {charge_at_pI:.4f}\n\n"
    
    # Interpretation
    charge_diff = charge_at_ph - charge_at_pI
    if charge_at_ph > 0:
        charge_status = "positively charged"
    elif charge_at_ph < 0:
        charge_status = "negatively charged"
    else:
        charge_status = "neutral"
        
    if abs(charge_diff) < 0.1:
        comparison = "very close to"
    elif charge_diff > 0:
        comparison = "higher than"
    else:
        comparison = "lower than"
    
    output += "Interpretation:\n"
    output += f"At pH {ph}, the protein is {charge_status}.\n"
    output += f"The charge at pH {ph} is {comparison} the charge at its pI.\n"
    output += f"The protein will {'not ' if abs(charge_diff) < 0.5 else ''}behave differently in electrophoresis at this pH compared to its pI.\n"
    
    # Write to file
    with open(output_file, 'w') as out_file:
        out_file.write(output)
    
    print(f"Analysis complete! Results saved to {output_file}")

if __name__ == "__main__":
    main()