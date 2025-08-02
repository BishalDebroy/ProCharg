from Bio.SeqUtils.ProtParam import ProteinAnalysis
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os

class ProteinAnalyzerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("ProCharg: Protein Charge Analyzer")
        self.root.geometry("800x600")
        self.root.configure(bg='#f0f8ff')
        
        # Configure styles
        self.style = ttk.Style()
        self.style.configure('TFrame', background='#f0f8ff')
        self.style.configure('TLabel', background='#f0f8ff', font=('Arial', 10))
        self.style.configure('Header.TLabel', background='#2c3e50', foreground='white', font=('Arial', 14, 'bold'))
        self.style.configure('TButton', font=('Arial', 10))
        self.style.configure('TEntry', font=('Arial', 10))
        
        # Create main frames
        header_frame = ttk.Frame(root)
        header_frame.pack(fill=tk.X, padx=10, pady=10)
        
        input_frame = ttk.Frame(root)
        input_frame.pack(fill=tk.X, padx=20, pady=10)
        
        button_frame = ttk.Frame(root)
        button_frame.pack(fill=tk.X, padx=20, pady=10)
        
        output_frame = ttk.Frame(root)
        output_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)
        
        # Header
        header_label = ttk.Label(header_frame, text="Protein Charge Analyzer", style='Header.TLabel')
        header_label.pack(fill=tk.X, ipady=10)
        
        # Input fields
        ttk.Label(input_frame, text="Protein FASTA File:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.file_path = tk.StringVar()
        file_entry = ttk.Entry(input_frame, textvariable=self.file_path, width=50)
        file_entry.grid(row=0, column=1, padx=5)
        ttk.Button(input_frame, text="Browse", command=self.browse_file).grid(row=0, column=2, padx=5)
        
        ttk.Label(input_frame, text="pH Value:").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.ph_value = tk.StringVar(value="7.0")
        ph_entry = ttk.Entry(input_frame, textvariable=self.ph_value, width=10)
        ph_entry.grid(row=1, column=1, sticky=tk.W, padx=5)
        
        # Buttons
        ttk.Button(button_frame, text="Analyze Protein", command=self.analyze_protein).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Clear Results", command=self.clear_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Save Report", command=self.save_report).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Exit", command=root.quit).pack(side=tk.RIGHT, padx=5)
        
        # Output area
        self.output_text = tk.Text(output_frame, wrap=tk.WORD, font=('Courier New', 10), bg='white', padx=10, pady=10)
        self.output_text.pack(fill=tk.BOTH, expand=True)
        
        # Status bar
        self.status_var = tk.StringVar(value="Ready")
        status_bar = ttk.Label(root, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W)
        status_bar.pack(side=tk.BOTTOM, fill=tk.X)
        
        # Initialize results
        self.results = {}
    
    def browse_file(self):
        file_path = filedialog.askopenfilename(
            title="Select Protein FASTA File",
            filetypes=[("FASTA files", "*.fasta;*.fa;*.txt"), ("All files", "*.*")]
        )
        if file_path:
            self.file_path.set(file_path)
            self.status_var.set(f"Selected file: {os.path.basename(file_path)}")
    
    def read_protein_sequence(self, filename):
        try:
            with open(filename, 'r') as file:
                sequence_lines = []
                for line in file:
                    if line.startswith('>'):
                        continue
                    sequence_lines.append(line.strip())
                sequence = ''.join(sequence_lines)
                
                # Validate sequence
                valid_aa = "ACDEFGHIKLMNPQRSTVWY"
                if any(aa.upper() not in valid_aa for aa in sequence):
                    raise ValueError("Sequence contains invalid amino acid characters")
                
                return sequence
        except Exception as e:
            messagebox.showerror("File Error", f"Error reading file: {str(e)}")
            return None
    
    def analyze_protein(self):
        # Clear previous results
        self.output_text.delete(1.0, tk.END)
        
        # Get inputs
        file_path = self.file_path.get()
        ph_str = self.ph_value.get()
        
        # Validate inputs
        if not file_path:
            messagebox.showerror("Input Error", "Please select a protein FASTA file")
            return
        
        try:
            ph = float(ph_str)
            if ph < 0 or ph > 14:
                raise ValueError("pH must be between 0 and 14")
        except ValueError:
            messagebox.showerror("Input Error", "Please enter a valid pH value (0-14)")
            return
        
        # Read sequence
        self.status_var.set("Reading protein sequence...")
        sequence = self.read_protein_sequence(file_path)
        if not sequence:
            return
        
        # Analyze protein
        try:
            self.status_var.set("Analyzing protein properties...")
            protein_analysis = ProteinAnalysis(sequence)
            
            # Calculate properties
            mw = protein_analysis.molecular_weight()
            pI = protein_analysis.isoelectric_point()
            aromaticity = protein_analysis.aromaticity()
            instability_index = protein_analysis.instability_index()
            gravy = protein_analysis.gravy()
            sec_struct = protein_analysis.secondary_structure_fraction()
            charge_at_ph = protein_analysis.charge_at_pH(ph)
            charge_at_pI = protein_analysis.charge_at_pH(pI)
            
            # Store results
            self.results = {
                'file_path': file_path,
                'sequence': sequence,
                'ph': ph,
                'mw': mw,
                'pI': pI,
                'aromaticity': aromaticity,
                'instability_index': instability_index,
                'gravy': gravy,
                'sec_struct': sec_struct,
                'charge_at_ph': charge_at_ph,
                'charge_at_pI': charge_at_pI
            }
            
            # Display results
            self.display_results()
            self.status_var.set("Analysis complete")
            
        except Exception as e:
            messagebox.showerror("Analysis Error", f"Error during analysis: {str(e)}")
            self.status_var.set("Error occurred")
    
    def display_results(self):
        results = self.results
        seq_preview = results['sequence'][:50] + "..." if len(results['sequence']) > 50 else results['sequence']
        
        # Create the report
        report = f"{' PROTEIN ANALYSIS REPORT ':=^80}\n\n"
        report += f"Source File: {os.path.basename(results['file_path'])}\n"
        report += f"Sequence Length: {len(results['sequence'])} amino acids\n"
        report += f"Sequence Preview: {seq_preview}\n\n"
        
        report += f"{' Protein Properties ':-^80}\n"
        report += f"Molecular Weight: {results['mw']:.2f} Da\n"
        report += f"Isoelectric Point (pI): {results['pI']:.4f}\n"
        report += f"Aromaticity: {results['aromaticity']:.4f}\n"
        report += f"Instability Index: {results['instability_index']:.2f}\n"
        report += f"GRAVY (Hydrophobicity): {results['gravy']:.4f}\n\n"
        
        report += f"Secondary Structure Fractions:\n"
        report += f"  Helix: {results['sec_struct'][0]:.4f}\n"
        report += f"  Turn: {results['sec_struct'][1]:.4f}\n"
        report += f"  Sheet: {results['sec_struct'][2]:.4f}\n\n"
        
        report += f"{' Charge Analysis ':-^80}\n"
        report += f"Charge at pH {results['ph']}: {results['charge_at_ph']:.4f}\n"
        report += f"Charge at pI (pH {results['pI']:.4f}): {results['charge_at_pI']:.4f}\n\n"
        
        # Interpretation
        charge_diff = results['charge_at_ph'] - results['charge_at_pI']
        if results['charge_at_ph'] > 0:
            charge_status = "positively charged"
        elif results['charge_at_ph'] < 0:
            charge_status = "negatively charged"
        else:
            charge_status = "neutral"
            
        if abs(charge_diff) < 0.1:
            comparison = "very close to"
        elif charge_diff > 0:
            comparison = "higher than"
        else:
            comparison = "lower than"
        
        report += f"Interpretation:\n"
        report += f"At pH {results['ph']}, the protein is {charge_status}.\n"
        report += f"The charge at pH {results['ph']} is {comparison} the charge at its pI.\n"
        report += f"The protein will {'not ' if abs(charge_diff) < 0.5 else ''}behave differently in electrophoresis at this pH compared to its pI.\n\n"
        
        report += f"{'=' * 80}"
        
        # Display in text widget
        self.output_text.delete(1.0, tk.END)
        self.output_text.insert(tk.END, report)
    
    def clear_results(self):
        self.output_text.delete(1.0, tk.END)
        self.file_path.set("")
        self.ph_value.set("7.0")
        self.results = {}
        self.status_var.set("Ready")
    
    def save_report(self):
        if not self.results:
            messagebox.showwarning("Save Error", "No results to save")
            return
        
        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            title="Save Analysis Report"
        )
        
        if not file_path:
            return
        
        try:
            # Get the current content of the output text widget
            report_content = self.output_text.get(1.0, tk.END)
            
            with open(file_path, 'w') as f:
                f.write(report_content)
            
            self.status_var.set(f"Report saved to {os.path.basename(file_path)}")
            messagebox.showinfo("Success", "Report saved successfully!")
        except Exception as e:
            messagebox.showerror("Save Error", f"Failed to save report: {str(e)}")

if __name__ == "__main__":
    root = tk.Tk()
    app = ProteinAnalyzerApp(root)
    root.mainloop()