import tkinter as tk
from tkinter import filedialog, messagebox
from Bio import SeqIO
from PIL import Image, ImageTk
from NCBIBlastParser import NCBIBlastParser
from io import StringIO
import os


class BLASTApp:
    def __init__(self, master):
        # Initialize the application
        self.master = master
        master.title("BLASTApp")  # Set the title of the window

        # Create text entry field for input FASTA sequences
        self.fasta_entry = tk.Text()

        # Set the background color and margin for UI elements
        self.background_color = "light grey"
        self.margin_x = 35

        # Initialize NCBI Blast Parser
        self.ncbi_parser = NCBIBlastParser()

        # Initialize file_path attribute
        self.file_path = None

        # Header
        self.header_frame = tk.Frame(master, bg="#007786", padx=15)
        self.header_frame.pack(fill=tk.X)

        # App name label
        self.app_name_label = tk.Label(self.header_frame, text="NCBI BLAST Parser", font=("livvic", 18, "bold"),
                                       bg="#007786", fg="light grey")
        self.app_name_label.pack(side=tk.LEFT)

        # Logo
        self.logo = Image.open("logo.png")
        self.logo = self.logo.resize((120, 50))
        self.logo = ImageTk.PhotoImage(self.logo)
        self.logo_label = tk.Label(self.header_frame, image=self.logo, bg="#007786")
        self.logo_label.pack(side=tk.RIGHT)

        # Spacer
        space_line = tk.Frame(master, height=20, bg=self.background_color)
        space_line.pack(fill=tk.X)

        # Main content frame
        self.content_frame = tk.Frame(master, padx=self.margin_x, bg=self.background_color)
        self.content_frame.pack(fill=tk.BOTH, expand=True)

        # Left Side
        self.left_frame = tk.Frame(self.content_frame, bg=self.background_color)
        self.left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Input FASTA sequence widget
        self.fasta_label = tk.Label(self.left_frame, text="Enter or Upload FASTA Sequences:", bg=self.background_color)
        self.fasta_label.pack(anchor=tk.W)

        # Entry field for FASTA sequences
        self.fasta_entry_frame = tk.Frame(self.left_frame, bg=self.background_color)
        self.fasta_entry_frame.pack(anchor=tk.W)

        self.fasta_entry = tk.Text(self.fasta_entry_frame, height=8, width=65)
        self.fasta_entry.pack(side=tk.LEFT)

        # Scrollbar for FASTA entry
        self.fasta_scrollbar = tk.Scrollbar(self.fasta_entry_frame, orient=tk.VERTICAL, command=self.fasta_entry.yview)
        self.fasta_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.fasta_entry.config(yscrollcommand=self.fasta_scrollbar.set)

        # Buttons for FASTA operations
        self.button_row_frame = tk.Frame(self.left_frame, bg=self.background_color)
        self.button_row_frame.pack(fill=tk.X)

        # Browse FASTA file button
        self.browse_fasta_button = tk.Button(self.button_row_frame, text="Choose File", command=self.browse_fasta_file)
        self.browse_fasta_button.pack(side=tk.LEFT)

        # Clear input button
        self.clear_inputs_button = tk.Button(self.button_row_frame, text="Clear",
                                             command=self.clear_fasta_inputs)
        self.clear_inputs_button.pack(side=tk.LEFT)

        # Run BLAST search button
        self.run_blast_button = tk.Button(self.button_row_frame, text="Run BLAST Search", command=self.run_blast_search)
        self.run_blast_button.pack(side=tk.LEFT)

        # Count sequence types button
        self.count_button = tk.Button(self.button_row_frame, text="Count Sequence Types",
                                      command=self.count_sequence_types)
        self.count_button.pack(side=tk.LEFT)

        # Right Side
        self.right_frame = tk.Frame(self.content_frame, bg=self.background_color)
        self.right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # Output panel label
        self.output_label = tk.Label(self.right_frame, text="Output Panel:", bg=self.background_color)
        self.output_label.pack(anchor=tk.W)

        # Output text widget
        self.output_text_frame = tk.Frame(self.right_frame, bg=self.background_color)
        self.output_text_frame.pack(anchor=tk.W)

        self.output_text = tk.Text(self.output_text_frame, height=30, width=68)
        self.output_text.pack(side=tk.LEFT)

        # Search box
        self.search_frame = tk.Frame(self.right_frame, bg=self.background_color)
        self.search_frame.pack(anchor=tk.W)

        # Spacer
        self.blank_line = tk.Label(self.search_frame, text="", height=1, bg=self.background_color)
        self.blank_line.pack(side=tk.TOP)

        # Label
        search_label = tk.Label(self.search_frame, text="Search in the Output Results:", bg=self.background_color)
        search_label.pack(side=tk.TOP, anchor="w")

        self.search_entry = tk.Entry(self.search_frame, width=50)
        self.search_entry.pack(side=tk.LEFT)

        self.search_button = tk.Button(self.search_frame, text="Search", command=self.search_output)
        self.search_button.pack(side=tk.LEFT)

        # Scrollbar for output text
        self.output_scrollbar = tk.Scrollbar(self.output_text_frame, orient=tk.VERTICAL, command=self.output_text.yview)
        self.output_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.output_text.config(yscrollcommand=self.output_scrollbar.set)

        # Search for motif section
        self.motif_frame = tk.Frame(self.left_frame, bg=self.background_color)
        self.motif_frame.pack(fill=tk.X)

        # Spacer
        self.blank_line = tk.Label(self.motif_frame, text="", height=1, bg=self.background_color)
        self.blank_line.pack(side=tk.TOP)

        # Browse motif file label
        self.browse_motif_file_label = tk.Label(self.motif_frame, text="Select a Valid XML File:",
                                                bg=self.background_color)
        self.browse_motif_file_label.pack(side=tk.LEFT)

        # Entry field for motif file path
        self.file_entry = tk.Entry(self.motif_frame, width=51)
        self.file_entry.pack(side=tk.LEFT)

        # Browse motif file button
        self.browse_motif_file_button = tk.Button(self.motif_frame, text="Choose File", command=self.browse_motif_file)
        self.browse_motif_file_button.pack(side=tk.LEFT)

        # Motif input section
        self.motif_input_frame = tk.Frame(self.left_frame, bg=self.background_color)
        self.motif_input_frame.pack(anchor=tk.W)

        # Spacer
        self.blank_line = tk.Label(self.motif_input_frame, text="", height=1, bg=self.background_color)
        self.blank_line.pack(side=tk.TOP)

        # Search motif label
        self.search_motif_label_frame = tk.Frame(self.motif_input_frame, bg=self.background_color)
        self.search_motif_label_frame.pack(side=tk.TOP, anchor=tk.W)
        self.motif_file_label = tk.Label(self.search_motif_label_frame, text="Search for a Sequence Motif:",
                                         bg=self.background_color)
        self.motif_file_label.pack(side=tk.LEFT)

        # Entry field for motif
        self.motif_label = tk.Label(self.motif_input_frame, text="Enter a Sequence Motif:", bg=self.background_color)
        self.motif_label.pack(side=tk.LEFT)
        self.motif_entry = tk.Entry(self.motif_input_frame, width=50)
        self.motif_entry.pack(side=tk.LEFT)

        # Search motif button
        self.search_motif_button = tk.Button(self.left_frame, text="Search", command=self.search_for_motif)
        self.search_motif_button.pack(anchor=tk.W)

        # Parse BLAST hits section
        self.parse_frame = tk.Frame(self.left_frame, bg=self.background_color)
        self.parse_frame.pack(fill=tk.X)

        # Spacer
        self.blank_line = tk.Label(self.parse_frame, text="", height=1, bg=self.background_color)
        self.blank_line.pack(side=tk.TOP)

        # Parse BLAST hits label
        self.parse_file_label = tk.Label(self.parse_frame, text="Parse BLAST Hits:",
                                         bg=self.background_color)
        self.parse_file_label.pack(side=tk.LEFT)

        # Parse BLAST hits button
        self.parse_button = tk.Button(self.left_frame, text="Parse", command=self.parse_blast_hits)
        self.parse_button.pack(anchor=tk.W)

        # Plot E-value distribution
        self.plot_frame = tk.Frame(self.left_frame, bg=self.background_color)
        self.plot_frame.pack(fill=tk.X)

        # Spacer
        self.blank_line = tk.Label(self.plot_frame, text="", height=1, bg=self.background_color)
        self.blank_line.pack(side=tk.TOP)

        # Plot file label
        self.plot_file_frame = tk.Frame(self.plot_frame, bg=self.background_color)
        self.plot_file_frame.pack(fill=tk.X)

        self.plot_file_label = tk.Label(self.plot_file_frame, text="Plot the E-value Distribution:",
                                        bg=self.background_color)
        self.plot_file_label.pack(side=tk.LEFT)

        # Output file entry
        self.output_frame = tk.Frame(self.plot_frame, bg=self.background_color)
        self.output_frame.pack(fill=tk.X)

        self.output_label = tk.Label(self.output_frame, text="Insert Output File Name:", bg=self.background_color)
        self.output_label.pack(side=tk.LEFT)
        self.output_entry = tk.Entry(self.output_frame, width=49)
        self.output_entry.pack(side=tk.LEFT)

        # Plot button
        self.plot_button = tk.Button(self.left_frame, text="Plot",
                                     command=self.plot_evalue_distribution)
        self.plot_button.pack(anchor=tk.W)

        # Filter BLAST hits by E-value
        self.filter_frame = tk.Frame(self.left_frame, bg=self.background_color)
        self.filter_frame.pack(fill=tk.X)

        self.filter_file_frame = tk.Frame(self.filter_frame, bg=self.background_color)
        self.filter_file_frame.pack(fill=tk.X)

        # Spacer
        self.blank_line = tk.Label(self.filter_file_frame, text="", height=1, bg=self.background_color)
        self.blank_line.pack(side=tk.TOP)

        # Filter file label
        self.filter_file_label = tk.Label(self.filter_file_frame, text="Filter BLAST Hits by E-value:",
                                          bg=self.background_color)
        self.filter_file_label.pack(side=tk.LEFT)

        # Threshold entry
        self.threshold_frame = tk.Frame(self.filter_frame, bg=self.background_color)
        self.threshold_frame.pack(fill=tk.X)

        self.filter_label = tk.Label(self.threshold_frame, text="Enter E-value Threshold:", bg=self.background_color)
        self.filter_label.pack(side=tk.LEFT)
        self.filter_threshold_entry = tk.Entry(self.threshold_frame, width=49)
        self.filter_threshold_entry.pack(side=tk.LEFT)

        # Filter button
        self.filter_button = tk.Button(self.left_frame, text="Filter",
                                       command=self.filter_blast_hits_by_evalue)
        self.filter_button.pack(anchor=tk.W)

        # Reset button
        self.clear_all_button = tk.Button(self.left_frame, text="Reset", bg="#007786", fg="white",
                                          command=self.clear_all_inputs)
        self.clear_all_button.pack(side=tk.LEFT, pady=5)

        # Footer frame
        self.footer_frame = tk.Frame(master, bg="grey", padx=10, pady=5)
        self.footer_frame.pack(fill=tk.X, side=tk.BOTTOM)

        # Footer label
        self.footer_label = tk.Label(self.footer_frame, text="© 2024 BLASTApp. All rights reserved.", bg="grey",
                                     fg="white")
        self.footer_label.pack()

    # Function to browse and load a FASTA file into the entry field
    def browse_fasta_file(self):
        # Clear previous content
        self.fasta_entry.delete(1.0, tk.END)
        # Open file dialog and store file path
        self.file_path = filedialog.askopenfilename()
        if self.file_path:
            with open(self.file_path, 'r') as file:
                fasta_sequence = file.read()
                # Insert FASTA sequence into entry field
                self.fasta_entry.insert(tk.END, fasta_sequence)

    # Function to clear file path entry field
    def clear_browse_entries(self):
        self.file_entry.delete(0, tk.END)

    # Function to browse and load a motif XML file
    def browse_motif_file(self):
        # Clear previous file path entry
        self.clear_browse_entries()
        # Clear previous output
        self.output_text.delete(1.0, tk.END)
        # Open file dialog
        file_path = filedialog.askopenfilename()
        if file_path:
            # Insert file path into entry field
            self.file_entry.insert(0, file_path)

    # Function to run BLAST search
    # Function to run BLAST search
    def run_blast_search(self):
        # Get FASTA sequence from entry field
        fasta_sequence = self.fasta_entry.get("1.0", tk.END).strip()
        if fasta_sequence:
            sequences = []
            # Parse FASTA sequence
            for record in SeqIO.parse(StringIO(fasta_sequence), "fasta"):
                sequences.append(record)
            # Set sequences for NCBI parser
            self.ncbi_parser.sequences = sequences
            try:
                # Set output directory
                output_dir = os.path.dirname(self.file_path) if self.file_path else os.getcwd()
                # Run BLAST search
                self.ncbi_parser.run_blast_search(output_dir)
                # Show success message
                messagebox.showinfo("Success", "BLAST search completed and results saved.")
            except Exception as e:
                # Show error message if an exception occurs
                messagebox.showerror("Error", f"Error: {str(e)}")
        else:
            # Show error message if no FASTA sequence is provided
            messagebox.showerror("Error", "No FASTA sequence provided.")

    # Function to count sequence types in the provided FASTA sequence
    def count_sequence_types(self):
        # Get FASTA sequence from entry field
        fasta_sequence = self.fasta_entry.get("1.0", tk.END).strip()
        # Initialize BLAST parser with the FASTA sequence
        blast_parser = NCBIBlastParser(StringIO(fasta_sequence))
        # Count sequence types
        type_counts, sequence_types = blast_parser.count_sequence_types()
        # Clear output text widget
        self.output_text.delete(1.0, tk.END)
        # Display counts of each sequence type
        for seq_type, count in type_counts.items():
            self.output_text.insert(tk.END, f"{seq_type}: {count}\n")
        self.output_text.insert(tk.END, "\n")
        # Display sequence types
        for seq_index, seq_type in sequence_types:
            self.output_text.insert(tk.END, f"{seq_index}: {seq_type}\n")

    # Function to search for a motif in BLAST hits
    def search_for_motif(self):
        # Clear output text widget
        self.output_text.delete(1.0, tk.END)
        # Get XML file path and motif from entry fields
        xml_file = self.file_entry.get().strip()
        motif = self.motif_entry.get().strip()
        # Check if XML file path and motif are provided
        if not xml_file:
            messagebox.showerror("Error", "Please select a valid XML file.")
            return
        if not motif:
            messagebox.showerror("Error", "Please enter a motif.")
            return
        # Initialize BLAST parser
        blast_parser = NCBIBlastParser()
        # Search for motif in BLAST hits
        motifs_found = blast_parser.search_for_motif(xml_file, motif, self.output_text)

        if not self.output_text.get(1.0, tk.END).strip():
            # Display message if no motif found in BLAST records
            self.output_text.insert(tk.END, "No motif found in the blast records.")
        else:
            if motifs_found:
                return

    # Function to parse BLAST hits from XML file
    def parse_blast_hits(self):
        # Get XML file path from entry field
        blast_file = self.file_entry.get().strip()
        if not blast_file:
            # Show error message if no XML file selected
            messagebox.showerror("Error", "Please select a valid XML file.")
            return
        # Initialize BLAST parser
        blast_parser = NCBIBlastParser()
        try:
            # Parse BLAST hits from XML file
            blast_parser.parse_blast_hits(blast_file)
            # Show success message
            messagebox.showinfo("Success", "BLAST hits parsed!")
        except Exception as e:
            # Show error message if an exception occurs
            messagebox.showerror("Error", f"Error: {str(e)}")

    # Function to plot E-value distribution
    def plot_evalue_distribution(self):
        # Get XML file path and output file name from entry fields
        xml_file = self.file_entry.get().strip()
        output_file_name = self.output_entry.get().strip()
        if not os.path.isfile(xml_file):
            # Show error message if XML file is invalid
            messagebox.showerror("Error", "Please select a valid XML file.")
            return
        if not output_file_name:
            # Show error message if no output file name provided
            messagebox.showerror("Error", "Please enter an output file name.")
            return
        try:
            # Plot E-value distribution
            self.ncbi_parser.plot_evalue_distribution(xml_file, output_file_name)
            # Show success message
            messagebox.showinfo("Success", "E-value distribution plotted successfully.")
        except Exception as e:
            # Show error message if an exception occurs
            messagebox.showerror("Error", f"Error: {str(e)}")

    # Function to filter BLAST hits by E-value threshold
    def filter_blast_hits_by_evalue(self):
        # Get XML file path and E-value threshold from entry fields
        xml_file = self.file_entry.get().strip()
        evalue_threshold_entry = self.filter_threshold_entry.get().strip()
        if not xml_file:
            # Show error message if no XML file selected
            messagebox.showerror("Error", "Please select a valid XML file.")
            return
        if not evalue_threshold_entry:
            # Show error message if no E-value threshold provided
            messagebox.showerror("Error", "Please enter an E-value threshold.")
            return
        try:
            # Convert E-value threshold to float
            evalue_threshold = float(evalue_threshold_entry)
            blast_parser = NCBIBlastParser()
            # Filter BLAST hits by E-value threshold
            blast_parser.filter_blast_hits_by_evalue(xml_file, evalue_threshold)
            # Show success message
            messagebox.showinfo("Success", "BLAST hits filtered!")
        except ValueError:
            # Show error message if E-value threshold is invalid
            messagebox.showerror("Error", "Please enter a valid numeric E-value threshold.")
        except Exception as e:
            # Show error message if an exception occurs
            messagebox.showerror("Error", f"Error: {str(e)}")

    # Function to search for a text in the output text widget
    def search_output(self):
        # Get search text from entry field
        search_text = self.search_entry.get().strip().lower()
        if search_text:
            # Remove previous highlights
            self.output_text.tag_remove('highlight', '1.0', tk.END)
            start_index = '1.0'
            found = False
            while True:
                # Search for text in output text widget
                start_index = self.output_text.search(search_text, start_index, stopindex=tk.END,
                                                      nocase=True)
                if not start_index:
                    break
                found = True
                end_index = f"{start_index}+{len(search_text)}c"
                # Add highlight tag to matching text
                self.output_text.tag_add('highlight', start_index, end_index)
                start_index = end_index
            if not found:
                # Display message if no matches found
                messagebox.showinfo("Search Results", "No matches found.")
            else:
                # Configure highlight tag
                self.output_text.tag_config('highlight', background='yellow')

    # Function to clear all input fields and output text widget
    def clear_all_inputs(self):
        # Clear all input fields and output text widget
        self.fasta_entry.delete(1.0, tk.END)
        self.file_entry.delete(0, tk.END)
        self.motif_entry.delete(0, tk.END)
        self.output_entry.delete(0, tk.END)
        self.filter_threshold_entry.delete(0, tk.END)
        self.output_text.delete(1.0, tk.END)
        self.search_entry.delete(0, tk.END)

    # Function to clear FASTA input text widget
    def clear_fasta_inputs(self):
        self.fasta_entry.delete(1.0, tk.END)
        self.output_entry.delete(0, tk.END)


def main():
    root = tk.Tk()
    root.geometry("1200x700")
    root.resizable(False, False)
    BLASTApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
