from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import matplotlib.pyplot as plt
import os


class NCBIBlastParser:
    def __init__(self, fasta_file=None):
        self.fasta_file = fasta_file
        if fasta_file:
            self.sequences = list(SeqIO.parse(fasta_file, 'fasta'))  # Read FASTA file if provided
        else:
            self.sequences = []  # Initialize an empty list if no FASTA file provided

    @staticmethod
    # Function to determine the type of sequence (DNA, RNA, or Amino Acid)
    def get_sequence_type(sequence):
        # Check for amino acid sequences
        if 'M' in sequence.seq:
            return 'Amino Acid'
        # Check for RNA sequences
        elif 'U' in sequence.seq:
            return 'RNA'
        else:  # Default to DNA
            return 'DNA'

    # Function to count the occurrences of each sequence type
    def count_sequence_types(self):
        # Initialize counts for each sequence type
        type_counts = {'DNA': 0, 'RNA': 0, 'Amino Acid': 0}
        sequence_types = []
        # Iterate over sequences and count their types
        for i, sequence in enumerate(self.sequences, start=1):
            # Determine sequence type
            sequence_type = self.get_sequence_type(sequence)
            # Increment count for the respective sequence type
            type_counts[sequence_type] += 1
            # Store sequence description and type
            description = sequence.description if hasattr(sequence, "description") else f">Sequence {i}"
            if not description.startswith(">"):
                description = f">{description}"
            sequence_types.append((description, sequence_type))
        return type_counts, sequence_types

    # Function to perform BLAST search for each sequence
    def run_blast_search(self, output_dir):
        # Iterate over sequences and perform BLAST search
        for sequence in self.sequences:
            # Determine sequence type
            sequence_type = self.get_sequence_type(sequence)
            # Run appropriate BLAST algorithm based on sequence type
            blast_result = self.run_appropriate_blast(sequence, sequence_type)
            # Save BLAST result
            self.save_blast_result(blast_result, sequence.id, output_dir)

    @staticmethod
    # Function to choose the appropriate BLAST algorithm based on sequence type
    def run_appropriate_blast(sequence, sequence_type):
        # Run BLASTp for amino acid sequences
        if sequence_type == 'Amino Acid':
            result_handle = NCBIWWW.qblast("blastp", "nr", sequence.seq)
        # Run BLASTn for RNA sequences
        elif sequence_type == 'RNA':
            result_handle = NCBIWWW.qblast("blastn", "nr", sequence.seq.transcribe())
        else:  # Run BLASTn for DNA sequences
            result_handle = NCBIWWW.qblast("blastn", "nr", sequence.seq)
        # Read and store BLAST result
        blast_result = result_handle.read()
        return blast_result

    @staticmethod
    # Function to save BLAST result to a file
    def save_blast_result(blast_result, sequence_id, output_dir):
        # Write BLAST result to a file
        with open(os.path.join(output_dir, f"{sequence_id}_blast_result.xml"), "w") as result_file:
            result_file.write(blast_result)

    @staticmethod
    # Function to search for a motif within BLAST results
    def search_for_motif(xml_file, motif, output_text):
        # Open and read XML file containing BLAST results
        result_handle = open(xml_file)
        blast_records = NCBIXML.read(result_handle)
        motifs_found = False
        # Iterate over BLAST records and alignments
        for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                # Check if motif exists in subject sequence
                if motif in hsp.sbjct:
                    motif_start = hsp.sbjct.find(motif) + 1  # Find start position of motif
                    motif_end = motif_start + len(motif) - 1  # Find end position of motif
                    # Output motif information
                    output_text.insert("end", f"Motif '{motif}' found in {alignment.title}, "
                                              f"Start: {motif_start}, End: {motif_end}\n")
                    output_text.insert("end", "-" * 68 + "\n")
                    motifs_found = True
        return motifs_found

    @staticmethod
    # Function to parse BLAST hits and save to a text file
    def parse_blast_hits(xml_file):
        # Open and read XML file containing BLAST results
        result_handle = open(xml_file)
        blast_records = NCBIXML.read(result_handle)
        output_directory = os.path.dirname(xml_file)  # Get directory path
        # Define output file path
        output_file_path = os.path.join(output_directory,
                                        f"{os.path.splitext(os.path.basename(xml_file))[0]}_parsed_output.txt")
        # Write parsed BLAST hits to the output file
        with open(output_file_path, "w") as output_file:
            for alignment in blast_records.alignments:
                output_file.write(f"Sequence: {alignment.title}\n")
                output_file.write(f"Length: {alignment.length}\n")
                output_file.write(f"E-value: {alignment.hsps[0].expect}\n")
                output_file.write(f"Alignment: {alignment.hsps[0].sbjct}\n\n")

    @staticmethod
    # Function to plot E-value distribution and save the plot
    def plot_evalue_distribution(xml_file, output_file_name):
        # Open and read XML file containing BLAST results
        result_handle = open(xml_file)
        blast_records = NCBIXML.read(result_handle)
        # Extract E-values from BLAST records
        evalues = [hit.expect for alignment in blast_records.alignments for hit in alignment.hsps]
        # Plot histogram of E-values
        plt.hist(evalues, bins=20, log=True)
        plt.xlabel('E-value')
        plt.ylabel('Frequency')
        plt.title('E-value Distribution')
        output_directory = os.path.dirname(xml_file)  # Get directory path
        # Define output file path
        output_file_path = os.path.join(output_directory, output_file_name)
        # Save plot to specified file path
        plt.savefig(output_file_path)
        plt.close()

    @staticmethod
    # Function to filter BLAST hits by E-value and save filtered hits to a text file
    def filter_blast_hits_by_evalue(xml_file, threshold):
        # Open and read XML file containing BLAST results
        result_handle = open(xml_file)
        blast_records = NCBIXML.read(result_handle)
        output_directory = os.path.dirname(xml_file)  # Get directory path
        # Define output file path
        output_file = os.path.join(output_directory,
                                   f"{os.path.splitext(os.path.basename(xml_file))[0]}_filtered_hits.txt")
        # Write filtered BLAST hits to the output file
        with open(output_file, "w") as filtered_file:
            for alignment in blast_records.alignments:
                for hit in alignment.hsps:
                    # Check if E-value is below threshold
                    if hit.expect < threshold:
                        filtered_file.write(f"Sequence: {alignment.title}\n")
                        filtered_file.write(f"Length: {alignment.length}\n")
                        filtered_file.write(f"E-value: {hit.expect}\n")
                        filtered_file.write(f"Alignment: {hit.sbjct}\n\n")
