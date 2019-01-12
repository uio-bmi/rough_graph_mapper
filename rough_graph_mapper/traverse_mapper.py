from util import run_hybrid_between_bwa_and_minimap
import logging


class TraverseMapper:
        def __init__(self, fasta_file_name, linear_reference_file_name, graph_dir, chromosomes):
            self.base_name = fasta_file_name.split(".")[0]

            # First map to linear references
            run_hybrid_between_bwa_and_minimap(linear_reference_file_name, fasta_file_name, self.base_name + ".sam")

            # Split by chromosome

            # Grahalign the reads that were mapped to linear references

            # Get the reads that did not map well to linear reference

            # Map these with traverser
