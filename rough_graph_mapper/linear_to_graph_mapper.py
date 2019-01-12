import os
from .util import split_sam_by_chromosomes, run_bwa_mem
from offsetbasedgraph import Graph, SequenceGraph, NumpyIndexedInterval
from multiprocessing import Process
import logging
from .sam_to_graph_aligner import SamToGraphAligner
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s: %(message)s')


def map_single_chromosome(file_base_name, graph_dir, chromosome):
    sam_file_name = file_base_name + "_chr" + chromosome + ".sam"
    aligner = SamToGraphAligner(graph_dir, chromosome, sam_file_name)
    aligner.align_sam()


class LinearToGraphMapper:
        def __init__(self, fasta_file_name, linear_reference_file_name, graph_dir, chromosomes):
            self.chromosomes = chromosomes
            self.graph_dir = graph_dir

            self.base_name = '.'.join(fasta_file_name.split(".")[:-1])
            # First align to linear reference
            run_bwa_mem(linear_reference_file_name, fasta_file_name, self.base_name + ".sam", arguments="-t 10")

            assert os.path.isfile(self.base_name + ".sam"), "No sam file generated. Did BWA MEM fail?"

            # Split sam by chromosome
            split_sam_by_chromosomes(self.base_name + ".sam", chromosomes)

            self.map_all_chromosomes()

        def map_all_chromosomes(self):
            processes = []
            for chromosome in self.chromosomes:
                process = Process(target=map_single_chromosome, args=(self.base_name, self.graph_dir, chromosome))
                process.start()

                processes.append(process)

            for p in processes:
                p.join()

            logging.info("Done mapping reads")

if __name__ == "__main__":
    mapper = LinearToGraphMapper("sim_2m.fa", "../data/human_full/testreference.fa", "../data/human_full/",
                                 [str(chrom) for chrom in range(1, 23)] + ["X"])
