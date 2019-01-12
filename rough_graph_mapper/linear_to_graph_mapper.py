import os
from .util import split_sam_by_chromosomes, run_bwa_mem, run_hybrid_between_bwa_and_minimap
from offsetbasedgraph import Graph, SequenceGraph, NumpyIndexedInterval
from multiprocessing import Process
import logging
from .sam_to_graph_aligner import SamToGraphAligner
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s: %(message)s')


def map_single_chromosome(file_base_name, graph_dir, chromosome, minimum_mapq_to_graphalign=60):
    sam_file_name = file_base_name + "_chr" + chromosome + ".sam"
    aligner = SamToGraphAligner(graph_dir, chromosome, sam_file_name, minimum_mapq_to_graphalign)
    aligner.align_sam()


class LinearToGraphMapper:
        def __init__(self, fasta_file_name, linear_reference_file_name, graph_dir, chromosomes, minimum_mapq_to_graphalign=60, write_final_alignments_to_file=None):
            self.chromosomes = chromosomes
            self.graph_dir = graph_dir
            self.minimum_mapq_to_graphalign = minimum_mapq_to_graphalign
            self.write_final_alignments_to_file=write_final_alignments_to_file

            self.base_name = '.'.join(fasta_file_name.split(".")[:-1])
            # First align to linear reference
            run_hybrid_between_bwa_and_minimap(linear_reference_file_name, fasta_file_name, self.base_name + ".sam",
                                               bwa_arguments="-t 10",
                                               minimap_arguments="-t 40 -k19 -w11 --sr --frag=yes -A2 -B8 -O12,32 -E2,1 -r50 -p.5 -N20 -f90000,180000 -n2 -m20 -s40 -g200 -2K50m --heap-sort=yes -N 7 -a")
            #run_bwa_mem(linear_reference_file_name, fasta_file_name, self.base_name + ".sam", arguments="-t 10")
            assert os.path.isfile(self.base_name + ".sam"), "No sam file generated. Did BWA MEM fail?"

            # Split sam by chromosome
            split_sam_by_chromosomes(self.base_name + ".sam", chromosomes)
            self.map_all_chromosomes()

        def map_all_chromosomes(self):
            processes = []
            for chromosome in self.chromosomes:
                process = Process(target=map_single_chromosome, args=(self.base_name, self.graph_dir, chromosome, self.minimum_mapq_to_graphalign))
                process.start()

                processes.append(process)

            for p in processes:
                p.join()

            # Merge results from all chromosomes
            out_file = None
            out_file_name = None
            if self.write_final_alignments_to_file is not None:
                out_file_name = self.write_final_alignments_to_file
                out_file = open(out_file_name, "w")

            for chromosome in self.chromosomes:
                with open(self.base_name + "_chr" + chromosome + ".sam", "r") as f:
                    for line in f:
                        if self.write_final_alignments_to_file is not None:
                            out_file.writelines([line])
                        else:
                            print(line.strip())

            if self.write_final_alignments_to_file is not None:
                logging.info("Merged all graphalignments into file %s" % out_file_name)
                out_file.close()

            logging.info("Done mapping reads")

if __name__ == "__main__":
    mapper = LinearToGraphMapper("sim_2m.fa", "../data/human_full/testreference.fa", "../data/human_full/",
                                 [str(chrom) for chrom in range(1, 23)] + ["X"])
