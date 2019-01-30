import numpy as np
from tqdm import tqdm
from .util import number_of_lines_in_file, read_sam
import logging
from offsetbasedgraph import Graph, SequenceGraph, NumpyIndexedInterval
import mappy as mp
from .single_read_aligner import SingleSequenceAligner


class SamToGraphAligner:
    def __init__(self, graph_dir, chromosome, sam_file_name, minimum_mapq_to_graphalign=60,
                 minimum_score_to_graphalign=0):
        self.graph_dir = graph_dir
        self.chromosome = chromosome
        self.sam_file_name = sam_file_name

        self.n_skipped_supplementary = 0
        self.n_skipped_low_mapq = 0
        self.n_skipped_low_score = 0
        self.n_aligned = 0
        self.n_did_not_align = 0
        self.minimum_mapq_to_graphalign=minimum_mapq_to_graphalign
        self.minimum_score_to_graphalign=minimum_score_to_graphalign
        self.reads_skipped_because_low_mapq_and_multimapping = set()
        self.out_file = open(self.sam_file_name + ".graphalignments", "w")
        self._read_graph_data()


        #with open(self.sam_file_name + ".multimapping.txt", "w") as f:
        #    f.writelines((name + "\n" for name in self.reads_skipped_because_low_mapq_and_multimapping))

        #logging.info("%d reads have low mapq and are multimapping. Wrote these to file %s." % (
        #             len(self.reads_skipped_because_low_mapq_and_multimapping), self.sam_file_name + ".multimapping.txt"))

    def _read_graph_data(self):
        chromosome = self.chromosome
        graph = Graph.from_file(self.graph_dir + chromosome + ".nobg")
        self.graph = graph
        self.sequence_graph = SequenceGraph.from_file(self.graph_dir + chromosome + ".nobg.sequences")
        linear_path = NumpyIndexedInterval.from_file(self.graph_dir + chromosome + "_linear_pathv2.interval")
        self.linear_path = linear_path

    def _align_sam_record(self, record):
        if record.mapq < self.minimum_mapq_to_graphalign:
            self.n_skipped_low_mapq += 1
            return

        if record.score < self.minimum_score_to_graphalign:
            self.n_skipped_low_score += 1
            return

        sequence = record.sequence
        linear_pos = int(record.start) - 0  # Appearently 0-indexed when read by pysam
        position = self.linear_path.position_at_offset(linear_pos)

        if record.alternative_alignments is not None:
            if record.mapq < 60 and len(record.alternative_alignments.split(";")) > 10:
                self.reads_skipped_because_low_mapq_and_multimapping.add(record.name)
                return

        if record.is_reverse:
            # Is mapped to reverse strand
            sequence = mp.revcomp(sequence)
            #logging.warning("Found reverse alignment. Not implemented now, ignoring")

        sequence = self.sequence_graph._letter_sequence_to_numeric(np.array(list(sequence.lower())))
        # logging.info("Seq: %s" % sequence)
        aligner = SingleSequenceAligner(self.graph, self.sequence_graph, position.region_path_id,
                                        int(position.offset), sequence, n_mismatches_allowed=7, print_debug=False)
        aligner.align()
        a = aligner.get_alignment()
        if a:
            self.n_aligned += 1
            self.out_file.writelines(["%s\t%s\t%d\n" % (record.name, a.to_file_line(), aligner.n_mismatches_so_far)])
        else:
            self.n_did_not_align += 1

    def align_sam(self):
        if self.chromosome == "X":
            progress_position = 23
        else:
            progress_position = int(self.chromosome)

        #logging.info("Aligning %s" % self.sam_file_name)
        for record in tqdm(read_sam(self.sam_file_name, check_sq=False), desc="Chromosome " + str(self.chromosome),
                            total=number_of_lines_in_file(self.sam_file_name), position=progress_position):
            self._align_sam_record(record)

        logging.info("%d reads not aligned" % self.n_did_not_align)
        logging.info("%d skipped because low mapq" % self.n_skipped_low_mapq)
        logging.info("%d skipped because low score" % self.n_skipped_low_score)



