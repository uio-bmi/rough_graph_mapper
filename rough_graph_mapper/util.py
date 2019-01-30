import subprocess
import logging
from tqdm import tqdm
import pysam
from collections import defaultdict
import numpy as np


def read_fasta(file_name):
    i = 0
    out = {}
    record_name = None
    f = open(file_name)
    for line in f:
        if i % 5000000 == 0:
            logging.info("%d fasta sequence processed" % i)

        i += 1
        if line.startswith(">"):
            record_name = line.strip().replace(">", "")
        else:
            sequence = line.strip()
            out[record_name] = sequence

    return out

def read_fasta_to_numeric_sequences(file_name, sequence_graph):
    i = 0
    out = {}
    record_name = None
    f = open(file_name)
    for line in f:
        if i % 5000000 == 0:
            logging.info("%d fasta sequence processed" % i)

        i += 1
        if line.startswith(">"):
            record_name = line.strip().replace(">", "")
        else:
            sequence = line.strip()
            numeric_sequence = sequence_graph._letter_sequence_to_numeric(np.array(list(sequence.lower())))
            if np.all(numeric_sequence == 0):
                logging.debug(sequence)
                logging.debug(numeric_sequence)
                raise Exception("Could not convert")
            out[record_name] = (sequence, numeric_sequence)

    return out


def split_sam_by_chromosomes(sam_file, chromosomes):
    chromosome_index = set(chromosomes)
    base_name = sam_file.split(".")[0]
    n_not_matched = 0
    outfiles = {}
    for chrom in chromosomes:
        outfiles[chrom] = open(base_name + "_chr" + chrom + ".sam", "w")

    n_lines = number_of_lines_in_file(sam_file)

    with open(sam_file) as f:
        for i, line in enumerate(tqdm(f, total=n_lines)):
            if line.startswith("@"):
                for chrom in chromosomes:
                    outfiles[chrom].writelines([line])
                continue

            chrom = line.split()[2]
            if chrom not in chromosome_index:
                n_not_matched += 1
                continue

            outfiles[chrom].writelines([line])

    for chrom in chromosomes:
        outfiles[chrom].close()

    logging.info("Done splitting by chromosome. %d lines did not match any chromosomes or were unaligned." % n_not_matched)


class Alignment:
    def __init__(self, name, chromosome, start, end, sequence, is_reverse, flag=0, mapq=0, score=0,
                    alternative_alignments=None, pysam_object=None):
        self.sequence = sequence
        self.name = name
        self.start = start
        self.end = end
        self.is_reverse = is_reverse
        self.chromosome = chromosome
        self.score = score
        self.mapq = mapq
        self.flag = flag
        self.alternative_alignments = alternative_alignments
        self.pysam_object = pysam_object

    def set_start(self, new_start):
        self.start = new_start
        self.pysam_object.reference_start = new_start

    def set_end(self, new_end):
        self.end = new_end
        self.pysam_object.reference_end = new_end

    def set_mapq(self, new_mapq):
        self.mapq = new_mapq
        self.pysam_object.mapping_quality = new_mapq


def read_sam(sam_file_name, chr=None, start=None, stop=None, skip_supplementary=True, check_sq=True):
    """ Wrapper around pysam"""
    f = pysam.AlignmentFile(sam_file_name, "r", check_sq=check_sq)
    for a in f.fetch(chr, start, stop, until_eof=True):
        if skip_supplementary and a.is_supplementary:
            continue

        alternative_alignments = None
        if a.has_tag("XA"):
            alternative_alignments = a.get_tag("XA")

        if a.has_tag("AS"):
            score = a.get_tag("AS")
        else:
            score = 0

        alignment_object = Alignment(a.query_name, a.reference_name, a.reference_start, a.reference_end, a.query_sequence,
                                     a.is_reverse, a.flag, a.mapping_quality, score, alternative_alignments, a)

        yield alignment_object


def number_of_lines_in_file(file_name):
    return sum(1 for _ in open(file_name, 'rb'))


def select_lowest_mapq_from_two_sam_files(sam1_file_name, sam2_file_name, output_file_name):
    sam2_n_good_alignments = defaultdict(int)
    mapq_sam2 = defaultdict(int)

    out_sam = pysam.AlignmentFile(output_file_name, "w", template=pysam.AlignmentFile(sam1_file_name, "r"))

    logging.info("Correcting mapq scores. First reading sam2")
    for alignment in tqdm(read_sam(sam2_file_name), total=number_of_lines_in_file(sam2_file_name)):
        sam2_n_good_alignments[alignment.name] += 1
        mapq_sam2[alignment.name] = max(alignment.mapq, mapq_sam2[alignment.name])  # Keep best mapq

    n_adjusted = 0

    logging.info("Reading sam file 1, checking if mapqs should be lowered")
    for alignment in tqdm(read_sam(sam1_file_name), total=number_of_lines_in_file(sam1_file_name)):
        if alignment.name not in mapq_sam2:
            logging.warning("Sam file 2 does not have %s which is in sam file 2." % alignment.name)

        mapq2 = mapq_sam2[alignment.name]

        # If sam2 has low mapq due to multiple good alignments, use sam2's mapq score (we trust that one)
        if alignment.mapq > mapq2 and sam2_n_good_alignments[alignment.name] > 1:
            alignment.set_mapq(mapq2)
            n_adjusted += 1

        out_sam.write(alignment.pysam_object)

    logging.info("%d mapqs were adjusted down" % n_adjusted)
    logging.info("Correct sam file written to %s" % output_file_name)


def run_hybrid_between_bwa_and_minimap(reference_file_name, fasta_file_name,
                                       output_file_name, bwa_arguments="", minimap_arguments="-a --sr"):
    logging.info("Running hybrid between bwa and minimap. Will run both and adjust mapq's")
    run_bwa_mem(reference_file_name, fasta_file_name, output_file_name + ".bwa.tmp", arguments=bwa_arguments)
    run_minimap2(reference_file_name, fasta_file_name, output_file_name + ".minimap.tmp", arguments=minimap_arguments)
    select_lowest_mapq_from_two_sam_files(output_file_name + ".bwa.tmp", output_file_name + ".minimap.tmp",
                                          output_file_name)


def run_bwa_mem(reference_file_name, fasta_file_name, output_file_name, arguments=""):
    command = "bwa mem " + arguments + " " + reference_file_name + " " + fasta_file_name
    logging.info("Running bwa mem with command %s" % command)
    with open(output_file_name, "w") as outfile:
        process = subprocess.Popen(command.split(), stdout=outfile, stderr=subprocess.PIPE)

    output = process.stderr.read().decode("utf-8")
    logging.info(output)
    if "fail" in output or "error" in output:
        logging.error("BWA MEM seems to have failed. Are file paths correct? Check BWA output.")

    logging.info("Done running BWA MEM. Wrote alignments to %s" % output_file_name)

def run_minimap2(reference_file_name, fasta_file_name, output_file_name, arguments="-ax sr"):
    command = "minimap2 " + arguments + " " + reference_file_name + " " + fasta_file_name
    logging.info("Running Minimap 2 with command %s" % command)
    with open(output_file_name, "w") as outfile:
        process = subprocess.Popen(command.split(), stdout=outfile, stderr=subprocess.PIPE)

    output = process.stderr.read().decode("utf-8")
    logging.info(output)
    if "fail" in output or "Error" in output:
        logging.error("Minimap 2 seems to have failed. Are file paths correct? Check Minimap 2 output.")

    logging.info("Done running Minimap 2. Wrote alignments to %s" % output_file_name)

def read_fasta_file(fasta_file_name):
    pass
