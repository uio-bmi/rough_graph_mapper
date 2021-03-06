import subprocess
import logging
from tqdm import tqdm
import pysam
from collections import defaultdict
import numpy as np
import pickle
import sys


def read_fasta(file_name):
    i = 0
    f = open(file_name)
    for line in f:
        if i % 5000000 == 0:
            logging.info("%d fasta sequence processed" % i)

        i += 1
        if line.startswith(">"):
            record_name = line.strip().replace(">", "")
        else:
            sequence = line.strip()
            yield record_name, sequence


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
        name = base_name + "_chr" + chrom + ".sam"
        outfiles[chrom] = open(name, "w")
        logging.info("Will write to %s" % name)

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
                    alternative_alignments=None, pysam_object=None, text_line=None):
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
        self.text_line = text_line
        self.pysam_object = pysam_object

    def set_start(self, new_start):
        self.start = int(new_start)
        if self.text_line != None:
            self.text_line[3] = str(self.start)
        else:
            self.pysam_object.reference_start = new_start

    def set_end(self, new_end):
        self.end = new_end
        self.pysam_object.reference_end = new_end

    def set_mapq(self, new_mapq):
        self.mapq = new_mapq
        self.pysam_object.mapping_quality = new_mapq


def read_sam_fast(sam_file_name):
    with open(sam_file_name) as f:
        for line in f:

            if line.startswith("@"):
                continue

            l = line.split()
            if len(l) <= 0:
                continue

            if int(l[1]) >= 2048:
                # supplementary (note: secondary alignments are still accepted)
                continue

            try:
                ref = l[2]
                position = int(l[3])
                name = l[0]
                mapq = int(l[4])
                score = int(l[13].replace("AS:i:", ""))
            except IndexError:
                logging.warning("Could not parse line %s" % line)
                continue

            yield Alignment(name, ref, position, None, None, None, None, mapq, score, text_line=l)


def read_sam(sam_file_name, chr=None, start=None, stop=None, skip_supplementary=True, check_sq=True,
             return_pysam_objects=False, skip_secondary=False):
    """ Wrapper around pysam"""
    f = pysam.AlignmentFile(sam_file_name, "r", check_sq=check_sq)
    for a in f.fetch(chr, start, stop, until_eof=True):
        if skip_supplementary and a.is_supplementary:
            continue

        if skip_secondary and a.is_secondary:
            continue

        if return_pysam_objects:
            yield a
            continue

        alternative_alignments = None
        if a.has_tag("XA"):
            alternative_alignments = a.get_tag("XA")

        if a.has_tag("AS"):
            score = a.get_tag("AS")
        else:
            score = 0

        try:
            reference_name = a.reference_name
        except ValueError:
            reference_name = "0"


        alignment_object = Alignment(a.query_name, reference_name, a.reference_start, a.reference_end, a.query_sequence,
                                     a.is_reverse, a.flag, a.mapping_quality, score, alternative_alignments, a)

        yield alignment_object


def number_of_lines_in_file(file_name):
    return sum(1 for _ in open(file_name, 'rb'))


def merge_single_line_sams(sam1_file_name, sam2_file_name):
    n_changed_to_sam2 = 0
    sam1 = open(sam1_file_name)
    sam2 = open(sam2_file_name)

    for i, line1 in enumerate(sam1):
        if line1.startswith("@"):
            print(line1.strip())
            continue

        if i % 100000 == 0:
            logging.info("%d lines processed." % i)

        l = line1.split()
        id1 = l[0]

        try:
            alignment_score = int(l[13].replace("AS:i:", ""))
        except IndexError:
            # Probably not aligned, skip
            print(line1.strip())
            continue

        for line2 in sam2:
            l2 = line2.split()
            if line2.startswith("@"):
                continue

            id2 = l2[0]

            assert id2 == id1

            try:
                alignment_score2 = int(l2[13].replace("AS:i:", ""))
            except IndexError:
                alignment_score2 = 0

            if alignment_score2 > alignment_score:
                print(line2.strip())
                n_changed_to_sam2 += 1
            else:
                print(line1.strip())

            break  # Only read 1 line in sam2

    logging.info("%d lines changed to sam2")


def merge_sams2(sam1_file_name, sam2_file_name, scores_are_double=False, only_score_lowering=False):
    n_changed_to_sam2 = 0
    n_minimap_better_score = 0
    n_mapq_lowered = 0
    sam1 = open(sam1_file_name)
    sam2 = open(sam2_file_name)

    lowered = []

    if only_score_lowering:
        logging.info("Will only change when scores are higher (and then choose the minimum mapq)")

    for i, line1 in enumerate(sam1):
        if line1.startswith("@"):
            print(line1.strip())
            continue

        
        if i % 100000 == 0:
            logging.info("%d lines processed. %d changed to sam 2, %d mapqs lowered" % (i, n_changed_to_sam2, n_mapq_lowered))

        l = line1.split()
        id1 = l[0]
        if int(l[1]) >= 256:
            continue  # Supplementary


        try:
            alignment_score = int(l[13].replace("AS:i:", ""))
            mapq = int(l[4])
            position1 = int(l[3])
        except IndexError:
            # Probably not aligned, skip
            alignment_score = 0
            mapq = 0
            position1 = 0
            #print(line1.strip())
            #continue

        #logging.info("Sam 1: %s" % (id1))

        for line2 in sam2:
            l2 = line2.split()
            if line2.startswith("@"):
                continue

            if int(l2[1]) >= 256:
                # secondary or supplementary alignment
                continue

            id2 = l2[0]
            #logging.info("   Sam 1: %s" % (id2))

            if id2 != id1:
                logging.error("Id1 %s != id2 %s. Is file not sorted, or are some alignments missing? Skipping for now." % (id1, id2))
                if int(id2) < int(id1):
                    continue
                else:
                    break

            assert id2 == id1

            try:
                alignment_score2 = int(l2[13].replace("AS:i:", ""))
                if scores_are_double:
                    alignment_score2 //= 2  # Divide by two, assuming this is minimap

                #alignment_score2 = int(l2[13].replace("AS:i:", ""))
                mapq2 = int(l2[4])
                position2 = int(l2[3])
            except IndexError:
                alignment_score2 = 0
                mapq2 = 0

            if alignment_score2 >= alignment_score and abs(position1 - position2) > 10:
                if not only_score_lowering:
                    l2[4] = str(0)  # Lower mapq, we have multiple good hits
                else:
                    # New mapq should be the minimum mapq
                    l2[4] = str(min(mapq, mapq2))

                print('\t'.join(l2).strip())
                n_changed_to_sam2 += 1
            elif alignment_score2 >= alignment_score * 1.0 and mapq2 < mapq and alignment_score2 > 140:
                # Same score, but lower mapq -- we want the lower mapq
                l[4] = str(mapq2)
                n_mapq_lowered += 1
                lowered.append(id1)
                print('\t'.join(l).strip())
            else:
                print(line1.strip())

            break  # Only read 1 line in sam2

    with open("lowered.txt", "w") as f:
        for l in lowered:
            f.writelines([l + "\n"])

    logging.info("%d lines changed to sam2")
    logging.info("%d mapqs lowered")


def run_hybrid_between_bwa_and_minimap(reference_file_name, fasta_file_name,
                                       output_file_name, bwa_arguments="", minimap_arguments="-a --sr"):
    logging.info("Running hybrid between bwa and minimap. Will run both and adjust mapq's")
    run_bwa_mem(reference_file_name, fasta_file_name, output_file_name + ".bwa.tmp", arguments=bwa_arguments)
    run_minimap2(reference_file_name, fasta_file_name, output_file_name + ".minimap.tmp", arguments=minimap_arguments)
    merge_sams2(output_file_name + ".bwa.tmp", output_file_name + ".minimap.tmp", output_file_name)


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
    command = "/home/ivar/dev/minimap2-2.17_x64-linux/minimap2 " + arguments + " " + reference_file_name + " " + fasta_file_name
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
