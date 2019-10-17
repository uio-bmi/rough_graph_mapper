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
                mapq2 = 0

            if alignment_score2 > alignment_score:
                print(line2.strip())
                n_changed_to_sam2 += 1
            else:
                print(line1.strip())

            break  # Only read 1 line in sam2

    logging.info("%d lines changed to sam2")


def improve_mapping_with_two_sams(sam1_file_name, sam2_file_name):

    n_changed_to_sam2 = 0
    n_mapq_lowered = 0
    n_untouched = 0

    sam1 = open(sam1_file_name)
    sam2 = open(sam2_file_name)

    prev_line1_id = -1
    sam2_best_alignment_score = 0
    sam2_best_alignment = None
    sam2_n_good_alignments = 0
    sam2_best_mapq = 0
    i2 = 0

    lowered_ids = set()
    last_line = None
    for i, line1 in enumerate(sam1):
        last_line = line1
        if line1.startswith("@"):
            print(line1.strip())
            continue

        if i % 100000 == 0:
            logging.info("%d lines processed. On line %d in sam2. %d changed to sam 2. %d lowered mapq. %d untouched" %
                         (i, i2, n_changed_to_sam2, n_mapq_lowered, n_untouched))

        l = line1.split()
        id1 = l[0]

        if int(l[1]) >= 256:
            # secondary or supplementary alignment
            continue

        mapq = int(l[4])
        try:
            alignment_score = int(l[13].replace("AS:i:", ""))
        except IndexError:
            # Probably not aligned, skip
            print(line1.strip())
            continue

        assert id1 != prev_line1_id, "Sam 1 contains duplicate entries with same ID"
        assert int(id1) > int(prev_line1_id), "prev id: %d, id now: %d. Line %s" % (int(prev_line1_id), int(id1), line1)
        prev_line1_id = id1

        #logging.info("%s at %s:%s with mapq %d and score %d" % (id1, l[2], l[3], mapq, alignment_score))

        # read sam2 as long as we have this id
        for line2 in sam2:
            if line2.startswith("@"):
                continue

            i2 += 1
            l2 = line2.split()
            id2 = l2[0]

            if int(l2[1]) >= 2048:
                # supplementary alignment
                continue

            if int(id2) < int(id1):
                logging.error("File 2 has id %s > file 1 id %s." % (id2, id1))
                logging.error("Line 1: %s" % line1)
                logging.error("Line 2: %s" % line2)
                continue

            try:
                alignment_score2 = int(l2[13].replace("AS:i:", ""))
                mapq2 = int(l2[4])
            except IndexError:
                # not valid alignment
                #logging.info("   NOT VALID ALIGNMENT %s" % line2)
                alignment_score2 = 0
                mapq2 = 0


            #logging.info("    %s at %s:%s with mapq %d and score %d" % (id2, l2[2], l2[3], mapq2, alignment_score2))

            # NOTE: This test does not work for the last line in sam1 (but probably not a big issue)
            if id2 != id1:
                # Done with id1, now we need to return the wanted alignment
                if sam2_best_alignment_score > alignment_score * 2:
                    print(sam2_best_alignment.strip())
                    replace_id = sam2_best_alignment.split()[0]

                    assert replace_id == id1, "Replace id %s != id1 %s" % (replace_id, id1)
                    n_changed_to_sam2 += 1
                    #logging.info("      CHANGING TO MINIMAP. Best score is %d" % sam2_best_alignment_score)
                    lowered_ids.add(id1)

                elif mapq > sam2_best_mapq and sam2_n_good_alignments > 1 and sam2_best_alignment_score / 2 > alignment_score * 0.99:
                    # Lower the mapq, sam2 has multiple good alignments and low mapq
                    l[4] = str(sam2_best_mapq)
                    print('\t'.join(l).strip())
                    #logging.info("      LOWERING MAPQ to %s" % l[4])
                    n_mapq_lowered += 1
                else:
                    # no adjustment, just print line
                    n_untouched += 1
                    print(line1.strip())

                # Done with one line in sam1, initlize scores etc for next line in sam1
                sam2_best_alignment_score = alignment_score2
                sam2_best_alignment = line2
                sam2_n_good_alignments = 1
                sam2_best_mapq = mapq2

                break  # break, and conitnue iterating sam1
            else:
                # New line in sam2 with same id, update scores
                if alignment_score2 > sam2_best_alignment_score:
                    sam2_best_alignment_score = alignment_score2
                    sam2_best_alignment = line2

                sam2_n_good_alignments +=1

                if mapq2 > sam2_best_mapq:
                    sam2_best_mapq = mapq2

            #logging.info(sam2_n_good_alignments)

    print(last_line)
    logging.info("%d changed to sam 2. %d lowered mapq. %d untouched. Total: %d" %
                 (n_changed_to_sam2, n_mapq_lowered, n_untouched, n_changed_to_sam2+n_mapq_lowered+n_untouched))

    with open("changed1.pckl", "wb") as f:
        logging.info("Writing %d changed ids to file" % len(lowered_ids))
        pickle.dump(lowered_ids, f)


def select_lowest_mapq_from_two_sam_files(sam1_file_name, sam2_file_name, output_file_name=None):
    sam2_n_good_alignments = defaultdict(int)
    mapq_sam2 = defaultdict(int)
    sam2_best_alignment = {}

    lowered_ids = set()
    n_changed_to_sam2 = 0

    if output_file_name is not None:
        out_sam = pysam.AlignmentFile(output_file_name, "wh", template=pysam.AlignmentFile(sam1_file_name, "r"))

    logging.info("Correcting mapq scores. First reading sam2")
    for alignment in tqdm(read_sam(sam2_file_name, skip_supplementary=True), total=number_of_lines_in_file(sam2_file_name)):
        sam2_n_good_alignments[alignment.name] += 1
        mapq_sam2[alignment.name] = max(alignment.mapq, mapq_sam2[alignment.name])  # Keep best mapq

        if alignment.name not in sam2_best_alignment:
            sam2_best_alignment[alignment.name] = alignment
        elif alignment.score > sam2_best_alignment[alignment.name].score:
            sam2_best_alignment[alignment.name] = alignment

    n_adjusted = 0
    n_untouched = 0

    logging.info("Reading sam file 1, checking if mapqs should be lowered")
    for alignment in tqdm(read_sam(sam1_file_name, skip_supplementary=True, skip_secondary=True), total=number_of_lines_in_file(sam1_file_name)):
        if alignment.name not in mapq_sam2:
            logging.warning("Sam file 2 does not have %s which is in sam file 2." % alignment.name)

        mapq2 = mapq_sam2[alignment.name]

        # if sam2 has higher score
        if sam2_best_alignment[alignment.name].score > alignment.score * 2:
            #print("sam2 has score %d > alignment score %d" % (sam2_best_alignment[alignment.name].score, alignment.score))
            alignment = sam2_best_alignment[alignment.name]
            n_changed_to_sam2 += 1
            lowered_ids.add(alignment.name)
        # If sam2 has low mapq due to multiple good alignments, use sam2's mapq score (we trust that one)
        elif alignment.mapq > mapq2 and sam2_n_good_alignments[alignment.name] > 1 and sam2_best_alignment[alignment.name].score / 2 > alignment.score * 0.8:
            alignment.set_mapq(mapq2)
            n_adjusted += 1
        else:
            n_untouched += 1

        if output_file_name is not None:
            out_sam.write(alignment.pysam_object)
        else:
            print(alignment.pysam_object.to_string())

    with open("changed2.pckl", "wb") as f:
        logging.info("Writing %d changed ids to file" % len(lowered_ids))
        pickle.dump(lowered_ids, f)

    logging.info("%d alignments were changed to sam2" % n_changed_to_sam2)
    logging.info("%d mapqs were adjusted down" % n_adjusted)
    logging.info("%d mapqs were untouched" % n_untouched)
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
