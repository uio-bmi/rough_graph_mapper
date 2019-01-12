import subprocess
import logging
from tqdm import tqdm

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
                continue

            chrom = line.split()[2]
            if chrom not in chromosome_index:
                n_not_matched += 1
                continue

            outfiles[chrom].writelines([line])

    for chrom in chromosomes:
        outfiles[chrom].close()

    logging.info("Done splitting by chromosome. %d lines did not match any chromosomes or were unaligned." % n_not_matched)


def number_of_lines_in_file(file_name):
    return sum(1 for _ in open(file_name, 'rb'))


def run_minimap2(reference_file_name, fasta_file_name, output_file_name, arguments=""):
    pass


def select_lowest_mapq_from_two_sam_files(sam1_file_name, sam2_file_name, output_file_name):
    pass


def run_hybrid_between_bwa_and_minimap(reference_file_name, fasta_file_name,
                                       output_file_name, bwa_arguments="", minimap_arguments=""):
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


def read_fasta_file(fasta_file_name):
    pass
