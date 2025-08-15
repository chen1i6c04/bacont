import csv
import gzip
import subprocess
from loguru import logger
from Bio import SeqIO


def syscall(cmd, stdout=False, stderr=False, raise_error=True):
    if stdout:
        stdout_str = subprocess.PIPE if isinstance(stdout, bool) else stdout
    else:
        stdout_str = None
    if stderr:
        stderr_str = subprocess.PIPE if isinstance(stderr, bool) else stderr
    else:
        stderr_str = None
    shell, executable = (True, '/bin/bash') if isinstance(cmd, str) else (False, None)
    child_process = subprocess.run(
        cmd, shell=shell, stdout=stdout_str, stderr=stderr_str, universal_newlines=True, executable=executable
    )
    if child_process.returncode and raise_error:
        logger.error(f"Command '{cmd}' failed")
    elif child_process.returncode and raise_error is False:
        logger.debug(f"Command '{cmd}' failed")
    return child_process


def validate_input(file):
    def is_gzip(f):
        with open(f, 'rb') as f:
            signature = f.read(2)
        return signature == b'\x1f\x8b'
    open_function = gzip.open if is_gzip(file) else open
    with open_function(file, 'rt') as handle:
        try:
            if any(SeqIO.parse(handle, 'fastq')):
                logger.info(f"FASTQ {file} checked")
                return 'fastq'
            else:
                raise ValueError
        except ValueError:
            handle.seek(0)  # Reset the file pointer to the beginning
            if any(SeqIO.parse(handle, 'fasta')):
                logger.info(f"FASTA {file} checked")
                return 'fasta'
            else:
                logger.error(f"Input file {file} is neither FASTA nor FASTQ format.")


def validate_medaka_model(model: str):
    output = syscall(['medaka', 'tools', 'list_models'], stdout=True).stdout
    available_models = set(output.splitlines()[0].replace('Available: ', '').split(', '))
    if not (model in available_models):
        logger.error(f"Medaka model {model} unavailable")


def sort_assembly(infile, outfile):
    records = SeqIO.parse(infile, 'fasta')
    records = sorted(records, key=lambda record: len(record.seq), reverse=True)
    SeqIO.write(records, outfile, 'fasta')


def annotate_assembly(assembly, out_assembly, assembly_info):
    coverage_dict = dict()
    circular_dict = dict()
    with open(assembly_info) as handle:
        spamreader = csv.reader(handle, delimiter='\t')
        next(spamreader)  # ignore header
        for row in spamreader:
            coverage_dict[row[0]] = row[2]
            circular_dict[row[0]] = row[3]
    annotated_records = []
    for record in SeqIO.parse(assembly, 'fasta'):
        length = len(record)
        coverage = coverage_dict[record.id]
        circular = circular_dict[record.id]
        record.description = f"length={length} coverage={coverage} circular={circular}"
        annotated_records.append(record)
    annotated_records = sorted(annotated_records, key=lambda record: len(record.seq), reverse=True)
    SeqIO.write(annotated_records, out_assembly, 'fasta')


def fasta2fastq(in_fasta, out_fastq):
    def is_gzip(f):
        with open(f, 'rb') as f:
            signature = f.read(2)
        return signature == b'\x1f\x8b'
    open_function = gzip.open if is_gzip(in_fasta) else open
    with open(out_fastq, 'w') as out_handle:
        for record in SeqIO.parse(open_function(in_fasta, 'rt'), 'fasta'):
            record.letter_annotations['phred_quality'] = [40] * len(record)
            out_handle.write(record.format('fastq'))