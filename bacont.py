import os
import re
import csv
import sys
import json
import shutil
import argparse
import subprocess
from tempfile import TemporaryDirectory, NamedTemporaryFile
from loguru import logger
from src.utils import (syscall, validate_input, validate_medaka_model, sort_assembly, annotate_assembly, fasta2fastq)
from src.database import check_database_installation


__location__ = os.path.dirname(os.path.abspath(__file__))


def filter_reads(input_file:str, output_file:str, keep_percent:int, min_length:int, min_quality:int, num_threads:int):
    logger.info(f"Keep {keep_percent} percentage of the best reads.")
    logger.info(f"Remove reads with length less than {min_length} or quality score less than {min_quality}.")
    cmd = ['filtlong', input_file]
    if keep_percent:
        cmd += ['--keep_percent', str(keep_percent)]
    if min_length:
        cmd += ['--min_length', str(min_length)]
    if min_quality:
        cmd += ['--min_mean_q', str(min_quality)]
    cmd = ' '.join(cmd) + ' | ' + f'pigz -6 -p {num_threads} > {output_file}'
    syscall(cmd)


def remove_adapters(input_file, output_file, threads=1):
    logger.info("Remove adapter sequence")
    cmd = f"porechop_abi -i {input_file} -o {output_file} --discard_middle --check_reads 1000 -t {threads}"
    syscall(cmd)


def preprocess_reads(reads, outdir, disable_trimming, excluded_target=None, keep_percent=95,
                     min_length=0, min_quality=0, threads=1):
    input_format = validate_input(reads)
    if input_format == 'fasta':
        converted_reads = os.path.join(outdir, 'READS.fastq')
        fasta2fastq(reads, converted_reads)
        reads = converted_reads
    if disable_trimming is False:
        trimmed_reads = os.path.join(outdir, f'READS.trim.fastq.gz')
        remove_adapters(reads, trimmed_reads, threads)
        reads = trimmed_reads
    if excluded_target:
        cleaned_reads = os.path.join(outdir, f'READS.clean.fastq.gz')
        exclude_target(input_fastq=reads, output_fastq=cleaned_reads, target=excluded_target, threads=threads)
        reads = cleaned_reads
    if any((keep_percent, min_length, min_quality)):
        filtered_reads = os.path.join(outdir, f'READS.filter.fastq.gz')
        filter_reads(reads, filtered_reads, keep_percent=keep_percent, min_length=min_length,
                     min_quality=min_quality, num_threads=threads)
        reads = filtered_reads
    report_file = os.path.join(outdir, 'nanoq.json')
    syscall(f"nanoq -j -s -f -r {report_file} -i {reads}")
    return reads


def run_dnaapler(flye_output, outdir, threads):
    assembly = os.path.join(flye_output, 'assembly.fasta')
    assembly_info = os.path.join(flye_output, 'assembly_info.txt')
    with open(assembly_info) as handle, NamedTemporaryFile('w') as tmpfile:
        spamreader = csv.reader(handle, delimiter='\t')
        next(spamreader)  # ignore header
        for row in spamreader:
            if row[3] == 'N':
                tmpfile.write(row[0] + '\n')
        tmpfile.flush()
        cmd = f"dnaapler all -i {assembly} -o {outdir} -t {threads} --ignore {tmpfile.name}"
        logger.info(f"Running : {cmd}")
        syscall(cmd)
    return os.path.join(outdir, 'dnaapler_reoriented.fasta')



def estimate_genome_size(reads, num_threads):
    current_location = os.path.dirname(os.path.abspath(__file__))
    program = os.path.join(current_location, 'bin', 'genome_size_raven.sh')
    with TemporaryDirectory() as tmpdir:
        tmpfile = os.path.join(tmpdir, os.path.basename(reads))
        syscall(f"rasusa reads -b 100M {reads} -o {tmpfile}")
        cmd = f"{program} {tmpfile} {num_threads}"
        output = syscall(cmd, stdout=True).stdout
    genome_size = int(output.strip())
    logger.info(f"Estimated genome size was {genome_size}bp.")
    return genome_size


def parse_genome_size(pattern):
    unit_map = {'K': 1e3, 'M': 1e6, 'G': 1e9, 'k': 1e3, 'm': 1e6, 'g': 1e9}
    result = re.fullmatch(r'^([\d.]+)([KMGkmg])', pattern)
    if result is None:
        logger.error(f"Couldn't parse {pattern}")
    else:
        value, unit = result.groups()
        return int(float(value) * unit_map[unit])


def run_flye(reads, outdir, flye_mode, threads, meta):
    cmd = f'flye --{flye_mode} {reads} -o {outdir} -t {threads}'
    if meta:
        cmd += ' --meta --keep-haplotypes --no-alt-contigs'
    logger.info(f"Running : {cmd}")
    syscall(cmd)


def run_medaka(assembly, reads, outdir, model, threads):
    cmd = f"medaka_consensus -i {reads} -d {assembly} -o {outdir} -m {model} -t {threads} -f"
    logger.info(f"Running : {cmd}")
    syscall(cmd)


def exclude_target(input_fastq, output_fastq, target, threads):
    cmd = (f"minimap2 --secondary=no -p 0.95 -A 2 -B 4 -O 6,24 -E 3,1 -t {threads} -ax map-ont {target} {input_fastq} | "
           f"samtools sort -@ {threads} -O BAM - | "
           f"samtools view -@ {threads} -f 4 -O BAM - | "
           f"samtools fastq -@ {threads} - | "
           f"pigz -9 -p {threads} > {output_fastq}")
    syscall(cmd)


def sequence_classify(fastq, kraken2_database, outdir, threads):
    kraken2_output = os.devnull
    kraken2_report = os.path.join(outdir, 'kraken2.txt')
    bracken_report = os.path.join(outdir, 'bracken.txt')

    # kraken2_cmd = (
    #     f"kraken2 --db {kraken2_database} --threads {threads} --memory-mapping --confidence 0.01 "
    #     f"--output {kraken2_output} --report {kraken2_report} {fastq}"
    # )
    ncbi_plasmid = os.path.join(__location__, "db", "ncbi_plasmids.mmi")
    kraken2_cmd = (f'minimap2 -ax map-ont --secondary=no -p 0.9 -A 2 -B 4 -O 4,16 -E 2,1 -t {threads} {ncbi_plasmid} {fastq} | '
                   f'samtools view -@ {threads} -f 4 | '
                   f'samtools fastq -@ {threads} | '
                   f'kraken2 --db {kraken2_database} --memory-mapping --confidence 0.05 '
                   f'--output {kraken2_output} --report {kraken2_report} /dev/fd/0 -t {threads}')
    syscall(kraken2_cmd)
    bracken_cmd = f"bracken -i {kraken2_report} -d {kraken2_database} -w /dev/null -o {bracken_report} -l S"
    try:
        subprocess.run(bracken_cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        bracken_cmd = f"bracken -i {kraken2_report} -d {kraken2_database} -w /dev/null -o {bracken_report} -l G"
        subprocess.run(bracken_cmd, shell=True, check=True)


def estimate_genome_completeness(input_file, output_dir, database, num_threads=1):
    specific_pattern = re.compile("^short_summary.specific.*.busco.json$")
    generic_pattern = re.compile("^short_summary.generic.*.busco.json$")
    with TemporaryDirectory() as tmpdir:
        cmd = f"busco -o busco --out_path {tmpdir} -m geno --offline --download_path {database} " \
              f"-i {input_file} -c {num_threads} --auto-lineage-prok"
        logger.info(f"Running: {cmd}")
        syscall(cmd, raise_error=False)
        out_path = os.path.join(tmpdir, "busco")
        for f in os.listdir(out_path):
            match = specific_pattern.fullmatch(f)
            if match:
                shutil.copyfile(
                    os.path.join(out_path, f),
                    os.path.join(output_dir, 'busco_specific_summary.json')
                )
            match = generic_pattern.fullmatch(f)
            if match:
                shutil.copyfile(
                    os.path.join(out_path, f),
                    os.path.join(output_dir, 'busco_generic_summary.json')
                )


def species_identify(input_file, output_dir, database, num_threads=1):
    cmd = f"skani search -q {input_file} -d {database} -o {output_dir} -n 10 -t {num_threads}"
    syscall(cmd, raise_error=False)


def check_dependency():
    version = {
        'rasusa': 'rasusa --version',
        'nanoq': 'nanoq --version',
        'filtlong': 'filtlong --version',
        'flye': 'flye -v',
        'dnaapler': 'dnaapler -V',
        'medaka': 'medaka --version',
        "kraken2": "kraken2 -v",
        "bracken": "bracken -v",
        "raven": "raven --version",
        "busco": "busco -v",
        "skani": "skani -V",
        "porechop_abi": "porechop_abi --version",
    }
    for program_name, cmd in version.items():
        child_process = syscall(cmd, stdout=True)
        if child_process.returncode:
            logger.error(f"Could not determine version of {program_name}")
            sys.exit("Abort")
        else:
            version = re.search('\d+(\.\d+)*', child_process.stdout).group()
            logger.info(f"Using {program_name:12} | {version}")


def setting_logger(outdir):
    fmt = "{time:YYYY-MM-DD HH:mm:ss} [{level}] {message}"
    logfile = os.path.join(outdir, 'bacont.log')
    logger.add(logfile, format=fmt, level='INFO', mode='w')
    logger.add(sys.stderr, format=fmt, level='ERROR')
    logger.add(lambda _: sys.exit(1), level="ERROR")


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required")
    required.add_argument("reads", help="ONT reads")
    required.add_argument("-o", "--outdir", required=True, help="Directory to store all the resulting files")

    optional = parser.add_argument_group("Optional")
    optional.add_argument("-t", "--threads", default=8, type=int, required=False, help="number of threads [Default 8]")
    optional.add_argument('-g', '--genome-size', default=None,
                          help='Estimated genome size eg. 3.2M, if empty, will bell AUTO.')
    optional.add_argument('-x', '--depth', default=100, type=int,
                          help='Sub-sample reads to this depth. Disable with --depth 0')
    optional.add_argument('-l', '--min-length', metavar='', default=0, type=int)
    optional.add_argument('-q', '--min-quality', metavar='', default=0, type=int)
    optional.add_argument('-p', '--keep-percent', metavar='', default=95, type=int,
                          help='Keep only this percentage of the best reads. '
                               'When this parameter is used, -l and -q will bell disabled. '
                               'Disable with -p 0')
    optional.add_argument('-a', '--disable_adapter_trimming', action='store_true',
                          help='Adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled.')
    optional.add_argument('-m', '--medaka_model', default='r1041_e82_400bps_bacterial_methylation',
                          help='The model to be used by Medaka (default: r1041_e82_400bps_bacterial_methylation)')
    optional.add_argument('--no_medaka', action='store_true', default=False,
                          help='Disable medaka polishing.')
    optional.add_argument('--exclude',
                          help='Fasta file or Minimap2 index of target to be excluded.')
    optional.add_argument('--meta', action='store_true', help='Enable meta-genome mode of flye.')
    optional.add_argument('--flye-mode', choices=['nano-raw', 'nano-hq'], default='nano-hq', help='Flye assembly mode.')

    database = parser.add_argument_group("Database")
    database.add_argument(
        "--kraken2_db",
        default='',
        required=False,
        help="Path of kraken2/bracken database. If not provide, won't run kraken2"
    )
    database.add_argument(
        "--busco_db",
        default='',
        required=False,
        help="Path of busco database. If not provide, won't run busco"
    )
    database.add_argument(
        "--skani_db",
        default='',
        required=False,
        help="Path of skANI database. If not provide, won't run skANI for identifying species."
    )

    args = parser.parse_args()
    threads = args.threads

    os.makedirs(args.outdir, exist_ok=True)
    setting_logger(args.outdir)

    logger.info("Checking dependencies")
    check_dependency()
    validate_medaka_model(args.medaka_model)
    check_database_installation(os.path.join(__location__, 'db'))

    reads = preprocess_reads(
        args.reads, args.outdir, disable_trimming=args.disable_adapter_trimming, excluded_target=args.exclude,
        keep_percent=args.keep_percent, min_length=args.min_length, min_quality=args.min_quality, threads=args.threads
    )
    if args.kraken2_db:
        logger.info("Running Kraken2 & Bracken")
        logger.info(f"Kraken2 database is {args.kraken2_db}")
        sequence_classify(reads, args.kraken2_db, args.outdir, args.threads)

    with open(os.path.join(args.outdir, 'nanoq.json')) as handle:
        total_bases = json.load(handle)['bases']
    logger.info(f"Total bases: {total_bases}bp")

    if args.genome_size:
        gsize = parse_genome_size(args.genome_size)
        logger.info(f"Using genome size was {gsize}bp.")
    else:
        gsize = estimate_genome_size(reads, threads)

    origin_depth = int(total_bases / gsize)
    logger.info(f"Estimated sequencing depth: {origin_depth}x.")

    if args.depth:
        if origin_depth > args.depth:
            logger.info(f"Subsampling reads from {origin_depth:.0f}x to {args.depth}x.")
            sub_reads = os.path.join(args.outdir, f'READS.sub.fastq')
            syscall(f"rasusa reads -b {gsize * args.depth} {reads} -o {sub_reads}")
        else:
            logger.info("No read depth reduction requested or necessary.")
            sub_reads = reads
    else:
        sub_reads = reads

    # Genome assemble
    flye_dir = os.path.join(args.outdir, 'flye')
    run_flye(sub_reads, flye_dir, args.flye_mode, args.threads, args.meta)
    shutil.copyfile(os.path.join(flye_dir, 'assembly_graph.gfa'), os.path.join(args.outdir, 'flye_unpolished.gfa'))
    shutil.copyfile(os.path.join(flye_dir, 'assembly_info.txt'), os.path.join(args.outdir, 'flye_info.txt'))
    # Reorients assembled
    dnaapler_dir = os.path.join(args.outdir, 'dnaapler')
    reoriented_asm = run_dnaapler(flye_output=flye_dir, outdir=dnaapler_dir, threads=args.threads)
    sort_assembly(reoriented_asm, os.path.join(args.outdir, 'flye.fasta'))
    # Sequence correction
    medaka_dir = os.path.join(args.outdir, 'medaka')
    if args.no_medaka:
        logger.info("Parameter `no_medaka` is activated, skip medaka polishing.")
        final_asm = os.path.join(args.outdir, 'flye.fasta')
    else:
        final_asm = os.path.join(medaka_dir, 'consensus.fasta')
        run_medaka(reoriented_asm, reads, medaka_dir, args.medaka_model, args.threads)
    annotated_asm = os.path.join(args.outdir, 'assembly.fasta')
    annotate_assembly(
        assembly=final_asm,
        out_assembly=annotated_asm,
        assembly_info=os.path.join(flye_dir, 'assembly_info.txt'),
    )
    # Evaluation of genome completeness
    if args.busco_db:
        estimate_genome_completeness(annotated_asm, args.outdir, args.busco_db, threads)

    if args.skani_db:
        logger.info("Running skANI.")
        species_identify(annotated_asm, os.path.join(args.outdir, 'skani.txt'), args.skani_db, threads)

    for dirpath in (flye_dir, medaka_dir, dnaapler_dir):
        try:
            shutil.rmtree(dirpath)
        except FileNotFoundError:
            pass
    syscall(f"rm {os.path.join(args.outdir, 'READS*')}")
    logger.info("Done")


if __name__ == '__main__':
    main()
