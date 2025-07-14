import os, glob
import pandas as pd
import concurrent.futures

from snakemake.shell import shell
import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logging.captureWarnings(True)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logging.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


# Install exception handler
sys.excepthook = handle_exception

# Start script

tmpdir = snakemake.resources.tmpdir
kmer_size = snakemake.config['kmer_size']

subspecies = glob.glob(f"{snakemake.input.prefetch}/*.txt")
subspecies = [subsp.split("/")[-1].rstrip(".txt") for subsp in subspecies]

os.makedirs(snakemake.output.specific_kmers, exist_ok=True)


def process_subspecies(subsp):
    # concat genomes of one subspecies
    shell(f"cat $(<{snakemake.input.prefetch}/{subsp}.txt) > {tmpdir}/{subsp}_combined.fasta")
    logging.info(f"Saved combined fasta file for subspecies {subsp}!")

    shell(f"sourmash signature kmers --check-sequence --force -k {kmer_size} --signatures {snakemake.input.specific_hashome}/{subsp}.sig.gz --save-kmers {snakemake.output.specific_kmers}/{subsp}_kmers.csv --sequences {tmpdir}/{subsp}_combined.fasta")


# Parallel execution
with concurrent.futures.ProcessPoolExecutor(max_workers = int(snakemake.threads)) as executor:
    futures = [executor.submit(process_subspecies, subsp) for subsp in subspecies]
    for future in concurrent.futures.as_completed(futures):
        future.result()  # This will raise any exceptions caught during the execution