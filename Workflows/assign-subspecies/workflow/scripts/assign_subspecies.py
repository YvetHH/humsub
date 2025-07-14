import pandas as pd
import numpy as np
import glob, os
import logging, traceback
from snakemake.shell import shell

import concurrent.futures

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

# Script starts here

tmpdir = snakemake.resources['tmpdir']

batch_file = pd.read_csv(f"output/batches/{snakemake.wildcards.batch}.csv")
input_signatures = snakemake.input.sketch

genomes = batch_file['name'].to_list()

def process_genome(genome):
    shell(f"sourmash signature extract {input_signatures} --name {genome} > {tmpdir}/{genome}.sig")
    shell(f"sourmash search {tmpdir}/{genome}.sig {snakemake.input.index} --containment --threshold 0 --output {tmpdir}/{genome}.csv")

with concurrent.futures.ThreadPoolExecutor(max_workers=snakemake.threads) as executor:
    executor.map(process_genome, genomes)

shell(f"csvtk concat -I -T -o {snakemake.output.search_result} {tmpdir}/*.csv")