import os, glob
import pandas as pd

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

specific_genes_fastas = glob.glob(f"{snakemake.input.input_fna_dir}/*.fasta")
tmpdir=snakemake.resources['tmpdir']

output_dir=snakemake.output.report
os.makedirs(output_dir, exist_ok=True)

existing_dbs = glob.glob(f"{snakemake.input.mmseqs_db}/*_h.index")
existing_dbs = [os.path.basename(subsp).rstrip("_h.index") for subsp in existing_dbs]

for fasta_file in specific_genes_fastas:
    query_subspecies = os.path.basename(fasta_file).rstrip("_specific_genes.fasta")

    subsp_tmp_dir = f"{tmpdir}/{query_subspecies}"
    os.makedirs(subsp_tmp_dir)
    for target_subspecies in existing_dbs:
        shell(f"mmseqs easy-search --search-type 2 {fasta_file} {snakemake.input.mmseqs_db}/{target_subspecies} {output_dir}/{query_subspecies}_vs_{target_subspecies}.m8 {subsp_tmp_dir} &> {snakemake.log[0]}")