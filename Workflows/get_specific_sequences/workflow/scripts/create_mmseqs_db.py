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

def concat_fasta_files(path_list, output_path):
    with open(output_path, 'w') as out:
        for file_path in path_list:
            file_path = file_path.strip()
            with open(file_path, 'r') as fasta:
                out.write(fasta.read() + "\n")


os.makedirs(snakemake.output.db, exist_ok=True)


prodigal_output_path = snakemake.config['prodigal_output']
prodigal_output = pd.read_csv(prodigal_output_path,sep='\t')

all_genomes = pd.read_csv(snakemake.input.genome_paths, sep='\t')
all_genomes['subspecies'] = all_genomes['subspecies'].astype(str).str.rjust(10,'0')
all_genomes['cluster95'] = all_genomes['cluster95'].astype(str).str.rjust(4,'0')

prodigal_output = prodigal_output.merge(all_genomes[['cluster95','subspecies','genome']], left_on='genome', right_on='genome')

cluster95 = snakemake.wildcards.cluster95
tmp_dir = f"{snakemake.resources['tmpdir']}/{cluster95}/mmseqs_db"
os.makedirs(tmp_dir, exist_ok=True)


cluster95_genes = prodigal_output.query("cluster95 == @cluster95")


output_dir = "output/genes"
os.makedirs(output_dir, exist_ok=True)

concat_fasta_files(cluster95_genes['path'].to_list(), f"{output_dir}/{cluster95}_prodigal.fna")

shell(f"seqkit translate --trim {output_dir}/{cluster95}_prodigal.fna > {snakemake.output.protein_output}")
shell(f"mmseqs createdb {snakemake.output.protein_output} {snakemake.output.db}/{cluster95}_full 2>> {snakemake.log[0]}")