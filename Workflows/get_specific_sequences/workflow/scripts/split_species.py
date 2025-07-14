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

output_dir = snakemake.output[0]
os.makedirs(output_dir, exist_ok=True)

genome_paths = pd.read_csv(snakemake.input.genome_paths, sep='\t')
genome_paths['subspecies'] = genome_paths['subspecies'].astype(str).str.rjust(10,'0')
genome_paths['cluster95'] = genome_paths['cluster95'].astype(str).str.rjust(4,'0')



for cluster95 in genome_paths['cluster95'].unique():
    genomes = genome_paths.query("cluster95 == @cluster95")
    genomes.columns = ['path','genome','subspecies','cluster95']
    genomes.to_csv(f"{output_dir}/{cluster95}.tsv",sep='\t',index=False)



