import os, glob
import pandas as pd
import numpy as np

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
os.makedirs(snakemake.output.batch_dir, exist_ok=True)

paths = snakemake.config["path_table"]

paths_df = pd.read_csv(paths, sep="\t", dtype={'cluster95':'str'})
paths_df['cluster95'] = paths_df['cluster95'].str.rjust(4,'0')

clusters95 = paths_df.cluster95.unique()

for cluster95 in clusters95:
    df = paths_df.loc[paths_df['cluster95'] == cluster95]
    names = (
        df.path.str.rsplit("/", 1, expand=True)[1]
        .str.rsplit(".", 2, expand=True)[0]
        .values.tolist()
    )
    output_df = pd.DataFrame()
    output_df["name"] = names
    output_df["genome_filename"] = ""
    output_df["protein_filename"] = df["path"].values.tolist()
    output_df = output_df.drop_duplicates(subset='name')
    output_df.to_csv(f"{snakemake.output.batch_dir}/batch_{cluster95}.csv", index=False)