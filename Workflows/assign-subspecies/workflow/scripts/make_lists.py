import pandas as pd
import numpy as np
import glob, os
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

os.makedirs(snakemake.output.batch_dir, exist_ok=True)

batch_size = snakemake.config["batch_size"]
if snakemake.config['input_dir'] != "":
    argument1 = f"{snakemake.config['input_dir']}/*{snakemake.config['format1']}"
    argument2 = f"{snakemake.config['input_dir']}/*{snakemake.config['format2']}"
    genomes = glob.glob(argument1)
    genomes.extend(glob.glob(argument2))

    genomes = [os.path.abspath(genome) for genome in genomes] # Get abs paths for all input genomes
    paths = pd.DataFrame(genomes, columns = ['path'])
else:   
    paths_file = snakemake.config["input_text_file"]
    paths = pd.read_csv(paths_file,names=['path'])

n_chunks = paths.shape[0] // batch_size
if paths.shape[0] % batch_size != 0:
    n_chunks += 1

for chunk_id, chunk in enumerate(np.array_split(paths, n_chunks), start=1):
    names = (
        chunk.path.str.rsplit("/", 1, expand=True)[1]
        .str.rstrip(".gz")
        .str.rsplit(".", 1, expand=True)[0]
        .values.tolist()
    )
    output_df = pd.DataFrame()
    output_df["name"] = names
    output_df["genome_filename"] = chunk["path"].values.tolist()
    output_df["protein_filename"] = ""
    output_df.to_csv(f"{snakemake.output.batch_dir}/{chunk_id}.csv", index=False)
    logging.info(f"Created batch file {snakemake.output.batch_dir}/{chunk_id}.csv")