import os, glob
import pandas as pd
import sourmash
from collections import defaultdict

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

cluster95 = snakemake.wildcards.cluster95
tmp_dir = snakemake.resources.tmpdir
specific_hashomes_path = f"{snakemake.config['signatures_path']}/{cluster95}"
kmer_size = snakemake.config['kmer_size']

batch = pd.read_csv(snakemake.input.batch, sep='\t')
batch['subspecies'] = batch['subspecies'].astype(str).str.rjust(10,'0')
batch['cluster95'] = batch['cluster95'].astype(str).str.rjust(4,'0')

cluster95_genomes = batch.query("cluster95 == @cluster95")

relevant_genomes = defaultdict()

os.makedirs(snakemake.output.prefetch, exist_ok=True)

for subsp in cluster95_genomes['subspecies'].unique():
    relevant_subsp_genomes = []
    subsp_input = cluster95_genomes.query("subspecies == @subsp")

    speicific_hashome_path = f"{specific_hashomes_path}/{subsp}.sig.gz"

    intersect_output_dir = f"{tmp_dir}/{subsp}"
    os.makedirs(intersect_output_dir, exist_ok=True)
    try:
        specific_hashome = sourmash.signature.load_one_signature(speicific_hashome_path, ksize=kmer_size)
    except:
        continue

    specific_hashes = set([hash1 for hash1 in specific_hashome.minhash.hashes])
    remaining_hashes = len(specific_hashes)
    for genome_path in subsp_input['path']:
        genome_name = subsp_input.query("path == @genome_path")['genome'].iloc[0]
        genome_sketch_path = f"{tmp_dir}/{genome_name}.sig.gz"
        shell(f"sourmash sketch dna -p k={kmer_size},scaled=1000,noabund -o {genome_sketch_path} {genome_path}")

        shell(f"sourmash signature intersect {speicific_hashome_path} {genome_sketch_path} -q -o {intersect_output_dir}/intersect_{genome_name}.sig")


        intersection = sourmash.signature.load_one_signature(f"{intersect_output_dir}/intersect_{genome_name}.sig", ksize=51)
        intersection_hashes = set([hash1 for hash1 in intersection.minhash.hashes])

        specific_hashes = specific_hashes - intersection_hashes
        current_remaining_hashes = len(specific_hashes)

        if remaining_hashes > current_remaining_hashes:

            relevant_subsp_genomes.append(genome_name)
            specific_hashes = specific_hashes - intersection_hashes
            logging.info(f"I have found {len(intersection_hashes)}! {len(specific_hashes)} remaining!")

        remaining_hashes = current_remaining_hashes

        if len(specific_hashes) == 0:
            relevant_genomes[subsp] = relevant_subsp_genomes
            break
        


genome_df = pd.DataFrame.from_dict(relevant_genomes, orient='index')
genome_df = genome_df.T.melt(var_name='subspecies', value_name='genome')
genome_df = genome_df.dropna()

genome_path_mapping = cluster95_genomes[['genome','path']].set_index('genome')['path']

genome_df['path'] = genome_df['genome'].map(genome_path_mapping)

for subsp in genome_df['subspecies'].astype(str).unique():
    subsp_genomes = genome_df.query('subspecies == @subsp')['path']
    subsp_genomes = subsp_genomes.drop_duplicates()
    subsp_genomes.to_csv(f"{snakemake.output.prefetch}/{subsp}.txt",header=None,index=None)
