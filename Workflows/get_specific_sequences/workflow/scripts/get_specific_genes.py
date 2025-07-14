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

subsp_specific_kmers_paths = glob.glob(f"{snakemake.input.specific_kmers}/*_kmers.csv")
tmpdir = snakemake.resources.tmpdir
prodigal_output_path = snakemake.config['prodigal_output']

prodigal_output = pd.read_csv(prodigal_output_path,sep='\t')

final_output = snakemake.output.specific_genes
os.makedirs(final_output, exist_ok=True)

for subsp_input in subsp_specific_kmers_paths:
    subsp = os.path.basename(subsp_input).rstrip("_kmers.csv")
    os.makedirs(f"{tmpdir}/tmp_seqs/{subsp}")

    kmer_input = pd.read_csv(subsp_input)[['sequence_name','kmer']]
    kmer_input = kmer_input.drop_duplicates(subset='kmer')
    kmer_input['sequence_name'] = kmer_input['sequence_name'].str.rsplit("_",n=1,expand=True)[0]


    relevant_genomes_paths = f"{snakemake.input.prefetch}/{subsp}.txt"
    relevant_genomes_df = pd.read_csv(relevant_genomes_paths,header=None)
    relevant_genomes_df.columns = ['genome_fasta']
    relevant_genomes_df['genome'] = relevant_genomes_df['genome_fasta'].str.split("/",expand=True).iloc[:,-1].str.rsplit(".",n=2,expand=True).iloc[:,0]

    kmer_input['kmer'].to_csv(f"{tmpdir}/{subsp}_kmers.txt",index=False,header=False)

    for genome in relevant_genomes_df['genome'].unique():

        prodigal_output_genome = prodigal_output.query("genome == @genome")['path'].iloc[0]

        if prodigal_output_genome:
            shell(f"cat {prodigal_output_genome} | seqkit grep -s -i -f {tmpdir}/{subsp}_kmers.txt > {tmpdir}/tmp_seqs/{subsp}/{genome}_kmers.fasta")
        else:
            logging.warning(f"No prodigal output file for {genome}!")

    
    shell(f"cat {tmpdir}/tmp_seqs/{subsp}/*_kmers.fasta > {final_output}/{subsp}_specific_genes.fasta")