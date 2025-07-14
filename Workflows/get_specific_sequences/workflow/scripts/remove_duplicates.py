import os, glob
import pandas as pd
from Bio import SeqIO

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

def eliminate_duplicate_pairs(df):
    # Create a set to store pairs we've seen
    seen_pairs = set()
    
    # Create a list to store indices to drop
    indices_to_drop = []
    
    for idx, row in df.iterrows():
        current_pair = (row['qseqid'], row['sseqid'])
        reverse_pair = (row['sseqid'], row['qseqid'])
        
        # Check if the reverse of this pair has been seen before
        if reverse_pair in seen_pairs:
            indices_to_drop.append(idx)
        else:
            seen_pairs.add(current_pair)
    
    # Drop the redundant rows
    df = df.drop(indices_to_drop)
    
    return df

def get_unique_genes(df,id_threshold):
    high_identity_rows = df[df['pident'] >= threshold]
    high_identity_rows = high_identity_rows.loc[high_identity_rows['qseqid'] != high_identity_rows['sseqid']]
    high_identity_rows = eliminate_duplicate_pairs(high_identity_rows)
    
    to_elimiante = []
    for qseqid in high_identity_rows['qseqid'].unique():
        similar_seqs = high_identity_rows.query("qseqid == @qseqid")['sseqid']
        to_elimiante.extend(similar_seqs)
    
    return high_identity_rows.loc[~high_identity_rows.qseqid.isin(to_elimiante)]['qseqid'].unique()


def has_similar_id(record_id, similar_ids):
    return any(sid in record_id for sid in similar_ids)


threshold = snakemake.config['deduplication_threshold']

subsp_specific_genes_paths = glob.glob(f"{snakemake.input.specific_genes_dup}/*_rmdup.fasta")

tmpdir=snakemake.resources['tmpdir']
final_output = snakemake.output.specific_genes_dedup
os.makedirs(final_output, exist_ok=True)

for subsp_input in subsp_specific_genes_paths:
    subsp = os.path.basename(subsp_input).rstrip("_specific_genes_rmdup.fasta")

    shell(f"makeblastdb -in {subsp_input} -dbtype nucl -out {tmpdir}/{subsp}_mydb")

    report_path = f"{tmpdir}/{subsp}_blast_report.txt"
    shell(f"blastn -db {tmpdir}/{subsp}_mydb -query {subsp_input} -outfmt 6 -out {report_path}")

    blast_report = pd.read_csv(report_path,sep='\t',header=None)
    blast_report.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']

    unique_genes = get_unique_genes(blast_report, threshold)

    output_fasta = f"{final_output}/{subsp}_specific.fasta"

    with open(subsp_input, "r") as handle, open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(handle, "fasta"):
            if has_similar_id(record.id, unique_genes):
                SeqIO.write(record, output_handle, "fasta")
