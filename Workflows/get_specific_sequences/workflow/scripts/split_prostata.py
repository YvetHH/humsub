import numpy as np
import pandas as pd
from Bio import SeqIO

from snakemake.shell import shell
import logging, traceback
import os, glob
from collections import defaultdict
from multiprocessing import Pool

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



def read_sequence_from_fasta(file_path, sequence_id):
    for record in SeqIO.parse(file_path, "fasta"):
        if record.id == sequence_id:
            return str(record.seq)
    return None

def calculate_median_seq_len(file_path, sequences_id):
    seq_len_list = []
    for seq_id in sequences_id:
        seq = read_sequence_from_fasta(file_path, seq_id)
        len_seq = len(seq)
        seq_len_list.append(len_seq)
    return np.median(seq_len_list)

def split_dataframe(df, batch_size):
    """
    Splits a DataFrame into smaller DataFrames with a maximum length of batch_size.

    :param df: Pandas DataFrame to be split.
    :param batch_size: Maximum number of rows in each batch.
    :return: List of DataFrames.
    """
    batches = [df[i:i+batch_size] for i in range(0, df.shape[0], batch_size)]
    return batches


# annotation_file = f"{snakemake.input.annotation}/annotations.tsv"
# annotation = pd.read_csv(annotation_file,sep='\t',index_col=0)
# for comparison in annotation['fasta'].unique():
#     annotation.index = annotation.index.str.replace(f"{comparison}_","")

# reps = annotation['ko_id']

# cluster_attribution_file = snakemake.input.cluster_attribution
# cluster_attr = pd.read_csv(cluster_attribution_file,sep='\t',header=None,names=['cluster_rep','member'])
# cluster_attr = cluster_attr.set_index('cluster_rep')['member']

# mapped_functions = {}
# for cluster_id, members in cluster_attr.groupby(level=0):
#     function = reps[cluster_id]
#     for member in members.values:
#         mapped_functions[member] = function

# # Convert the mapping to a Pandas Series
# seq_id_function_mapping = pd.DataFrame.from_dict(mapped_functions, orient='index').reset_index()
# seq_id_function_mapping.columns = ['seq_id','function']

msa_output = snakemake.input.msa_output
msa = pd.read_csv(msa_output,sep='\t')
#msa = msa.drop('function',axis=1)
msa['Group'] = msa['Group'].astype(str).str.rjust(10,'0')
#msa = msa.merge(seq_id_function_mapping, left_on='Reference Sequence ID', right_on='seq_id').drop("seq_id",axis=1)
msa['Position'] = msa['Position'].astype(int)

existing_subspecies = msa['Group'].unique()

max_number_of_variants = int(snakemake.config['max_number_of_variants'])

#keep_functions = []



def process_function(function):
    function_df = msa.query(f'function == "{function}"')
    n_variants = function_df['Position'].nunique()

    seqs = function_df['Reference Sequence ID'].tolist()
    # Assuming calculate_median_seq_len and necessary data are accessible
    median_len = calculate_median_seq_len(snakemake.input.protein_seqs, seqs)

    if n_variants / median_len < snakemake.config['diverse_threshold']:
        return function
    else:
        return None

# Function to parallelize the processing
def parallelize_processing(msa):
    unique_functions = msa['function'].unique()
    
    # Create a pool of workers
    with Pool(processes=snakemake.threads) as pool:  # Adjust the number of processes as needed
        results = pool.map(process_function, unique_functions)
    
    # Filter None results
    keep_functions = [function for function in results if function is not None]
    
    return keep_functions

# Assuming msa is defined
#keep_functions = parallelize_processing(msa)


# for function in msa['function'].unique():
#     function_df = msa.query('function == @function')
#     n_variants = function_df['Position'].nunique()

#     seqs = function_df['Reference Sequence ID'].tolist()
#     median_len = calculate_median_seq_len(snakemake.input.protein_seqs, seqs)

#     if n_variants/median_len < snakemake.config['diverse_threshold']:
#         keep_functions.append(function)

#keep_functions_df = msa.loc[msa['function'].isin(keep_functions)]

df_list = split_dataframe(msa, max_number_of_variants)

os.makedirs(snakemake.output.tmp_output, exist_ok=True)

for idx, df in enumerate(df_list):
    df.to_csv(f"{snakemake.output.tmp_output}/batch_{idx}.tsv",sep='\t')
