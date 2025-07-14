import os, glob
import pandas as pd
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



def filter_func(group):
    # Count subspecies with percent > upper_threshold
    high_percent_count = sum(group['percent'] > snakemake.config['upper_threshold'])
    # Count subspecies with percent < lower_threshold
    low_percent_count = sum(group['percent'] < snakemake.config['lower_threshold'])

    # Check if one subspecies has percent > 0.8 and all others have percent < 0.2
    if high_percent_count == 1 and low_percent_count == (len(group) - 1):
        return True
    return False


cluster95 = snakemake.wildcards.cluster95

annotation_file = f"{snakemake.input.annotations}/annotations.tsv"
annotation = pd.read_csv(annotation_file, sep='\t',index_col=0)

cluster_attribution_file = snakemake.input.cluster_attribution
cluster_attr = pd.read_csv(cluster_attribution_file,sep='\t',header=None,names=['cluster_rep','member'])

subsp_def_file = snakemake.config['subsp_def']
subsp_def = pd.read_csv(subsp_def_file,sep='\t')
subsp_def.subspecies = subsp_def.subspecies.astype(str).str.rjust(10,'0')
subsp_def['cluster95'] = subsp_def.subspecies.str[:4]
subsp_mapping= dict(zip(subsp_def['genome'],subsp_def['subspecies']))

fasta_header_file = snakemake.config['fasta_header_path']
fasta_headers = pd.read_csv(fasta_header_file,sep='\t')
fasta_headers_mapping = dict(zip(fasta_headers['fasta_header'],fasta_headers['genome']))

for comparison in annotation['fasta'].unique():
    annotation.index = annotation.index.str.replace(f"{comparison}_","")

grouped_functions = defaultdict()
for function in annotation['ko_id'].unique():
    if function:
        shared_function = annotation.query("ko_id == @function").index
        grouped_functions[function] = shared_function

subsp_size = subsp_def.query("cluster95 == @cluster95").subspecies.value_counts()

function_df = pd.DataFrame()
for function in grouped_functions:
    genes = grouped_functions[function]
    all_genes = cluster_attr.loc[cluster_attr.cluster_rep.isin(genes)].member.to_list()
    headers = ["_".join(gene.split("_")[:-1]) for gene in all_genes]
    genomes = []
    for header in headers:
        if "GUT_GENOME" in header:
            genome = "_".join(header.split("_")[:-1])
            genomes.append(genome)
        else:
            genome = fasta_headers_mapping[header]
            genomes.append(genome)
    
    genomes = set(genomes)
    present_subspecies = [subsp_mapping[genome] for genome in genomes]
    
    tmp_subsp_df = pd.DataFrame(pd.Series([subsp_mapping[genome] for genome in genomes]).value_counts() / subsp_size).reset_index()
    tmp_subsp_df.columns = ['subspecies','percent']
    tmp_subsp_df['function'] = function
    function_df = pd.concat([function_df,tmp_subsp_df])

function_df = function_df.dropna(subset='function')
function_df['percent'] = function_df['percent'].fillna(0)

map_ko = annotation[['ko_id','kegg_hit']].reset_index(drop=True).drop_duplicates().dropna()
function_df = function_df.merge(map_ko, left_on='function',right_on='ko_id').drop('ko_id',axis=1)


specific_functions = function_df.groupby('function').filter(filter_func)
shared_functions = function_df.loc[~function_df['function'].isin(specific_functions['function'])]

specific_functions.to_csv(snakemake.output.specific_functions, sep ='\t', index=False)
shared_functions.to_csv(snakemake.output.shared_functions, sep ='\t', index=False)