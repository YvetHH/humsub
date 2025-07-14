import os, glob
import pandas as pd
from Bio import AlignIO
from collections import defaultdict, Counter
from multiprocessing import Pool


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

# MSA parsing functions

def read_group_mapping(file_path):
    """
    Reads a TSV file using pandas and returns a mapping of sequence names to groups.
    """
    df = pd.read_csv(file_path, sep='\t', header=None, names=['sequence', 'group'])
    df['group'] = df['group'].astype(str).str.rjust(10,'0')
    return dict(zip(df['sequence'], df['group']))

def process_fasta_header(record, name_map):
    """
    Processes the FASTA header to extract the gene/genome and map it using name_map.
    """
    fasta_header = record.id
    gene = fasta_header.split(" ")[0].strip(">")
    if "GUT_GENOME" in gene:
        genome = "_".join(gene.split("_")[:2])
    else:
        genome = "_".join(gene.split("_")[:2])
    new_id = name_map.get(genome, record.id)  # Default to original ID if not found in name_map
    record.id = new_id
    record.description = new_id
    return record


def find_variants(alignment, group_mapping):
    """
    Finds and reports the most common variant for each group at a given position if at least
    one variant at that position satisfies the condition of being present in more than 50% 
    of one group and less than 10% in the other groups.
    """
    # Organize sequences by group
    group_sequences = defaultdict(list)
    original_names = defaultdict(list)
    for record in alignment:
        processed_record = process_fasta_header(record, group_mapping)
        group = str(processed_record.id)
        group_sequences[group].append(str(processed_record.seq))  # Convert Seq to str for easier manipulation
        original_names[group].append(record.name)

    # Calculate the thresholds for reporting a variant in a group
    group_thresholds = {group: (len(sequences) // 2, max(1, len(sequences) // 10))
                        for group, sequences in group_sequences.items()}

    variants_to_report = defaultdict(list)
    alignment_length = alignment.get_alignment_length()

    # First pass to find qualifying variants
    for position in range(alignment_length):
        position_variants = {group: Counter(seq[position] for seq in sequences if seq[position] != '-')
                             for group, sequences in group_sequences.items()}
        for group, variants in position_variants.items():
            for variant, count in variants.items():
                if count > group_thresholds[group][0]:
                    # This variant qualifies, check against other groups
                    in_less_than_10_percent_other_groups = all(
                        position_variants[other_group][variant] <= group_thresholds[other_group][1]
                        for other_group in group_sequences if other_group != group
                    )
                    if in_less_than_10_percent_other_groups:
                        # This position has a qualifying variant, so we'll report for all groups
                        variants_to_report[position].append(group)
                        break

    final_variants = []
    for position, groups in variants_to_report.items():
        for group in group_sequences:
            most_common_variant, count = Counter(seq[position] for seq in group_sequences[group] if seq[position] != '-').most_common(1)[0]
            total_genes = len(group_sequences[group])
            percentage = count / total_genes  # Calculate the percentage of genes with the variant
            if group in groups or most_common_variant != '-':
                variant_info = {
                    'Position': position + 1,  # 1-based indexing
                    'Reference Sequence ID': original_names[group][0],  # Assuming the first sequence can be the reference
                    'Group': group,
                    'Variant': most_common_variant,
                    'Percentage': percentage  # Add the percentage of genes with the variant
                }
                final_variants.append(variant_info)
    if len(final_variants) > 0:
        df = pd.DataFrame(final_variants).sort_values(by='Position')
        return df
    else:
        return pd.DataFrame()




# Start script

cluster95 = snakemake.wildcards.cluster95

tmpdir = f"{snakemake.resources.tmpdir}/{cluster95}"
tmp_headers = f"{tmpdir}/headers"
tmp_faa = f"{tmpdir}/faa"
tmp_msa = f"{tmpdir}/msa"

os.makedirs(tmp_headers)
os.makedirs(tmp_msa)
os.makedirs(tmp_faa)

# Read shared_function summary file
shared_functions_file = snakemake.input.shared_functions
shared_functions_df = pd.read_csv(shared_functions_file, sep='\t')
shared_functions_df.subspecies = shared_functions_df.subspecies.astype(str).str.rjust(10,'0')

# Read subspecies mapping file
subsp_def_file = snakemake.config['subsp_def']
subsp_def = pd.read_csv(subsp_def_file,sep='\t')
subsp_def.subspecies = subsp_def.subspecies.astype(str).str.rjust(10,'0')
subsp_def['cluster95'] = subsp_def.subspecies.str[:4]
n_genomes = subsp_def.query("cluster95 == @cluster95").shape[0]
subsp_mapping= dict(zip(subsp_def['genome'],subsp_def['subspecies']))

# Read cluster attribution file
cluster_attribution_file = snakemake.input.cluster_attribution
cluster_attr = pd.read_csv(cluster_attribution_file,sep='\t',header=None,names=['cluster_rep','member'])

# Size of each cluster
cluster_size = pd.DataFrame(cluster_attr['cluster_rep'].value_counts())

# Read functional annotation of each cluster representative
annotation_file = f"{snakemake.input.annotations}/annotations.tsv"
annotation = pd.read_csv(annotation_file, sep='\t',index_col=0)
for comparison in annotation['fasta'].unique():
    annotation.index = annotation.index.str.replace(f"{comparison}_","")

rep_function_mapping = annotation['ko_id'].dropna()

cluster_size['function'] = cluster_size.index.map(rep_function_mapping)
cluster_size = cluster_size.dropna()
cluster_size = cluster_size.merge(pd.DataFrame(cluster_size['function'].value_counts()), left_on='function', right_index=True)
cluster_size.columns = ['base_cluster_size','function','function_cluster_size']

function_base_cluster = cluster_size.groupby('function').sum('base_cluster_size')['base_cluster_size']
function_base_cluster = pd.DataFrame(function_base_cluster)

function_base_cluster = function_base_cluster.merge(cluster_size[['function','function_cluster_size']].reset_index(drop=True), left_index=True, right_on='function')
function_base_cluster = function_base_cluster.drop_duplicates().set_index("function")

function_base_cluster['full_function_size'] = function_base_cluster['base_cluster_size'] + function_base_cluster['function_cluster_size']

n_functions_prefilter = function_base_cluster.index.nunique()

filtered_functions = function_base_cluster.loc[function_base_cluster['full_function_size'] <= n_genomes * snakemake.config['duplicates_allowed']].index

logging.info(f"I have filtered {n_functions_prefilter - filtered_functions.nunique()} functions based on duplicates.")

group_mapping = read_group_mapping(snakemake.config['header_subspecies_mapping'])


keep_functions = shared_functions_df.loc[shared_functions_df['percent'] > snakemake.config['msa_threshold']]['function'].value_counts().loc[shared_functions_df.loc[shared_functions_df['percent'] > snakemake.config['msa_threshold']]['function'].value_counts() > 1].index
shared_functions_df = shared_functions_df.loc[shared_functions_df['function'].isin(keep_functions)]

shared_functions_df = shared_functions_df.loc[shared_functions_df['function'].isin(filtered_functions)]


functions = shared_functions_df['function'].unique()

# Maximum number of unique 'function' values per DataFrame
N = snakemake.config['maximal_functions_per_batch']

# Get unique 'function' values and split into chunks of size N
chunks = [functions[i:i + N] for i in range(0, len(functions), N)]

os.makedirs(snakemake.output.msa, exist_ok=True)
# Create a DataFrame for each chunk and append to the list
for batch_id, batch in enumerate(chunks):
    batch_df = shared_functions_df[shared_functions_df['function'].isin(batch)]
    batch_df.to_csv(f"{snakemake.output.msa}/batch_{batch_id}.tsv",sep='\t')
    