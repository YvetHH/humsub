import os, glob
import pandas as pd
from Bio import AlignIO, SeqIO
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

def parse_fasta_header(record):
    """
    Simplified processing based on the provided examples in the group_mapping explanation.
    This function will extract the base ID and map it to the format used in the group_mapping.
    """
    fasta_header = record.id
    parts = fasta_header.split("_")
    if "GUT_GENOME" in fasta_header:
        # Handles GUT_GENOME format
        base_id = f"{parts[0]}{parts[1]}"
    else:
        # Handles other formats like NZ_CP011307.1
        base_id = parts[0]
    return base_id


def dereplicate_fasta_by_group(fasta_path, group_mapping, output_path="clean_fasta.faa"):
    # Load sequences from fasta file
    sequences = list(SeqIO.parse(fasta_path, "fasta"))
    
    # Preprocess sequences and group them
    group_to_seqs = {}
    for record in sequences:
        base_id = parse_fasta_header(record)
        for key, value in group_mapping.items():
            if base_id in key:
                group = value
                if group not in group_to_seqs:
                    group_to_seqs[group] = []
                group_to_seqs[group].append(record)
                break
    
    # Deduplicate sequences within each group
    deduped_seqs = []
    seq_mapping = {}  # Map of original seq_id to the representative seq_id
    seen_seqs = {}  # Map of sequence str to seq_id for deduplication
    for group, records in group_to_seqs.items():
        for record in records:
            seq_str = str(record.seq)
            if seq_str not in seen_seqs:
                seen_seqs[seq_str] = record.id
                deduped_seqs.append(record)
            seq_mapping[record.id] = seen_seqs[seq_str]
    
    # Save the dereplicated fasta file
    SeqIO.write(deduped_seqs, output_path, "fasta")
    
    return seq_mapping

def fill_out_msa_with_duplicates(msa_path, seq_mapping, output_path="final_msa.afa"):
    # Load the MSA
    alignment = list(SeqIO.parse(msa_path, "fasta"))
    
    # Create a new alignment including duplicates
    final_alignment = []
    for record in alignment:
        # Find all original seq_ids that map to this sequence
        original_ids = [seq_id for seq_id, mapped_id in seq_mapping.items() if mapped_id == record.id]
        for orig_id in original_ids:
            new_record = record[:]
            new_record.id = orig_id
            final_alignment.append(new_record)
    
    # Save the final MSA output
    SeqIO.write(final_alignment, output_path, "fasta")

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

    #Calculate the thresholds for reporting a variant in a group
    group_thresholds = {group: (len(sequences) // 2, max(1, len(sequences) // 10))
                        for group, sequences in group_sequences.items()}

    #group_thresholds = {group: (len(sequences) * 0.2, max(1, len(sequences) // 10)) for group, sequences in group_sequences.items()}

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
tmp_clean_faa = f"{tmpdir}/clean_faa"
tmp_msa = f"{tmpdir}/msa"
tmp_msa_short = f"{tmp_msa}/short"

os.makedirs(tmp_headers)
os.makedirs(tmp_faa)
os.makedirs(tmp_clean_faa)
os.makedirs(tmp_msa)
os.makedirs(tmp_msa_short)

os.makedirs("output/msa/batch_output/{cluster95}",exist_ok=True)

shared_functions_file = snakemake.input.shared_functions
shared_functions_df = pd.read_csv(shared_functions_file, sep='\t')
shared_functions_df.subspecies = shared_functions_df.subspecies.astype(str).str.rjust(10,'0')

subsp_def_file = snakemake.config['subsp_def']
subsp_def = pd.read_csv(subsp_def_file,sep='\t')
subsp_def.subspecies = subsp_def.subspecies.astype(str).str.rjust(10,'0')
subsp_def['cluster95'] = subsp_def.subspecies.str[:4]
subsp_mapping= dict(zip(subsp_def['genome'],subsp_def['subspecies']))

cluster_attribution_file = snakemake.input.cluster_attribution
cluster_attr = pd.read_csv(cluster_attribution_file,sep='\t',header=None,names=['cluster_rep','member'])

annotation_file = f"{snakemake.input.annotations}/annotations.tsv"
annotation = pd.read_csv(annotation_file, sep='\t',index_col=0)
for comparison in annotation['fasta'].unique():
    annotation.index = annotation.index.str.replace(f"{comparison}_","")

group_mapping = read_group_mapping(snakemake.config['header_subspecies_mapping'])


keep_functions = shared_functions_df.loc[shared_functions_df['percent'] > snakemake.config['msa_threshold']]['function'].value_counts().loc[shared_functions_df.loc[shared_functions_df['percent'] > snakemake.config['msa_threshold']]['function'].value_counts() > 1].index
shared_functions_df = shared_functions_df.loc[shared_functions_df['function'].isin(keep_functions)]

functions = shared_functions_df['function'].unique()


def process_function(function):
    reps = annotation.query("ko_id == @function").index.to_list()
    genes = cluster_attr.loc[cluster_attr['cluster_rep'].isin(reps)].member
    if len(genes) > 1:
        genes.to_csv(f"{tmp_headers}/{function}.txt", index=False, header=False)

        shell(f"seqkit grep -f {tmp_headers}/{function}.txt {snakemake.input.protein_seqs} > {tmp_faa}/{function}.faa")

        seq_mapping = dereplicate_fasta_by_group(f"{tmp_faa}/{function}.faa", group_mapping, f"{tmp_clean_faa}/{function}.faa")

        shell(f"muscle -align {tmp_clean_faa}/{function}.faa  -output {tmp_msa}/short/{function}.afa")
        fill_out_msa_with_duplicates(f"{tmp_msa}/short/{function}.afa", seq_mapping, f"{tmp_msa}/{function}.afa")
        
        alignment = AlignIO.read(f"{tmp_msa}/{function}.afa", 'fasta')
        variants = find_variants(alignment, group_mapping)
        if type(variants) == pd.core.frame.DataFrame:
            variants['function'] = function
            return variants

def process_function_parallel(function):
    try:
        return process_function(function)
    except Exception as e:
        logging.info(f"Error processing function {function}: {str(e)}")
        return pd.DataFrame()  # Return an empty DataFrame in case of an error


if __name__ == "__main__":
    # Define the number of processes you want to use for parallel execution
    num_processes = snakemake.threads  # You can adjust this number according to your needs

    # Initialize a Pool of worker processes
    with Pool(num_processes) as pool:
        # Execute process_function_parallel for each function in parallel
        results = pool.map(process_function_parallel, functions)

    # Concatenate the results into the output_table DataFrame
    output_table = pd.concat(results, axis=0)

    # Save the output DataFrame to a CSV file
    output_table.to_csv(snakemake.output.msa, sep='\t', index=False)