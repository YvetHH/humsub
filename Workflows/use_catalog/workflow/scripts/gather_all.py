import os, glob

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
import pandas as pd

# Receive inputs from Snakemake
sample_paths = snakemake.input
output_path = snakemake.output[0]


def get_relab(sample_path):
    sample_file = os.path.basename(sample_path)
    sample = os.path.splitext(sample_file)[0]
    df = pd.read_csv(sample_path)
    if "potential_false_negative" in df.columns:
        df = df[df.potential_false_negative == False]
    df = df[df.intersect_bp >= 25000]
    if df.empty:
        return pd.DataFrame()
    df[sample] = df["f_unique_weighted"] / df["f_unique_weighted"].sum()
    df = df[["name", sample]].set_index("name")
    df = df.T.groupby(axis=1, level=0).sum().T
    return df


# Aggregate all dataframes
all_dfs = []
for path in sample_paths:
    try:
        df = get_relab(path)
        if not df.empty:
            all_dfs.append(df)
    except Exception as e:
        print(f"Failed to process {path}: {e}", file=sys.stderr)

if all_dfs:
    final_df = pd.concat(all_dfs, axis=1).fillna(0)
    final_df.index = final_df.index.astype(str).str.rjust(10, "0")
    final_df.to_csv(output_path)
else:
    print("No valid dataframes to concatenate.", file=sys.stderr)
    pd.DataFrame().to_csv(output_path)
