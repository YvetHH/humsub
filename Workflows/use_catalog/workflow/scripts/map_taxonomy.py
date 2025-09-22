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

import pandas as pd

subsp_relab_file = snakemake.input.subsp_relab
taxonomy_table_file = snakemake.input.taxonomy_table
taxonomy_version = snakemake.params.taxonomy_version


def map_taxonomy(subsp_relab_file, taxonomy_table_file, taxonomy_version):
    # Load the subspecies relative abundance data
    subsp_relab = pd.read_csv(subsp_relab_file, index_col=0)
    subsp_relab.index = subsp_relab.index.astype(str).str.rjust(10, "0")
    # Load the taxonomy table
    taxonomy_table = pd.read_csv(taxonomy_table_file, sep="\t", index_col=0)
    taxonomy_table["cluster95"] = (
        taxonomy_table["cluster95"].astype(str).str.rjust(4, "0")
    )

    present_subspecies = pd.DataFrame(subsp_relab.index)
    present_subspecies.columns = ["subspecies"]
    present_subspecies["cluster95"] = present_subspecies["subspecies"].str[:4]

    # Merge the subspecies with the taxonomy table
    merged = present_subspecies.merge(taxonomy_table, on="cluster95", how="left")
    merged = merged.set_index("subspecies")

    # Filter out rows where 'taxonomy' is NaN
    merged = merged[merged[taxonomy_version].notna()]
    merged = merged[~merged.index.duplicated(keep="first")]
    # Reorder columns to match the expected output
    merged["subspecies_mapped"] = merged[taxonomy_version] + "_" + merged.index.str[4:]
    subsp_relab = subsp_relab.reset_index()
    subsp_relab["name"] = subsp_relab["name"].map(merged["subspecies_mapped"])
    subsp_relab = subsp_relab.set_index("name")
    return subsp_relab


mapped_relab = map_taxonomy(subsp_relab_file, taxonomy_table_file, taxonomy_version)
mapped_relab.to_csv(snakemake.output[0], sep="\t")
