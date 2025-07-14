import pandas as pd

from snakemake.shell import shell
import logging, traceback
import os, glob
from collections import defaultdict

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

output_df = pd.DataFrame()
for batch_input in snakemake.input:
    if os.path.getsize(batch_input) > 0:
        try:
            tmp_df = pd.read_csv(batch_input, sep='\t')
            output_df = pd.concat([output_df, tmp_df])
        except:
            logging.info(f"The batch input {batch_input} is empty, continuing . . .")

output_df['Position'] = output_df['Position'].astype(int)
output_df.to_csv(snakemake.output.msa, sep ='\t',index=False)