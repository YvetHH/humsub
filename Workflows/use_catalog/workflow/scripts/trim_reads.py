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
os.makedirs(snakemake.output[0], exist_ok=True)

if snakemake.config["paired_reads"]:
    read_1 = snakemake.input[0].split("/")[-1]
    read_2 = snakemake.input[1].split("/")[-1]
else:
    read_1 = snakemake.input[0].split("/")[-1]

if snakemake.config["trim_reads"] == True:
    shell(
        f"trim-low-abund.py -C 3 -Z 18 -V -M {snakemake.params.mem} -T {snakemake.params.tmp_dir} --gzip {snakemake.input[0]} {snakemake.input[1]}"
    )
    shell(f"mv {snakemake.wildcards.sample}*.abundtrim {snakemake.output[0]}")
    shell(f"mv {snakemake.output[0]}/{read_1}.abundtrim {snakemake.output[0]}/{read_1}")
    shell(f"mv {snakemake.output[0]}/{read_2}.abundtrim {snakemake.output[0]}/{read_2}")
elif snakemake.config["trim_reads"] == False:
    if snakemake.config["paired_reads"]:
        os.symlink(snakemake.input[0], f"{snakemake.output[0]}/{read_1}")
        os.symlink(snakemake.input[1], f"{snakemake.output[0]}/{read_2}")
    else:
        os.symlink(snakemake.input[0], f"{snakemake.output[0]}/{read_1}")
else:
    logger.error(f"Value for 'trim_reads' can be either True or False!")
