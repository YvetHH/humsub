import os, glob
import shutil
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

# start script


import pandas as pd
from snakemake.shell import shell

import sys, os

from sklearn.model_selection import train_test_split, KFold
from scipy.cluster import hierarchy as hc
import fastcluster as fc
from scipy import spatial as sp
import numpy as np

# import utils
from tqdm import trange

from collections import defaultdict

sys.path.append(f"{snakemake.params.script_dir}/utils/")

sys.path.append(f"{snakemake.params.script_dir}/")
from prediction_strength import *

import os, glob


cluster95 = int(snakemake.wildcards.cluster95)
prediction_strength_cutoff = snakemake.config["prediction_strength_cutoff"]
lower_threshold_n_genomes = int(snakemake.config["lower_threshold_n_genomes"])
n_iterations = int(snakemake.config["n_iterations"])


input_tsv_file = snakemake.input.distance_matrix
temporary_output = f"{snakemake.config['tmp_dir']}/{snakemake.wildcards.cluster95}"
os.makedirs(temporary_output, exist_ok=True)

if os.path.getsize(input_tsv_file) > 0:
    DistMatrix = 1 - pd.read_csv(input_tsv_file)
    DistMatrix = DistMatrix.set_index(DistMatrix.columns)
    if len(DistMatrix) > 3:
        all_genomes = DistMatrix.index.to_list()

        random_state = None
        K_max = min(len(all_genomes) // 2, 15)
        clusters = range(1, K_max + 1)

        strengths = np.zeros((n_iterations * 2, len(clusters)))

        for i in trange(n_iterations):
            random_state = None

            index_train, index_test = train_test_split(
                np.arange(DistMatrix.shape[0]),
                test_size=0.5,
                shuffle=True,
                random_state=random_state,
            )

            dist_matrix_train = DistMatrix.iloc[index_train, index_train].values

            # hc (deterministic)
            linkage_train = fc.linkage(dist_matrix_train, method="single")

            dist_matrix_test = DistMatrix.iloc[index_test, index_test].values
            linkage_test = fc.linkage(dist_matrix_test, method="single")

            for j, k in enumerate(clusters):
                labels_train = hc.fcluster(linkage_train, k, criterion="maxclust")
                labels_test = hc.fcluster(linkage_test, k, criterion="maxclust")

                # HC may not be able to cluster with k clusters, in such cases we set 0 as prediction strength.

                if (labels_test.max() < k) or (labels_train.max() < k):
                    strengths[i, j] = 0
                    strengths[i + n_iterations, j] = 0

                else:
                    dist_matrix_cross = DistMatrix.iloc[index_train, index_test]
                    transfered_labels = (
                        dist_matrix_cross.groupby(labels_train).mean().idxmin().values
                    )
                    pred_str = calculate_prediction_strength_labels(
                        labels_test, transfered_labels, k
                    )
                    strengths[i, j] = pred_str

                    transfered_labels = (
                        dist_matrix_cross.T.groupby(labels_test).mean().idxmin().values
                    )
                    pred_str = calculate_prediction_strength_labels(
                        labels_train, transfered_labels, k
                    )
                    strengths[i + n_iterations, j] = pred_str
        try:
            kopt = np.array(clusters)[
                strengths.mean(axis=0) > prediction_strength_cutoff
            ][-1]
        except:
            kopt = 1
        full_linkage = fc.average(sp.distance.squareform(DistMatrix))
        labels = hc.fcluster(full_linkage, kopt, criterion="maxclust")
        leader_index, _ = hc.leaders(full_linkage, labels)
        threshold = (
            full_linkage[(leader_index.max() - DistMatrix.shape[0]), 2].min() * 1.01
        )

        HQ_clusters = pd.Series(labels, index=DistMatrix.index)
        output_tsv = (
            pd.DataFrame(HQ_clusters)
            .reset_index()
            .rename(columns={"index": "genome", 0: "subspecies"})
        )

        counts_per_cluster = HQ_clusters.value_counts()

        small_cluters = counts_per_cluster.index[counts_per_cluster < 3]
        HQ_clusters.loc[HQ_clusters.isin(small_cluters)] = pd.NA
        counts_per_cluster = HQ_clusters.value_counts()

        HQ_clusters = HQ_clusters.map(
            {
                oldnr: f"{cluster95} sub{i+1}"
                for i, oldnr in enumerate(counts_per_cluster.index)
            }
        )

        Subspecies = HQ_clusters

        if Subspecies.isnull().any():
            try:
                Subspecies.loc[Subspecies.isnull()] = (
                    DistMatrix.loc[:, Subspecies.isnull()]
                    .groupby(Subspecies)
                    .mean()
                    .idxmin()
                )
            except:
                pass

        output_tsv = output_tsv.rename(columns={"index": "genome"})

        dictionary = defaultdict()
        for subspecies in output_tsv.subspecies.unique():
            df = output_tsv.query("subspecies == @subspecies")
            for genome in df.genome:
                hq_genome = genome
            dictionary[subspecies] = hq_genome
        values = list(dictionary.values())

        output_tsv.subspecies = output_tsv.subspecies.apply(
            lambda x: f"{cluster95:04}001{x:03}"
        )

        output_tsv.to_csv(snakemake.output.subsp_definition, sep="\t", index=False)

        Info = {
            "cluster95": str(cluster95).rjust(4, "0"),
            "cluster965": "001",
            "N_genomes": len(output_tsv),
            "Prediction_strengths": list(strengths.mean(0)),
            "N_clusters": kopt,
            "N_subspecies": len(Subspecies.unique()),
            "Ngenomes_per_cluster": list(
                output_tsv.subspecies.value_counts().sort_index().values
            ),
            "Representatives": list(dictionary.values()),
            "Threshold": threshold,
        }

        Info = pd.DataFrame(pd.Series(Info)).T
        Info.to_csv(snakemake.output.info, sep="\t", index=False)
    else:
        logger.warning(
            f"Cluster {cluster95} doesn't have enough genomes, not clustering it . . ."
        )
        Info = {
            "cluster95": str(cluster95).rjust(4, "0"),
            "cluster965": "001",
            "N_genomes": np.nan,
            "Prediction_strengths": np.nan,
            "N_clusters": np.nan,
            "N_subspecies": 1,
            "Ngenomes_per_cluster": np.nan,
            "Representatives": np.nan,
            "Threshold": np.nan,
            "Min_distance_rep": np.nan,
        }
        Info = pd.DataFrame(pd.Series(Info)).T
        Info.to_csv(snakemake.output.info, sep="\t", index=False)
        genomes = DistMatrix.index.to_list()

        output_tsv = pd.DataFrame(genomes, columns=["genome"])
        output_tsv["subspecies"] = f"{cluster95:04}001001"
        output_tsv["subspecies"] = (
            output_tsv["subspecies"].astype("str").str.rjust(10, "0")
        )
        output_tsv.to_csv(snakemake.output.subsp_definition, sep="\t", index=False)

else:
    logger.warning(
        f"The input file for cluster95 {cluster95} is empty! Touching the output, in order to continue . . ."
    )
    shell(f"touch {snakemake.output.subsp_definition}")
    shell(f"touch {snakemake.output.info}")
