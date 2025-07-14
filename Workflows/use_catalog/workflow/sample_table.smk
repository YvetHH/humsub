import pandas as pd


ADDITIONAL_SAMPLEFILE_HEADERS = []


def load_sample_table(sample_table="samples.tsv"):
    sampleTable = pd.read_csv(sample_table, index_col=0, sep="\t")
    validate_sample_table(sampleTable)
    return sampleTable


def get_sample_table():
    """Delayed evaluation to avoid early FileNotFoundError"""
    sample_table_path = config.get("sample_table", "samples.tsv")
    if not os.path.exists(sample_table_path):
        raise FileNotFoundError(f"Sample table {sample_table_path} not found.")
    return load_sample_table(sample_table_path)


def validate_sample_table(sampleTable):
    Expected_Headers = ADDITIONAL_SAMPLEFILE_HEADERS
    for h in Expected_Headers:
        if not (h in sampleTable.columns):
            raise Exception(f"Expect '{h}' to be found in samples.tsv")

    if not sampleTable.index.is_unique:
        duplicated_samples = ", ".join(
            sampleTable.index[sampleTable.index.duplicated()]
        )
        raise Exception(
            f"Expect Samples to be unique. Found {duplicated_samples} more than once"
        )

    if not (sampleTable.columns.str.startswith("Reads_QC_").sum() >= 1):
        raise IOError(
            "QC reads need to be in the sample table. "
            "Didn't find any columns with Reads_QC_<fraction>"
        )


class FileNotInSampleTableException(Exception):
    def __init__(self, message):
        super(FileNotInSampleTableException, self).__init__(message)


def get_files_from_sampleTable(sample, Headers):
    sampleTable = get_sample_table()
    if sample not in sampleTable.index:
        raise FileNotInSampleTableException(
            f"Sample name {sample} is not in sampleTable"
        )

    if isinstance(Headers, str):
        Headers = [Headers]

    NheadersFound = sampleTable.columns.isin(Headers).sum()

    if NheadersFound == 0:
        raise FileNotInSampleTableException(
            f"None of the Files are in sampleTable. They should be added later."
            f"\nSample: {sample}\nFiles: {Headers}"
        )
    elif NheadersFound < len(Headers):
        raise IOError(
            f"Not all of the Headers are in sampleTable (found {NheadersFound})."
            f"\nSample: {sample}\nFiles: {Headers}"
        )

    files = sampleTable.loc[sample, Headers]

    if files.isnull().all():
        raise FileNotInSampleTableException(
            f"The following files were not available for sample {sample} in the SampleTable"
        )
    elif files.isnull().any():
        raise IOError(f"Not all files present for sample {sample}. Some are missing.")

    return list(files)


def get_quality_controlled_reads(wildcards):
    sampleTable = get_sample_table()
    if sampleTable.columns.str.contains("R2").any():
        MULTIFILE_FRACTIONS = ["R1", "R2"]
    else:
        MULTIFILE_FRACTIONS = ["se"]

    QC_Headers = ["Reads_QC_" + f for f in MULTIFILE_FRACTIONS]
    sample_dir_path = "/".join(config["sample_table"].split("/")[:-1]) + "/"
    return [
        sample_dir_path + s
        for s in get_files_from_sampleTable(wildcards.sample, QC_Headers)
    ]
