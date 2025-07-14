# HuMSub Catalogue  
Subspecies of the human gut microbiota carry implicit information for in-depth microbiome research

This repository provides all data and code needed to reproduce the results from our publication.

## Repository Structure

The code is organized into:

### Workflows (`/Workflows`)
Snakemake pipelines used in the study:

1. `Quality-filtering` – MAG filtering using GUNC and BUSCO  
2. `sourmash-compare` – Subspecies delineation  
3. `sourmash-hash-selection` – Generation of subspecies-specific sketch databases  
4. `use_catalog` – Subspecies quantification from metagenomes  
5. `assign-subspecies` – Assigning new genomes to HuMSub clusters  
6. `get_specific_sequences` – Identifying subspecies-specific genes

### Scripts (`/Scripts`)
Jupyter notebooks for:
- Statistical and machine learning analyses
- Figure generation

Notebooks are grouped by topic, matching the publication’s structure. Each includes all necessary data for re-execution.

## HuMSub Catalogue Data

The catalogue is available in two prebuilt sourmash SBT index formats:

| File                        | Use Case                                |
|-----------------------------|------------------------------------------|
| `HuMSub_51_1000.sbt.zip`    | General subspecies quantification (k=51) |
| `HuMSub_21_1000.sbt.zip`    | Mastiff database queries (k=21)          |

Download from Zenodo:  
https://zenodo.org/records/15862096

## Simulated Metagenomic Data

For benchmarking and testing:

| File                      | Description                                                  |
|---------------------------|--------------------------------------------------------------|
| `humgut_samples.tar.gz`   | Simulated paired-end reads from HumGut genomes               |
| `new_samples.tar.gz`      | Simulated paired-end reads from genomes outside of HumGut    |

Corresponding taxonomic distributions are available in `Scripts/benchmark/`.

## Disclaimer

The HuMSub catalogue includes genomes from some non-gut-associated phyla (e.g., Elusimicrobiota, Eremiobacteriota, Patescibacteria) retained from the original HumGut reference.  
Although these were not detected in CRC datasets, they were preserved for completeness based on genome quality scores. Use with caution in downstream interpretation.

## Getting Started

1. Install Snakemake (https://snakemake.readthedocs.io) version 7
2. Clone this repository
3. Modify the appropriate `config.yaml` for each workflow
4. Run with:

```bash
snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --cores 4

# External Resources

This directory lists essential external files required to reproduce the results of the HuMSub study. These resources are hosted externally due to their size and licensing constraints.

## 1. HumGut Genome Metadata

- **File**: `All_genomes.tsv`  
- **Description**: Metadata and download links for all genomes in the HumGut catalog
- **Source**:  
  https://arken.nmbu.no/~larssn/humgut/All_genomes.tsv

## 2. HuMSub Catalogue and Benchmark Data

- **Zenodo Record**:  
  https://zenodo.org/records/15862096

This includes:

- `HuMSub_51_1000.sbt.zip` – k=51, scaled=1000 (subspecies quantification)
- `HuMSub_21_1000.sbt.zip` – k=21, scaled=1000 (Mastiff queries)
- `humgut_samples.tar.gz` – simulated reads from HumGut genomes
- `new_samples.tar.gz` – simulated reads from genomes outside HumGut

After download, place files into the appropriate subdirectories (e.g. `resources/`, `test_data/`).

## Note

Some files may be automatically downloaded by Snakemake rules if not found in the expected locations. Refer to the main `README.md` for pipeline instructions.