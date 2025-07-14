# HuMSub catalogue - Human gut microbiota subspecies carry implicit information for in-depth microbiome research

This repository provides access to all data and code required to replicate findings from the publication.

Since the HumGut catalog is too big to upload here, you can access its metadata at: https://arken.nmbu.no/~larssn/humgut/All_genomes.tsv.
There, you can use the links provided to download all of the genomes from their home repositories.

DISCLAIMER: The HuMSub catalogue contains genomes beloning to phyla not usually associated with the human gut (Elusimicrobiota; Eremiobacteriota; Patescibacteria), which were kept from the initial HumGut reference. While their true association with the human gut remains uncertain, we retained them to ensure completeness, relying on genome quality scores. These phyla were not detected in the follow up analyses of CRC datasets, and should be interpreted with caution.

We have split the code into "Workflows" and "Scripts". In the former, you can find all relevant SnakeMake workflows used for the publication, including:

1) Quality filtering workflow that used GUNC and BUSCO for enhanced MAGs filtering ("Quality-filtering")

2) Subspecies delineation workflow ("sourmash-compare")

3) Subspecies-specific panhashome generation workflow ("sourmash-hash-selection")

4) Subspecies quantification workflow ("use_catalog")

5) Assign new genomes to existing subspecies workflow ("assign-subspecies")

6) Subspecies-specific genes workflow ("get_specific_sequences")

 To run them, please install SnakeMake, create/modify required input files and configure the corresponding "config" files.

 

In the "Scripts" directory, you will find Jupyter Notebooks used for doing statistical/ML analysis and creating figures. They are sorted into different directories, following the structure of the publication. All data to run them is available in the corresponding subdirectories.

 

The catalog itself is available in two versions:

1) HuMSub_51_1000.sbt.zip: Subspecies-specific hashes created based on kmer length of 51 and scaling factor of 1000. It is used for general subspecies quantification;

2) HuMSub_21_1000.sbt.zip:  Subspecies-specific hashes created based on kmer length of 21 and scaling factor of 1000. It is used for quering Mastiff database for studying biogeographical distribution of subspecies.

 

Finally, you can find simulated metagenomic samples:

1) humgut_samples.tar.gz: Contains ten pair-end simulated sequencing data created using MAGs from the HumGut catalog;

2) new_samples.tar.gz: Contains ten pair-end simulated sequencing data created using MAGs outside the HumGut catalog.

The correpoding taxonomical distributions are available within the "Scripts/benchmark".