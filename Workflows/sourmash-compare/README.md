This is a SnakeMake workflow that is used to cluster individual genomes per species level to subspecies groups.

To run it, it requires an input table with two columns: "path" and "cluster95".
The column "path" denotes a local path to genes called in amino-acid format.
The column "cluster95" denotes a species attribution. Effectively, the scirpt will use this group as base groups for clustering.

In the first rule, it separates paths based on cluster95 information. Then, it sketches that input and then compares using sourmash.

Finally, it takes the output distance matrices and clusters based on them.
The final output are text files for each for the cluster95 (species) with annotated subspecies information.