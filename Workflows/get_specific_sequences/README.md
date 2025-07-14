This is a SnakeMake workflow designed for extracting subspecies-specific genes based on subspecies-specific hashes.

In short, the workflow:
1)	Splits the input to cluster95 (species) level;
2)	Prefetches the required sketches of genomes. The idea is to find the minimal set of input genomes of a subspecies that contain all of the subspecies-specific hashes;
3)	Extracts kmer sequences from those genomes;
4)	Looks at the prodigal output and greps all genes where those kmers are found;
5)	Creates MMSeq2 database from those genes and clusters them on 90% ID. The idea is not to functionally annotate all genes separately, but just cluster representatives;
6)	Annotates those cluster representatives using DRAM;
7)	It uses that annotation to group genes based on KEGG function. Then, it separates genes that are fully specific to a subspecies and genes that are shared;
8)	Does MSA using muscle on genes that are shared
9)	Parses the MSA output and filters so that it keeps functions that are not too diverse outside of the subspecies context and extracts shared variants;
10)	Runs PROSTATA to predict a variant effect on the shared variants.  
