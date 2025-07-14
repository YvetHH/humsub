This is a SnakeMake workflow designed for using subspecies-specific hashes to assign new genomes to existing subspecies.
To run it, you need:
1)	Input genomes – either all in one directory or a text file with specifying paths;
2)	‘specific_hashomes’ – HuMSub catalog

In short, the workflow:
1)	Creates batches of input genomes for parallel processing;
2)	Sketches the input genomes. Ensure that the sketching parameters (kmer length and scaling factor match the HuMSub catalog you are using);
3)	Compares the genome sketch with the subspecies-specific hashes. 

The output is a tsv table with indicated overlap between the input genome and each of the subspecies in the catalog.
