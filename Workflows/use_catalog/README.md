This is a SnakeMake workflow designed for using subspecies-specific hashes to quantify subspecies abundances from metagenomic sequencing data.

For that, you need:
1)	HuMSub catalog – indexed sourmash signature file
2)	sample_table – a text file that can be generated using metagenome-atlas init command. The workflow requires that you first do QC on sequencing data using atlas run qc.

In short, the workflow:
1)	Optionally, you can choose if you want to trim sequencing reads, as sometimes suggested by sourmash. See https://sourmash.readthedocs.io/en/latest/using-sourmash-a-guide.html#how-should-i-prepare-my-data for details;
2)	Sketches the input sequencing reads using sourmash;
3)	Finally, it runs the “gather” command from sourmash to quantify specific hashes.