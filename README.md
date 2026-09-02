# Pipeline To Get Orthologs From Genome Assemblies #

## Reference-Based Genome Assembly Using bwa-mem2

Assembles genomes from Illumina short-read sequencing data (e.g., low-coverage ~10X) using a reference-based mapping approach with `bwa-mem2`. Reads are first cleaned and trimmed with ```AdapterRemoval```. Assembly completeness and quality are assessed using `BUSCO`.

## Masking Repeats in Assembled Contigs

Identifies repeat families de novo with `RepeatModeler`, then soft-masks repeats in the assembly using `RepeatMasker` and `WindowMasker`.

## Pairwise Genome Alignment Chains Using make_lastz_chains

Generates pairwise genome alignment chains between the genome of interest and a chosen reference genome (default: human hg38, Ensembl) using the `LASTZ` pairwise DNA sequence aligner, orchestrated via `make_lastz_chains` and run via `Nextflow`.

## Get Orthologs from Genome Alignments Using TOGA2

Infers orthologous genes from genome alignments using `TOGA2`, run via `Nextflow`. Outputs include a multifasta of orthologous sequences, a list of missing orthologs, query genome annotation in BED/GTF format, an orthology classification table (intact/partially intact/lost/uncertain loss per transcript), gene loss and duplication calls with underlying inactivating mutations, protein and codon alignments, and a BED file of processed pseudogenes.
