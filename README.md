2018_ecoli_pathogenicity
===

GWAS on E. coli intrinsic virulence phenotype

Dependencies
------------

- Input genomes are available through [FigShare](https://figshare.com/articles/Escherichia_coli_pathogenicity_GWAS_input_genome_sequences/8866259), and the archive should be uncompressed into `data/genomes`
- All necessary software can be installed through `conda` or `mamba`: `conda create -n 2018_ecoli_pathogenicity pyseer prokka harvesttools gubbins roary unitig-counter bwa bedtools mash blast ncbi-genome-download ete3 jupyterlab snakemake staramr` followed by `conda activate 2018_ecoli_pathogenicity`
- A copy of the [UniRef50](ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/) database should be made into a blast database and available as `db/uniref50`

Usage
-----

To run the full analysis you can use the provided [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline, which will produce all the plots showed in the manuscript and all intermediate files:

```
snakemake -p --cores CPU all 
```

where `CPU` is the number of cores available. Please bare in mind that the whole pipeline may require more than 16G or RAM and take more than a day to complete.

Reference
---------

[Major role of the high-pathogenicity island (HPI) in the intrinsic extra-intestinal virulence of *Escherichia coli* revealed by a genome-wide association study](https://www.biorxiv.org/content/10.1101/712034v1), biorxiv, doi: 10.1101/712034

More reproducibility
---------

The `conda.txt` and `pip.txt` files contain the exact software version used to run this analysis.

Copyright
---------

Copyright 2018 EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

