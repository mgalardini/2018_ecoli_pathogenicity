#!/bin/bash

mkdir -p out/plots/
src/plot_tree out/gubbins/tree.nwk data/phenotypes/phenotypes.tsv out/associations/kmer_counts_lmm.txt out/associations/kmer_gene_counts_lmm.txt out/associations/rtab_gene_counts_lmm.txt data/phylogroups.tsv out/plots/ --only-phenotype --width 5 --dpi 300 --extension pdf --prefix 1_
src/plot_tree out/gubbins/tree.nwk data/phenotypes/phenotypes.tsv out/associations/kmer_counts_lmm.txt out/associations/kmer_gene_counts_lmm.txt out/associations/rtab_gene_counts_lmm.txt data/phylogroups.tsv out/plots/ --width 5 --dpi 300 --extension pdf --prefix 2_
src/plot_tree out/gubbins/tree.nwk data/phenotypes/phenotypes.tsv out/associations/kmer_counts_lmm.txt out/associations/kmer_gene_counts_lmm.txt out/associations/rtab_gene_counts_lmm.txt data/phylogroups.tsv out/plots/ --only-phenotype --width 10 --dpi 300 --extension pdf --prefix 3_ --labels
src/plot_tree out/gubbins/tree.nwk data/phenotypes/phenotypes.tsv out/associations/kmer_counts_lmm.txt out/associations/kmer_gene_counts_lmm.txt out/associations/rtab_gene_counts_lmm.txt data/phylogroups.tsv out/plots/ --width 10 --dpi 300 --extension pdf --prefix 4_ --labels
src/plot_tree out/gubbins/tree.nwk data/phenotypes/phenotypes.tsv out/associations/kmer_counts_lmm.txt out/associations/kmer_gene_counts_lmm.txt out/associations/rtab_gene_counts_lmm.txt data/phylogroups.tsv out/plots/ --only-phenotype-phylogroups --width 5 --dpi 300 --extension pdf --prefix 6_
src/plot_virulence out/gubbins/tree.nwk data/phenotypes/phenotypes.tsv data/phylogroups.tsv out/roary/gene_presence_absence.Rtab out/virulence_genes.tsv out/plots/ --width 5 --dpi 300 --extension pdf --prefix 7_
