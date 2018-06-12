#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
 
parser$add_argument("samples",
                    help="RNA samples file")
parser$add_argument("directory",
                    help="kallisto output directory")
parser$add_argument("--reference",
                    default="IAI55", 
                    help="Reference strain for contrasts [default %(default)s]")
parser$add_argument("--cores",
                    type="integer",
                    default=1, 
                    help="Number of cores to use [default %(default)s]")
parser$add_argument("--pvalue",
                    default=0.01,
                    type="double",
                    help="p-value threshold [default %(default)s]")
parser$add_argument("--foldchange",
                    default=0.0,
                    type="double",
                    help="fold-change threshold [default %(default)s]")

args <- parser$parse_args()

suppressPackageStartupMessages(library("tximport"))
suppressPackageStartupMessages(library('DESeq2'))
suppressPackageStartupMessages(library('BiocParallel'))

# prepare samples table
s2c <- read.table(file.path(args$samples),
                  header=TRUE,
                  stringsAsFactors=TRUE)
kal_dirs <- sapply(s2c$sample,
                   function(id) file.path(args$directory,
                                          id,
                                          "abundance.h5"))
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# import kallisto's counts
txi.kallisto <- tximport(s2c$path,
                         type="kallisto",
                         txOut=TRUE,
                         countsCol='est_counts')
counts <- as.matrix(txi.kallisto$counts)
colnames(counts) <- s2c$sample
for (col in colnames(counts)){
    counts[, col] <- sapply(counts[, col], as.integer)
}
rownames(s2c) <- colnames(counts)

# run DESeq2
dds <- DESeqDataSetFromMatrix(counts,
                              s2c,
                              ~strain)
multicoreParam <- MulticoreParam(workers=as.numeric(args$cores))
dds <- DESeq(dds,
             parallel=TRUE,
             BPPARAM=multicoreParam)

# overall results
res <- results(dds,
               alpha=args$pvalue,
               lfcThreshold=args$foldchange,
               altHypothesis='greaterAbs',
               parallel=TRUE,
               BPPARAM=multicoreParam)
write.csv(as.data.frame(res),
          file=file.path(args$directory,
                         'overall.csv'))
# strain's contrasts
strains <- unique(s2c$strain)
strains <- strains[strains != args$reference]
for(i in 1:length(strains))
{
    res <- results(dds,
                   contrast=c('strain',
                              toString(strains[i]),
                              args$reference),
                   alpha=args$pvalue,
                   lfcThreshold=args$foldchange,
                   altHypothesis='greaterAbs',
                   parallel=TRUE,
                   BPPARAM=multicoreParam)
    write.csv(as.data.frame(res),
              file=file.path(args$directory,
                             paste(toString(strains[i]),
                                   '.csv',
                                   sep=''))) 
}
