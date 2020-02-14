#!/usr/bin/env python
# Copyright 2017 Marco Galardini and John Lees

'''Summarise k-mer annotation at the gene level'''


def get_options():
    import argparse

    description = 'Summarise k-mer annotation at the gene level'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('annotations',
                        help='Annotated k-mer file from annotate_hits.py')
    parser.add_argument('pangenome',
                        help='Roary 0/1 file')

    parser.add_argument('--min-size',
                        default=30,
                        type=int,
                        help='Minimum kmer size to be considered specific [Default: 30]')
    parser.add_argument('--lineage',
                        default=None,
                        help='Only consider kmers belonging to a specific lineage')
    parser.add_argument('--nearby',
                        action='store_true',
                        help='Use up/downstream annotation, if not in a gene')

    return parser.parse_args()


def update_summary(summary, gene, pval, af, beta, length, gaf, specific):
    if summary[gene] != {}:
        summary[gene]['count'] += 1
        if specific:
            summary[gene]['specific'] += 1
        summary[gene]['af'] += af
        summary[gene]['beta'] += beta
        if log10p > summary[gene]['maxp']:
            summary[gene]['maxp'] = log10p
        summary[gene]['length'] += length
        summary[gene]['gaf'] = gaf
    else:
        summary[gene]['count'] = 1
        if specific:
            summary[gene]['specific'] = 1
        else:
            summary[gene]['specific'] = 0
        summary[gene]['af'] = af
        summary[gene]['beta'] = beta
        summary[gene]['maxp'] = log10p
        summary[gene]['length'] = length
        summary[gene]['gaf'] = gaf


if __name__ == "__main__":
    options = get_options()

    import sys
    import collections
    import pandas as pd
    from math import log10

    p = pd.read_csv(options.pangenome, sep='\t', index_col=0)
    psum = p.T.sum() / p.shape[1]
    
    summary = collections.defaultdict(dict)
    with open(options.annotations, 'r') as anot_file:
        for line in anot_file:
            anot_fields = line.rstrip().split("\t")
            variant = anot_fields[0]
            if len(variant) >= options.min_size:
                specific = True
            else:
                specific = False
            af = float(anot_fields[1])
            pvalue = float(anot_fields[3])
            beta = abs(float(anot_fields[4]))
            if options.lineage is not None:
                lineage = anot_fields[-2]
                if lineage != options.lineage:
                    continue
            annotations = anot_fields[-1].split(",")

            if pvalue > 0:
                log10p = -log10(pvalue)
                for annotation in annotations:
                    (position, down, inside, up) = annotation.split(";")
                    if inside != "":
                        update_summary(summary,
                                       inside,
                                       log10p,
                                       af,
                                       beta,
                                       len(variant),
                                       psum.loc[inside],
                                       specific)
                    elif options.nearby:
                        if down != "":
                            update_summary(summary,
                                           down,
                                           log10p,
                                           af,
                                           beta,
                                           len(variant),
                                           psum.loc[down],
                                           specific)
                        if up != "":
                            update_summary(summary,
                                           up,
                                           log10p,
                                           af,
                                           beta,
                                           len(variant),
                                           psum.loc[up],
                                           specific)

    # write output
    print("\t".join(["gene", "hits", "specific_hits", 'length', 'OG_af', "maxp",
                     "avg_af", "avg_maf", "avg_beta"]))
    for gene in summary:
        af = summary[gene]['af']/summary[gene]['count']
        if af > 0.5:
            maf = 1 - af
        else:
            maf = af

        print("\t".join([gene,
                         str(summary[gene]['count']),
                         str(summary[gene]['specific']),
                         str(summary[gene]['length']),
                         str(summary[gene]['gaf']),
                         str(summary[gene]['maxp']),
                         str(af),
                         str(maf),
                         str(summary[gene]['beta']/summary[gene]['count'])]))

