#!/usr/bin/env python


import sys
import pandas as pd
import argparse


def get_options():
    description = 'From annotated kmers to presence/absence matrix'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('annotation',
                        help='Annotated kmers file')
    parser.add_argument('kmers',
                        help='K-mers file (uncompressed)')
    parser.add_argument('mash',
                        help='Mash distance matrix (to extract the samples list)')

    parser.add_argument('--min-size',
                        default=30,
                        type=int,
                        help='Minimum kmer size to be considered specific [Default: 30]')
    parser.add_argument('--nearby',
                        action='store_true',
                        help='Use up/downstream annotation, if not in a gene')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    samples = [x for x in next(open(options.mash)).rstrip().split('\t')[1:]]

    d = {}
    with open(options.annotation, 'r') as anot_file:
        for line in anot_file:
            anot_fields = line.rstrip().split("\t")
            variant = anot_fields[0]
            if len(variant) < options.min_size:
                continue
            d[variant] = set()
            annotations = anot_fields[-1].split(",")
            for annotation in annotations:
                (position, down, inside, up) = annotation.split(";")
                if inside != "":
                    d[variant].add(inside)
                elif options.nearby:
                    if down != "":
                        d[variant].add(down)
                    if up != "":
                        d[variant].add(up)

    res = {}
    for l in open(options.kmers):
        variant, s = l.rstrip().split('|')
        variant = variant.rstrip()
        s = {x.split(':')[0] for x in s.lstrip().split(' ')}
        if variant not in d:
            continue
        for gene in d[variant]:
            res[gene] = res.get(gene, {})
            for x in samples:
                if x in s:
                    res[gene][x] = 1
                elif x not in res[gene]:
                    res[gene][x] = 0

    pd.DataFrame.from_dict(res, orient='index').to_csv(sys.stdout, sep='\t')
