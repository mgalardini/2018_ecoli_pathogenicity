#!/usr/bin/env python


import os
import sys
import argparse
import numpy as np
import pandas as pd


def get_options():
    description = 'Collate the annotated SNPs tables'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('snippy',
                        help='Snippy directory')

    parser.add_argument('--fname',
                        default='annotated_snps.tsv',
                        help='Filename for each individual table (default: %(default)s)')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    m = None
    for strain in os.listdir(options.snippy):
        try:
            n = pd.read_csv(os.path.join(options.snippy, strain, options.fname), sep='\t')
        except:
            continue
        n[strain] = 1
        if m is None:
            m = n.set_index(['acc', 'pos', 'ref', 'alt',
                             'locus', 'type', 'score',
                             'median_ic', 'n_aa', 'n_seq'])
        else:
            m = m.join(n.set_index(['acc', 'pos', 'ref', 'alt',
                             'locus', 'type', 'score',
                             'median_ic', 'n_aa', 'n_seq']), how='outer')

    m[np.isnan(m)] = 0.

    m = m.reset_index()
    m['index'] = ['_'.join([str(x) for x in y])
                  for y in m[['acc', 'pos', 'ref', 'alt',
                              'locus', 'type', 'score',
                              'median_ic', 'n_aa', 'n_seq']].values]

    m = m.drop(columns=['acc', 'pos', 'ref', 'alt',
                        'locus', 'type', 'score',
                        'median_ic', 'n_aa', 'n_seq'])

    m = m.set_index('index')

    m.astype(int).to_csv(sys.stdout, sep='\t')
