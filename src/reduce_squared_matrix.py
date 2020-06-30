#!/usr/bin/env python


import sys
import argparse
import pandas as pd


def get_options():
    description = 'Reduce squared matrix using another dataframe'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('df',
                        help='Samples dataframe')
    parser.add_argument('matrix',
                        help='Matrix to reduce')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    m = pd.read_csv(options.df, sep='\t', index_col=0)
    s = pd.read_csv(options.matrix, sep='\t', index_col=0)
    idx = s.index.intersection(m.index)
    s.loc[idx, idx].to_csv(sys.stdout, sep='\t')
