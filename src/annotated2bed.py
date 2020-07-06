#!/usr/bin/env python


import argparse


def get_options():
    description = 'Annotated kmers to BED (only GOIs)'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('annotated',
                        help='Annotated kmers file')
    parser.add_argument('genes',
                        help='Genes Of Interest file (tab-separated, second column)')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    genes = {x.rstrip().split()[-1]
             for x in open(options.genes)}

    for l in open(options.annotated):
        strings = l.rstrip().split()[-1]
        for string in strings.split(','):
            loc, _, gene, _ = string.split(';')
            if gene not in genes:
                continue
            # note: removing the .x from the chromosome id
            print('%s\t%s\t%s\t%s' % (loc.split(':')[0].split('.')[0],
                                      loc.split(':')[1].split('-')[0],
                                      loc.split(':')[1].split('-')[1],
                                      gene))
