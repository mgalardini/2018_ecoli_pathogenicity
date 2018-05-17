#!/usr/bin/env python
'''
Run and save a jupyter notebook by providing arguments from the command line

Uses the nbparameterise package, and requires the first cell to have variables
definitions.

Command line parsing based on this Stack Overflow answer:
https://stackoverflow.com/a/42355279/1237531
'''

from nbparameterise import extract_parameters, parameter_values, replace_definitions
import nbformat

def get_options():
    import argparse

    args_dict = {}
    class StoreDictKeyPair(argparse.Action):
         def __init__(self, option_strings, dest, nargs=None, **kwargs):
             self._nargs = nargs
             super(StoreDictKeyPair, self).__init__(option_strings, dest, nargs=nargs, **kwargs)
         def __call__(self, parser, namespace, values, option_string=None):
             for kv in values:
                 k,v = kv.split("=")
                 args_dict[k] = v
             setattr(namespace, self.dest, args_dict)

    description = "Run a Jupyter notebook by passing arguments to it"
    parser = argparse.ArgumentParser(description=description,
                                     prog='run_notebook.py')
    parser.add_argument('template', action='store',
                        help='Template notebook')
    parser.add_argument('output', action='store',
                        help='Output notebook')
    
    parser.add_argument('-k', '--key-value',
                        dest='args_dict',
                        action=StoreDictKeyPair,
                        nargs='+',
                        metavar='KEY=VAL')

    return parser.parse_args()


if __name__ == '__main__':
    options = get_options()

    if options.args_dict is None:
        options.args_dict = {}

    with open(options.template) as template:
        nb = nbformat.read(template, as_version=4)

    orig_parameters = extract_parameters(nb)
    params = parameter_values(orig_parameters, **options.args_dict)

    print('Old parameters: %s' % orig_parameters)
    print('New parameters: %s' % params)

    new_nb = replace_definitions(nb, params, execute=False)

    with open(options.output, 'w') as output:
        nbformat.write(new_nb, output)
