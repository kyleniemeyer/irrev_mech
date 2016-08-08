#! /usr/bin/env python
from __future__ import absolute_import

# Standard libraries
import sys
from argparse import ArgumentParser

# Local import
from .irrev_mech import convert_mech_irrev

def get_parser():
    # command line arguments
    parser = ArgumentParser(description = 'Generates chemical kinetic '
                                          'reaction mechanism with only '
                                          'irreversible reactions.',
                            prog='irrev_mech',
                            )
    parser.add_argument('-c', '--chem',
                        type = str,
                        required = True,
                        help = 'Input mechanism filename (e.g., mech.dat).'
                        )
    parser.add_argument('-t', '--thermo',
                        type = str,
                        default = None,
                        help = 'Thermodynamic database filename (e.g., '
                               'therm.dat), or nothing if in mechanism.'
                        )
    parser.add_argument('-r', '--range',
                        type = float, nargs=2,
                        default = [300.0, 5000.0],
                        help = 'Temperature range for fit in Kelvin '
                               '(e.g., 300 5000).'
                        )
    parser.add_argument('-o', '--output',
                        type = str,
                        default = 'mech_irrev.txt',
                        help = 'Output file name, default "mech_irrev.txt"')
    parser.add_argument('-n', '--numprocs',
                        type = int,
                        default = None,
                        help = 'Number of processes to use. Default: Number of CPUs on the machine.')

    return parser.parse_args()

def main(args=None):
    if args is None:
        args = get_parser()

    convert_mech_irrev(args.chem, args.thermo, args.range, args.output)

if __name__ == '__main__':
    sys.exit(main())
