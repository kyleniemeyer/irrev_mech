"""Driver function for package
"""

import sys
import os
from argparse import ArgumentParser

from ._version import __version__
from .irrev_mech import convert_mech_irrev
from .remove_plog import remove_plog_reactions


def mech_util(argv):
    """
    """
    parser = ArgumentParser(description='mech_util: Kinetic model utility.')

    parser.add_argument(
        '-m', '--model',
        help='input model filename for conversion (e.g., "mech.dat").',
        type=str,
        required=True
        )
    parser.add_argument(
        '--thermo',
        help='thermodynamic data filename (only necessary for Chemkin files).',
        type=str,
        default=None
        )
    
    parser.add_argument(
        '--remove_irrev',
        help='Flags conversion of all reversible reactions to irreversible only',
        action='store_true',
        default=False,
        )

    parser.add_argument(
        '--remove_plog',
        help='Flags removal of PLOG reactions',
        action='store_true',
        default=False,
        )
    
    parser.add_argument(
        '--method',
        help='Method for removing PLOG reactions (interpolate, or use nearest)',
        choices=['interpolate', 'nearest'],
        type=str,
        default='interpolate',
        )

    parser.add_argument(
        '--pressure',
        help='Specifies pressure for removal of PLOG reactions',
        type=str,
        default='1.0 atm',
        )

    parser.add_argument(
        '--temp_range',
        help='Specifies lower and upper temperature limits for fitting irreversible reactions',
        nargs=2,
        type=float,
        default=[300., 5000.],
        )

    parser.add_argument(
        '--num_threads',
        help=('Number of CPU cores to use for working in parallel. '
              'If no number, then use available number of cores minus 1.'
              ),
        nargs='?',
        const=0,
        default=1,
        type=int
        )
    
    parser.add_argument(
        '--output',
        help='output file for new model',
        type=str,
        default='mech.out'
        )

    parser.add_argument(
        '-V', '--version',
        action='store_true',
        help='Show the version of mech_util and quit'
        )

    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if args.version:
        print('mech_util {version} from {path} ()'.format(
            version=__version__,
            path=os.path.abspath(os.path.dirname(__file__))
            ))
        sys.exit(0)
    
    if not args.remove_irrev and not args.remove_plog:
        parser.error('Need to specify either remove_irrev or remove_plog')

    if args.remove_irrev:
        convert_mech_irrev(
            args.model, args.thermo, args.temp_range, 
            args.output, args.num_threads
            )
    
    if args.remove_plog:
        remove_plog_reactions(
            args.model, args.thermo, args.pressure, args.output, args.method
            )
