"""Main module"""
import sys

from .mech_util import mech_util


def main(args=None):
    if args is None:
        args = sys.argv[1:]
        mech_util(args)


if __name__ == '__main__':
    sys.exit(main())
