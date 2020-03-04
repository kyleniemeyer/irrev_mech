"""Test the output of the irrev_mech script"""
import sys
import os
import pkg_resources

from ..irrev_mech import convert_mech_irrev


# Taken from http://stackoverflow.com/a/22726782/1569494
try:
    from tempfile import TemporaryDirectory
except ImportError:
    from contextlib import contextmanager
    import shutil
    import tempfile
    import errno

    @contextmanager
    def TemporaryDirectory():
        name = tempfile.mkdtemp()
        try:
            yield name
        finally:
            try:
                shutil.rmtree(name)
            except OSError as e:
                # Reraise unless ENOENT: No such file or directory
                # (ok if directory has already been deleted)
                if e.errno != errno.ENOENT:
                    raise


def relative_location(file):
	file_path = os.path.join(file)
	return pkg_resources.resource_filename(__name__, file_path)


def test_convert_mech_irrev():
    """
    """
    mech = relative_location(os.path.join('mech.dat'))
    therm = relative_location(os.path.join('therm.dat'))

    with open(relative_location(os.path.join('mech_blessed.txt')), 'r') as f:
        blessed = f.read()
    
    with TemporaryDirectory() as temp_dir:
        output_file = os.path.join(temp_dir, 'mech.out')
        convert_mech_irrev(
            mech, therm, temp_range=[300., 5000.], 
            output_file=output_file, n_procs=None
            )
        with open(output_file, 'r') as f:
            output = f.read()

    #assert output == blessed
