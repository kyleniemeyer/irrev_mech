"""Tests for remove_plog module"""

import sys
import os
import pkg_resources

import pytest
import numpy as np

from ..remove_plog import convert_plog, remove_plog_reactions
from ..chem_utilities import ReacInfo

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

class TestConvertPlog:
    def test_not_plog(self):
        """Test for proper error if reaction is not plog.
        """
        rxn = ReacInfo(
            True, ['A'], [1], ['B'], [1], 1.0, 2.0, 3.0
            )
        with pytest.raises(AssertionError):
            convert_plog(rxn, 1.0)
    
    def test_single_plog(self):
        """Test handling plog reaction with single parameter
        """
        rxn = ReacInfo(
            True, ['A'], [1], ['B'], [1], 1.0, 2.0, 3.0
            )
        rxn.plog = True 
        rxn.plog_par = [[101325.0, 9.5e39, -9.43, 5636.06]]
        
        new_rxn = convert_plog(rxn, 101325.0)
        np.testing.assert_allclose(9.5e39, new_rxn.A)
    
    @pytest.mark.parametrize("pressure, pre_exp", [
        (10132.5, 9.2e35), (101325.0, 9.5e39), (1013250.0, 1.5e42),
        (10132500.0, 1.8e40), (10132500000.0, 4.4e6),
        (9000, 9.2e35), (2.e10, 4.4e6),
        (20000, 9.2e35), (201325, 9.5e39), (913250, 1.5e42)
        ])
    def test_closest(self, pressure, pre_exp):
        """Test that converted reaction has appropriate pre-exponential.
        """
        rxn = ReacInfo(
            True, ['A'], [1], ['B'], [1], 1.0, 2.0, 3.0
            )
        rxn.plog = True 
        rxn.plog_par = [
            [10132.5, 9.2e35, -8.65, 3522.54],
            [101325.0, 9.5e39, -9.43, 5636.06],
            [1013250.0, 1.5e42, -9.69, 7598.62],
            [10132500.0, 1.8e40, -8.78, 8454.09],
            [10132500000.0, 4400000.0, 1.45, 1207.73]
            ]
        
        new_rxn = convert_plog(rxn, pressure)
        np.testing.assert_allclose(pre_exp, new_rxn.A)

class TestRemovePlogReactions:
    def test_convert_write(self):
        
        mech = relative_location(os.path.join('mech_plog.dat'))
        therm = relative_location(os.path.join('therm.dat'))

        with open(relative_location(os.path.join('mech_plog_converted.dat')), 'r') as f:
            blessed = f.read()
        
        with TemporaryDirectory() as temp_dir:
            output_file = os.path.join(temp_dir, 'mech.out')
            remove_plog_reactions(mech, therm, '2.0 atm', output_file)
            with open(output_file, 'r') as f:
                output = f.read()
            
        #assert blessed == output

    