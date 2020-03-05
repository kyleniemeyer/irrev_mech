"""Tests for remove_plog module"""

import sys
import os
import pkg_resources

import pytest
import numpy as np

from ..remove_plog import find_nearest_parameters, interpolate_parameters
from ..remove_plog import convert_plog, remove_plog_reactions
from ..chem_utilities import ReacInfo, Q_

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

class TestFindNearestParameters:
    def test_single_plog(self):
        """Test handling plog reaction with single parameter
        """
        pars = find_nearest_parameters(
            Q_('101325.0 Pa'), [[Q_('101325.0 Pa'), 9.5e39, -9.43, 5636.06]]
            )
        np.testing.assert_allclose(9.5e39, pars[0])
    
    @pytest.mark.parametrize("pressure, pre_exp", [
        (Q_('10132.5 Pa'), 9.2e35), (Q_('101325.0 Pa'), 9.5e39), 
        (Q_('1013250.0 Pa'), 1.5e42), (Q_('10132500.0 Pa'), 1.8e40), 
        (Q_('10132500000.0 Pa'), 4.4e6), (Q_('9000 Pa'), 9.2e35), 
        (Q_('2.e10 Pa'), 4.4e6), (Q_('20000 Pa'), 9.2e35), 
        (Q_('201325 Pa'), 9.5e39), (Q_('913250 Pa'), 1.5e42)
        ])
    def test_closest(self, pressure, pre_exp):
        """Test that converted reaction has appropriate pre-exponential.
        """
        par_list = [
            [Q_('10132.5 Pa'), 9.2e35, -8.65, 3522.54],
            [Q_('101325.0 Pa'), 9.5e39, -9.43, 5636.06],
            [Q_('1013250.0 Pa'), 1.5e42, -9.69, 7598.62],
            [Q_('10132500.0 Pa'), 1.8e40, -8.78, 8454.09],
            [Q_('10132500000.0 Pa'), 4400000.0, 1.45, 1207.73]
            ]
        
        pars = find_nearest_parameters(pressure, par_list)
        np.testing.assert_allclose(pre_exp, pars[0])


class TestInterpolateParameters:
    def test_below(self):
        """Test when pressure is below all given.
        """
        par_list = [
            [Q_('10132.5 Pa'), 9.2e35, -8.65, 3522.54],
            [Q_('101325.0 Pa'), 9.5e39, -9.43, 5636.06],
            [Q_('1013250.0 Pa'), 1.5e42, -9.69, 7598.62],
            [Q_('10132500.0 Pa'), 1.8e40, -8.78, 8454.09],
            [Q_('10132500000.0 Pa'), 4400000.0, 1.45, 1207.73]
            ]
        pars = interpolate_parameters(Q_('1000.0 Pa'), par_list)
        np.testing.assert_allclose(pars, [9.2e35, -8.65, 3522.54])

    def test_above(self):
        """Test when pressure is above all given.
        """
        par_list = [
            [Q_('10132.5 Pa'), 9.2e35, -8.65, 3522.54],
            [Q_('101325.0 Pa'), 9.5e39, -9.43, 5636.06],
            [Q_('1013250.0 Pa'), 1.5e42, -9.69, 7598.62],
            ]
        pars = interpolate_parameters(Q_('10000000.0 Pa'), par_list)
        np.testing.assert_allclose(pars, [1.5e42, -9.69, 7598.62])


class TestConvertPlog:
    def test_not_plog(self):
        """Test for proper error if reaction is not plog.
        """
        rxn = ReacInfo(
            True, ['A'], [1], ['B'], [1], 1.0, 2.0, 3.0
            )
        with pytest.raises(AssertionError):
            convert_plog(rxn, 1.0)


class TestRemovePlogReactions:
    def test_convert_write(self):
        """Tests conversion and write for pressure at a given value
        """
        
        mech = relative_location(os.path.join('mech_plog.dat'))
        therm = relative_location(os.path.join('therm.dat'))

        with open(relative_location(os.path.join('mech_plog_converted.dat')), 'r') as f:
            blessed = f.read()
        
        with TemporaryDirectory() as temp_dir:
            output_file = os.path.join(temp_dir, 'mech.out')
            remove_plog_reactions(mech, therm, '1.0 atm', output_file)
            with open(output_file, 'r') as f:
                output = f.read()
            
        assert blessed == output
    