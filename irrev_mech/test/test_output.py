"""Test the output of the irrev_mech script"""
import os
from irrev_mech.irrev_mech import convert_mech_irrev


def test_convert_mech_irrev():
    pth = os.path.dirname(os.path.realpath(__file__))
    args = {'mech_name': os.path.join(pth, 'mech.dat'), 'therm_name': os.path.join(pth, 'therm.dat'), 'temp_range': [300.0, 5000.0],
            'output_file': os.path.join(pth, 'mech_output.txt'), 'n_procs': 1}
    convert_mech_irrev(**args)
    blessed = open(os.path.join(pth, 'mech_blessed.txt'), 'r').read()
    output = open(os.path.join(pth, 'mech_output.txt'), 'r').read()
    assert output == blessed
    if os.path.exists(os.path.join(pth, 'mech_output.txt')):
        os.remove(os.path.join(pth, 'mech_output.txt'))
