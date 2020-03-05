mech_util
=========

[![DOI](https://zenodo.org/badge/5316746.svg)](https://zenodo.org/badge/latestdoi/5316746)

This utility manipulates a Chemkin-format chemical kinetic mechanism in one of two ways:
- converts reversible reactions to irreversible reactions
- removes PLOG (pressure-log) reactions

`irrev_mech` requires [Python] 3.x and the [SciPy] stack, but no additional libraries (i.e., Chemkin, Cantera). It has only been tested on [Python] 3.5 and 3.6.

Usage
-------

`mech_util` can be used locally via

    $ python -m mech_util [options]

or installed as a package using `pip install --user mech_util` or via `python setup.py install`, and called using

    $ mech_util [options]

Use the option `-h` or `--help` to see the full usage instructions. To generate an irreversible mechanism, from the command line, use `python -m mech_util -c mech_name -t thermname` where `mech_name` and `therm_name` are the names of the chemical kinetics reaction mechanism file and thermodynamic database, e.g.:

    $ python -m mech_util --model mech.dat --thermo therm.dat --remove_irrev

The new model file has the name `mech_irrev.txt`.

The default temperature range used for parameter fitting is 300 K to 5000 K. This can be changed by specifying the `-r` or `--range` command line option, e.g.:

    $ python -m mech_util --model mech.dat --thermo therm.dat --remove_irrev --temp_range 1000 3000

License
-------

`mech_util` is released under the MIT license, see LICENSE for details.

Citation
--------
If you use this software as part of a scholarly publication, please cite the software directly using the DOI in the badge at the top of this file (and found here: https://zenodo.org/badge/latestdoi/5316746).


Further Reading
---------------

`mech_util` converts reversible reactions into irreversible reactions by fitting reverse Arrhenius coefficients using a nonlinear least-squares minimization.
The three reverse Arrhenius coefficients are determined using a nonlinear least-squares minimization, using the [SciPy] function `scipy.optimize.leastsq`. The fit is performed using a large number of calculated reverse rate coefficients, densely sampled over the specified temperature range. The initial guess for this minimization comes from an analytical fit to three temperatures (the high and low values of the range, plus a midpoint). See the appendix of Niemeyer and Sung's paper for more details:

* KE Niemeyer and CJ Sung. "Accelerating moderately stiff chemical kinetics in reactive-flow simulations using GPUs." *J. Comput. Phys.*, 256:854-871, 2014. doi:[10.1016/j.jcp.2013.09.025](http://dx.doi.org/10.1016/j.jcp.2013.09.025)

The use of the least-squares minimization to obtain better fits was suggested by Taylor et al.:

 * BD Taylor, DA Schwer, and A Corrigan. "Implementation of Thermochemistry and Chemical Kinetics in a GPU-based CFD Code." 53rd AIAA Aerospace Sciences Meeting, Kissimmee, FL, January 2015. doi:[10.2514/6.2015-0842](http://dx.doi.org/10.2514/6.2015-0842)

`mech_util` "de-plogs" models by using a given pressure to either choose the Arrhenius
parameters at the nearest pressure, or by interpolating parameter values.

Author
------

Created by [Kyle Niemeyer](https://niemeyer-research-group.github.io). Email address: [kyle.niemeyer@oregonstate.edu](mailto:kyle.niemeyer@oregonstate.edu)


[Python]: http://python.org/
[SciPy]: http://scipy.org/
