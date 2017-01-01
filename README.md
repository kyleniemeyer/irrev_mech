irrev_mech
=======

[![DOI](https://zenodo.org/badge/5316746.svg)](https://zenodo.org/badge/latestdoi/5316746)

This utility converts a Chemkin-format reaction mechanism with reversible reactions to one with only irreversible reactions.

It does this by fitting reverse Arrhenius coefficients using a nonlinear least-squares minimization.

`irrev_mech` requires [Python] (2.x or 3.x) and the [SciPy] stack, but no additional libraries (i.e., Chemkin, Cantera). It has only been tested on [Python] 2.7, but probably works on earlier versions.

Usage
-------

`irrev_mech` can be used locally via

    $ python -m irrev_mech [options]

or installed as a package using `pip install --user irrev_mech` or via `python setup.py install`, and called using

    $ irrev_mech [options]

Use the option `-h` or `--help` to see the full usage instructions. To generate an irreversible mechanism, from the command line, use `python -m irrev_mech -c mech_name -t thermname` where `mech_name` and `therm_name` are the names of the chemical kinetics reaction mechanism file and thermodynamic database, e.g.:

    $ python irrev_mech.py -c mech.dat -t therm.dat

You can also run `irrev_mech` without a thermodynamic database if the information is held in the chemistry mechanism file (after the species are declared), e.g.:

    $ python -m irrev_mech -c mech.dat

The new model file has the name `mech_irrev.txt`.

The default temperature range used for parameter fitting is 300 K to 5000 K. This can be changed by specifying the `-r` or `--range` command line option, e.g.:

    $ python -m irrev_mech -c mech.dat -t therm.dat -r 1000 3000

License
-------

`irrev_mech` is released under the MIT license, see LICENSE for details.

Citation
--------
If you use this software as part of a scholarly publication, please cite the software directly using the DOI in the badge at the top of this file (and found here: https://zenodo.org/badge/latestdoi/5316746).


Further Reading
---------------

The three reverse Arrhenius coefficients are determined using a nonlinear least-squares minimization, using the [SciPy] function `scipy.optimize.leastsq`. The fit is performed using a large number of calculated reverse rate coefficients, densely sampled over the specified temperature range. The initial guess for this minimization comes from an analytical fit to three temperatures (the high and low values of the range, plus a midpoint). See the appendix of Niemeyer and Sung's paper for more details:

* KE Niemeyer and CJ Sung. "Accelerating moderately stiff chemical kinetics in reactive-flow simulations using GPUs." *J. Comput. Phys.*, 256:854-871, 2014. doi:[10.1016/j.jcp.2013.09.025](http://dx.doi.org/10.1016/j.jcp.2013.09.025)

The use of the least-squares minimization to obtain better fits was suggested by Taylor et al.:

 * BD Taylor, DA Schwer, and A Corrigan. "Implementation of Thermochemistry and Chemical Kinetics in a GPU-based CFD Code." 53rd AIAA Aerospace Sciences Meeting, Kissimmee, FL, January 2015. doi:[10.2514/6.2015-0842](http://dx.doi.org/10.2514/6.2015-0842)

Author
------

Created by [Kyle Niemeyer](http://kyleniemeyer.com). Email address: [kyle.niemeyer@gmail.com](mailto:kyle.niemeyer@gmail.com)


[Python]: http://python.org/
[SciPy]: http://scipy.org/
