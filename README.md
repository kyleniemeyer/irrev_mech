irrev_mech
=======

This utility converts a Chemkin-format reaction mechanism with reversible reactions to one with only irreversible reactions.

It does this by fitting reverse Arrhenius coefficients using a nonlinear least-squares minimization.

`irrev_mech` requires [Python] 2.x and the [SciPy] stack, but no additional libraries (i.e., Chemkin, Cantera). It has only been tested on [Python] 2.7, but probably works on earlier versions.

Usage
-------

Use the command `python irrev_mech.py -h` to see the full usage instructions. To generate an irreversible mechanism, from the command line, use `python irrev_mech.py -m mechname -t thermname` where `mechname` and `thermname` are the names of the mechanism file and thermodynamic database, e.g.:

    $ python irrev_mech.py -m mech.dat -t therm.dat

You can also run `irrev_mech` without a thermodynamic database if the information is held in the mechanism file (after the species are declared), e.g.:

    $ python irrev_mech.py -m mech.dat

The new mechanism has the name `mech_irrev.txt`.

The default temperature range used for parameter fitting is 300 K to 5000 K. This can be changed by specifying the `-r` or `--range` command line option, e.g.:

    $ python irrev_mech.py -m mech.dat -t therm.dat -r 1000 3000

Misc
-------

The most up-to-date version of `irrev_mech` can be found at the [GitHub repository](https://github.com/kyleniemeyer/irrev_mech) on GitHub. 

License
-------

`irrev_mech` is released under the modified BSD license, see LICENSE for details.

If you use this software as part of a scholarly publication, please cite the following paper in addition to the GitHub repository:

 * KE Niemeyer and CJ Sung. "Accelerating moderately stiff chemical kinetics in reactive-flow simulations using GPUs." *J. Comput. Phys.*, 256:854-871, 2014. doi:[10.1016/j.jcp.2013.09.025](http://dx.doi.org/10.1016/j.jcp.2013.09.025) 
 
Further Reading
---------------

The three reverse Arrehenius coefficients are determined using a nonlinear least-squares minimization, using the [SciPy] function `scipy.optimize.leastsq`. The fit is performed using a large number of calculated reverse rate coefficents, densely sampled over the specified temperature range. The initial guess for this minimization comes from an analytical fit to three temperatures (the high and low values of the range, plus a midpoint). See the appendix of Niemeyer and Sung's paper cited above for more details. The use of the least-squares minimization to obtain better fits was suggested by Taylor et al.:

 * BD Taylor, DA Schwer, and A Corrigan. "Implementation of Thermochemistry and Chemical Kinetics in a GPU-based CFD Code." 53rd AIAA Aerospace Sciences Meeting, Kissimmee, FL, January 2015. doi:[10.2514/6.2015-0842](http://dx.doi.org/10.2514/6.2015-0842)

Author
------

Created by [Kyle Niemeyer](http://kyleniemeyer.com). Email address: [kyle.niemeyer@gmail.com](mailto:kyle.niemeyer@gmail.com)


[Python]: http://python.org/
[SciPy]: http://scipy.org/