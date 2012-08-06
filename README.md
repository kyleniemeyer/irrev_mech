irrev_mech
=======

This utility converts a Chemkin-format reaction mechanism with reversible reactions to one with only irreversible reactions.

It does this by fitting reverse Arrhenius coefficients based on reverse reaction rates at three temperature points.

`irrev_mech` requires [Python] 2.x, but no additional libraries (i.e. Chemkin, Cantera). It has only been tested on [Python] 2.7, but probably works on earlier versions.

[Python]: http://python.org/

Usage
-------

From the command line, use `python irrev_mech.py mechname thermname` where `mechname` and `thermname` are the names of the mechanism file and thermodynamic database, e.g.:

    $ python irrev_mech.py mech.dat therm.dat

You can also run `irrev_mech` without a thermodynamic database if the information is held in the mechanism file (after the species are declared), e.g.:

    $ python irrev_mech.py mech.dat

The new mechanism has the name `mech_irrev.txt`.

The default temperatures used for parameter fitting are 1000 K, 1750 K, and 2500 K. These can be changed by modifying the `Tfit` variable in `convert_mech_irrev()`.

Misc
-------

The most up-to-date version of `irrev_mech` can be found at the [GitHub repository](https://github.com/kyleniemeyer/irrev_mech) on GitHub.

License
-------

`irrev_mech` is released under the modified BSD license, see LICENSE for details.