#!/usr/bin/env python
# -*- encoding: utf-8 -*-

"""
"""

from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'irrev_mech', '_version.py')) as version_file:
    exec(version_file.read())

# Get the long description from the relevant file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='irrev_mech',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=__version__,

    description='Chemical kinetic model irreversible reaction converter',
    long_description=long_description,

    # Author details
    author='Kyle Niemeyer',
    author_email='kyle.niemeyer@gmail.com',

    # The project's main homepage.
    url='https://github.com/kyleniemeyer/irrev_mech',

    license='MIT',

    packages=find_packages(),

    entry_points={
        'console_scripts': [
            'irrev_mech=irrev_mech.__main__:main',
        ],
    },

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 5 - Production/Stable',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],

    # What does your project relate to?
    keywords='chemical_kinetics',

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    # install_requires=['peppercorn'],
    install_requires=['numpy',
                      'scipy'],

    platforms='OS independent'
)
