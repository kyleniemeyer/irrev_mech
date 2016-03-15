from setuptools import setup, find_packages
setup(
    name='irrev_mech',
    version='0.1.0',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'irrev_mech = irrev_mech.irrev_mech:main',
        ],
    },
)
