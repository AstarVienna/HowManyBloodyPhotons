#!/usr/bin/env python3
"""
How Many Bloody Photons
=======================

    $ pip install wheel twine

How to compile and put these on pip::

    $ python setup.py sdist bdist_wheel
    $ twine upload dist/*

Errors
------

- 'long_description_content_type not found':
  Can occur because the licence string is too long.
  Consider just referencing the GNU licences rather than including the full
  thing in the licence section.

"""

from setuptools import setup, find_packages

with open('README.md') as f:
    __readme__ = f.read()

with open('LICENCE') as f:
    __license__ = f.read()

with open('hmbp/version.py') as f:
    __version__ = f.readline().split("'")[1]


def setup_package():
    setup(name='HowManyPhotons',
          version=__version__,
          description="Simple photon count utility for astronomical fluxes",
          long_description=__readme__,
          long_description_content_type='text/markdown',
          author="Kieran Leschinski",
          author_email="kieran.leschinski@unive.ac.at",
          url='https://github.com/AstarVienna/HowManyBloodyPhotons',
          license="GNU General Public License v3",
          package_dir={'hmbp': 'hmbp'},
          include_package_data=True,
          packages=find_packages(exclude=('docs', 'tests')),
          install_requires=['numpy', 'astropy', 'scopesim',
                            'pyyaml', 'synphot', 'skycalc_ipy'],
          classifiers=["Programming Language :: Python :: 3",
                       "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
                       "Operating System :: OS Independent",
                       "Intended Audience :: Science/Research",
                       "Topic :: Scientific/Engineering :: Astronomy", ]
          )


if __name__ == '__main__':
    setup_package()
