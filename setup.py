#!/usr/bin/env python

#
# Setup script for atomneb
#

import os
import codecs
try:
      from setuptools import setup
except ImportError:
      from distutils.core import setup

import atomneb

with codecs.open('README.rst', 'r', 'utf-8') as fd:
    setup(name='atomneb',
          version=atomneb.__version__,
          description = 'atomneb: Python Package for Atomic Data of Ionized Nebulae',
          long_description=fd.read(),
          author='Ashkbiz Danehkar',
          author_email='ashkbiz.danehkar@students.mq.edu.au',
          url='https://atomneb.github.io/AtomNeb-py/',
          download_url = 'https://github.com/atomneb/AtomNeb-py',
          keywords = ['atomneb', 'Chianti atomic database', 'atomic datasets', 'recombination lines', 'collisionally excited lines', 'recombination coefficients'],
          license='http://www.gnu.org/licenses/gpl.html',
          platforms=['any'],
          packages=['atomneb'],
          #package_data={'atomneb': ['*.txt', 'text/*.readme']},
          data_files = [("", ["LICENSE"])],
          install_requires=['numpy','pandas','astropy'],
         )

