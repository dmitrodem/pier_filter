#!/usr/bin/env python2

from setuptools import setup

setup(name = 'mode_matching',
      version = '0.1',
      description = 'Mode matching method to evaluate resonant filters',
      url = 'http://localhost',
      author = 'Dmitriy Dyomin',
      author_email = 'dmitrodem@gmail.com',
      license = 'MIT',
      packages = ['mode_matching'],
      install_requires = [
          'numpy',
          'scipy',
          'pint',
          'matplotlib',
          'progressbar'
          ],
      entry_points = {
          'console_scripts': [
              'filter_example=mode_matching.Filter:main'
          ],
      },
      zip_safe = False)
