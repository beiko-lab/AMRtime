#!/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

readme = open('README.rst').read()
doclink = """
Documentation
-------------

The full documentation is at http://AMRtime.rtfd.org."""
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

setup(
    name='AMRtime',
    version='0.1.0',
    description='Metagenomic AMR detection using hierarchical machine learning models',
    long_description=readme + '\n\n' + doclink + '\n\n' + history,
    author='Finlay Maguire',
    author_email='finlaymaguire@gmail.com',
    url='https://github.com/fmaguire/AMRtime',
    packages=[
        'AMRtime',
    ],
    package_dir={'AMRtime': 'amrtime'},
    include_package_data=True,
    install_requires=[
    ],
    license='MIT',
    zip_safe=False,
    keywords='AMRtime',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
)
