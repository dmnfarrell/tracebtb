from setuptools import setup
import sys,os

with open('tracebtb/description.txt') as f:
    long_description = f.read()

setup(
    name = 'tracebtb',
    version = '0.5.0',
    description = 'TracebTB tool',
    long_description = long_description,
    url='https://github.com/dmnfarrell/tracebtb',
    license='GPL v3',
    author = 'Damien Farrell',
    author_email = 'damien.farrell@ucd.ie',
    packages = ['tracebtb'],
    package_data={'tracebtb': ['data/*.*','logos/*.*',
                  'description.txt']
                 },
    install_requires=['numpy',
                      'pandas>=2.0',
                      'matplotlib==3.8.3',
                      'polars',
                      'biopython',
                      'geopandas',
                      'networkx',
                      'pyfaidx',
                      'matplotlib_scalebar',
                      'bokeh==3.4.2',
                      'panel==1.4.4'
                      ],
    entry_points = {
        'gui_scripts': [
            'tracebtb=tracebtb.weblauncher:main',
            'tracebtb-prep=tracebtb.prepare:main']
            },
    classifiers = ['Operating System :: OS Independent',
            'Programming Language :: Python :: 3.7',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics'],
    keywords = ['bioinformatics','biology','genomics']
)
