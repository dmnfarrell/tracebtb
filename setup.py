from setuptools import setup
import sys,os

with open('tracebtb/description.txt') as f:
    long_description = f.read()

setup(
    name = 'tracebtb',
    version = '0.4.0-alpha',
    description = 'TracebTB tool',
    long_description = long_description,
    url='https://github.com/dmnfarrell/tracebtb',
    license='GPL v3',
    author = 'Damien Farrell',
    author_email = 'damien.farrell@ucd.ie',
    packages = ['tracebtb'],
    package_data={'tracebtb': ['data/*.*','logo.png',
                  'description.txt']
                 },
    install_requires=['numpy',
                      'pandas',
                      'matplotlib>=3.0',
                      'biopython',
                      'geopandas',
                      'pyqt5',
                      'PyQtWebEngine',
                      'toytree',
                      'contextily'
                      'pygraphviz'
                      ],
    entry_points = {
        'console_scripts': [
            'tracebtb=tracebtb.gui:main']
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
