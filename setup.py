from setuptools import setup
import sys,os

with open('btbwgstool/description.txt') as f:
    long_description = f.read()

setup(
    name = 'btbwgstool',
    version = '0.1.0',
    description = 'BTBGenIE tool',
    long_description = long_description,
    url='https://github.com/dmnfarrell/btbwgstools',
    license='GPL v3',
    author = 'Damien Farrell',
    author_email = 'farrell.damien@gmail.com',
    packages = ['btbwgstool'],
    package_data={'btbwgstool': ['data/*.*','logo.png',
                  'description.txt']
                 },
    install_requires=['numpy',
                      'pandas',
                      'matplotlib>=3.0',
                      'pyqt5,
                      'toytree',
                      ],
    entry_points = {
        'console_scripts': [
            'btbwgstool=btbwgstool.gui:main']
            },
    classifiers = ['Operating System :: OS Independent',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.7',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics'],
    keywords = ['bioinformatics','biology','genomics']
)
