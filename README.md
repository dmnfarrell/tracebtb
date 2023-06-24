# BTBwgstool

<img align="right" src=btbwgstool/logo.svg width=180px>

This is a desktop application designed to allow analysis of M.bovis strain types and associated meta data. It integrates genetic,  location and cattle movement data in a single tool.

This software is written in Python. It was developed on Ubuntu linux and should run on any linux desktop. The GUI is made using the Qt toolkit using PySide2/PyQt.

## Usage

Run the command `btbwgstool` from the command line.

## Installation

`pip install -e git+https://github.com/dmnfarrell/btbwgstool.git#egg=btbwgstool`

Notes: You may need to use pip3 on Ubuntu to ensure you use Python 3. Running this also requires you have **git** installed. The same command can be used to update to the latest version.

## Requirements 

* matplotlib
* pandas
* biopython
* toytree
* pyside2 or pyqt5

## Links

* https://toytree.readthedocs.io/
