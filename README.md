# bTBWGStool

<img align="right" src=btbwgstool/logo.svg width=180px>

This is a desktop application designed to allow analysis of M.bovis strain types and associated meta data. It integrates genetic,  location and cattle movement data in a single tool.

This software is written in Python. It was developed on Ubuntu linux and should run on any linux desktop. The GUI is made using the Qt toolkit using PySide2/PyQt.

## Usage

Run the command `btbwgstool` from the command line.

## Data formats

Three primary kinds of data are used in this tool. You need at least spatial and genetic distance data. They can all be provided as text files with a specific format.

### spatial and other meta data

This is a table of sample names with lat and long values for the spatial locations of each sample. Normally these are the locations of farm herds (centroids of land parcels) or sett locations. They could also be capture locations of wildlife.

```
sample,lat,lon,country
A1,54.1150448498352,-6.88253943775043,Ireland
A2,54.2118128119312,-6.847691705698,Ireland
B1,54.1822520792772,-7.09087486232739,Ireland
B2,54.184003143,-7.0817350411,Ireland
D1,54.411216388944,-7.00761048343717,Ireland
...
```

### genetic distances

these should be provided as a sequence alignment of snp positions derived from a variant calling tool like SNiPgenie.

```

...
```

### movement data

Movement data can be provided as a table with the following format:


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
