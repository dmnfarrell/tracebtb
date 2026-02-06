# TracebTB

<img align="right" src=tracebtb/logos/logo_alt.svg width=180px>

This is an experimental web application designed to allow analysis of *M.bovis* strain types and associated meta data. It integrates genetic, location and cattle movement data.

This software is written in Python. It was developed on Ubuntu. The web interface was developed with the [Panel](https://panel.holoviz.org/) framework and [Bokeh](https://bokeh.org/) plotting library.

## Usage

Run the command `tracebtb` from the command line to launch the web application. You should provide a saved project file by using the -p option: `tracebtb -p file.tracebtb`. This will launch a web browser window at the address localhost:5010.

See the [wiki](https://github.com/dmnfarrell/tracebtb/wiki) page for further documentation.

## Installation

### Linux

`pip install -e git+https://github.com/dmnfarrell/tracebtb.git#egg=tracebtb`

Notes: You may need to use pip3 on Ubuntu to ensure you use Python 3. Running this also requires you have **git** installed. The same command can be used to update to the latest version.

### Windows

It is recommended to use WSL and follow Linux instructions.

## Requirements

* matplotlib
* pandas
* geopandas
* biopython
* panel
* bokeh
* networkx
* polars
* deltalake

## Funding

This tool was developed as part of projects funded by the Irish Department of Agriculture and the Marine (DAFM).

## References

* O’Shaughnessy, J., Harvey, N., Byrne, B. et al. The genomic diversity and spatial patterns of Mycobacterium bovis in Ireland revealed by whole genome sequencing. Ir Vet J 79, 6 (2026). https://doi.org/10.1186/s13620-025-00324-0
* Harvey, N., McGrath, G., O’Shaughnessy, J. et al. Integrating whole-genome sequencing and epidemiology to characterise Mycobacterium bovis transmission in Ireland: a proof of concept. Ir Vet J 79, 3 (2026). https://doi.org/10.1186/s13620-025-00321-3