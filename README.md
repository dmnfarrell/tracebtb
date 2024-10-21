# TracebTB

<img align="right" src=tracebtb/logo.svg width=180px>

This is a set of applications designed to allow analysis of M.bovis strain types and associated meta data. It integrates genetic, location and cattle movement data for visualisation in a single tool.

This software is written in Python. It was developed on Ubuntu linux and should run on any linux desktop. The GUI was made using the Qt toolkit using PySide2/PyQt. The web interface was developed with the [Panel](https://panel.holoviz.org/) framework and [Bokeh](https://bokeh.org/) plotting library.

## Usage

### Web application

Run the command `tracebtb-web` from the command line to launch the web application. You should provide a saved project file by using the -p option: `tracebtb-web -p file.tracebtb`. This will launch a web browser window at the address localhost:5010.

See the [wiki](https://github.com/dmnfarrell/tracebtb/wiki) page for further documentation.

## Installation

### Linux

`pip install -e git+https://github.com/dmnfarrell/tracebtb.git#egg=tracebtb`

Notes: You may need to use pip3 on Ubuntu to ensure you use Python 3. Running this also requires you have **git** installed. The same command can be used to update to the latest version.

### Windows

Install python as normal. Also install the [git](https://gitforwindows.org/) program if you are to use the githb version. Then use pip from the command shell as above to install. You shouldn't need conda.

## Requirements

* matplotlib
* pandas
* geopandas
* biopython
* panel
* bokeh
* pyside2 or pyqt5/PyQtWebEngine (only for gui)

## Funding

This tool was developed as part of projects funded by the Irish Department of Agriculture and the Marine (DAFM).
