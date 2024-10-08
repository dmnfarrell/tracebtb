# TracebTB

<img align="right" src=tracebtb/logo.svg width=180px>

This is a desktop application designed to allow analysis of M.bovis strain types and associated meta data. It integrates genetic,  location and cattle movement data in a single tool.

This software is written in Python. It was developed on Ubuntu linux and should run on any linux desktop. The GUI is made using the Qt toolkit using PySide2/PyQt.

## Usage

Run the command `tracebtb` from the command line to launch the application. You can provide a saved project file by using the -p option: `tracebtb -p file.tracebtb`

See the [wiki](https://github.com/dmnfarrell/tracebtb/wiki) page for documentation.

## Installation

### Linux

`pip install -e git+https://github.com/dmnfarrell/tracebtb.git#egg=tracebtb`

Notes: You may need to use pip3 on Ubuntu to ensure you use Python 3. Running this also requires you have **git** installed. The same command can be used to update to the latest version.

### WSL on Windows

The tool can be run on Windows using Windows Subsystem for Linux (WSL). This is [easy to install](https://www.omgubuntu.co.uk/how-to-install-wsl2-on-windows-10) from the Micrsoft store. (You can also run `wsl --install' from a terminal window). You then install Ubuntu 22.04 LTS from the store. Just make sure your version of Windows 10 is also reasonably up to date.

Open Ubuntu from the start menu and run it. It will bring up a terminal. This is just like being in normal Linux. You should first install some required packages:

```sudo apt install python3-pip x11-apps libqt5extras5 graphviz-dev fasttree```

You can then install the package using pip as above.

## Requirements 

* matplotlib
* pandas
* geopandas
* biopython
* pyside2 or pyqt5/PyQtWebEngine
* pygraphviz
* contextily

## Funding

This tool was developed as part of projects funded by the Irish Department of Agriculture and the Marine (DAFM).
