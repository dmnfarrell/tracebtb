#!/usr/bin/env python

"""
    tracebtb core module
    Created Jan 2022
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

import sys,os,subprocess,glob,re,json
import time, datetime
import platform
import geopandas as gpd

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
data_path = os.path.join(module_path,'data')

if platform.system() == 'Windows':
    config_path = os.path.join(os.environ['APPDATA'], 'tracebtb')
else:
    config_path = os.path.join(home, '.config','tracebtb')
if not os.path.exists(config_path):
    os.makedirs(config_path, exist_ok=True)

configfile = os.path.join(config_path, 'settings.json')
defaults = {'dashboard':{'lpis_master_file':'','tree_file':None}}
if not os.path.exists(configfile):
    with open(configfile, "w") as outfile:
        json.dump(defaults, outfile)
    treefile = None

"""Config class to keep settings globally"""
class Config:
    def __init__(self):
        self.settings = {}

    def load(self, filename):
        """Check settings file"""

        defaults = {'dashboard':{'lpis_master_file':'','tree_file':None}}
        with open(filename) as f:
            self.settings = json.load(f)['dashboard']
        for key in defaults['dashboard']:
            if not key in self.settings:
                self.settings[key] = defaults['dashboard'][key]
        return

# Instantiate once here
config = Config()
config.load(configfile)

defaultfont = 'Lato'
defaults = {
            'FONT' :defaultfont,
            'FONTSIZE' : 10,
            'TIMEFORMAT' :'%m/%d/%Y',
            'ICONSIZE' : 28,
            'DPI' : 100,
            'THREADS': 4,
            'FACECOLOR': '#FAFAF6'
         }
#populate current class variable
for k in defaults:
    vars()[k] = defaults[k]

county_colors = {
    "Antrim": "#A0522D",   # Sienna
    "Armagh": "#6B8E23",   # Olive Drab
    "Carlow": "#4682B4",   # Steel Blue
    "Cavan": "#CD5C5C",    # Indian Red
    "Clare": "#DAA520",    # Goldenrod
    "Cork": "#9370DB",     # Medium Purple
    "Derry": "#5F9EA0",    # Cadet Blue
    "Donegal": "#D2691E",  # Chocolate
    "Down": "#556B2F",     # Dark Olive Green
    "Dublin": "#8FBC8F",   # Dark Sea Green
    "Fermanagh": "#B0C4DE",# Light Steel Blue
    "Galway": "#BC8F8F",   # Rosy Brown
    "Kerry": "#B8860B",    # Dark Goldenrod
    "Kildare": "#8B4513",  # Saddle Brown
    "Kilkenny": "#778899", # Light Slate Gray
    "Laois": "#4682B4",    # Steel Blue
    "Leitrim": "#8FBC8F",  # Dark Sea Green
    "Limerick": "#DAA520", # Goldenrod
    "Longford": "#B0C4DE", # Light Steel Blue
    "Louth": "#9370DB",    # Medium Purple
    "Mayo": "#8B4513",     # Saddle Brown
    "Meath": "#A0522D",    # Sienna
    "Monaghan": "#6B8E23", # Olive Drab
    "Offaly": "#BC8F8F",   # Rosy Brown
    "Roscommon": "#D2691E",# Chocolate
    "Sligo": "#556B2F",    # Dark Olive Green
    "Tipperary": "#B8860B",# Dark Goldenrod
    "Tyrone": "#778899",   # Light Slate Gray
    "Waterford": "#CD5C5C",# Indian Red
    "Westmeath": "#5F9EA0",# Cadet Blue
    "Wexford": "#A0522D",  # Sienna
    "Wicklow": "#6B8E23"   # Olive Drab
}

clade_colors = {
    1: "#FF5733",  # Fiery Red
    2: "#33FF57",  # Lime Green
    3: "#3357FF",  # Royal Blue
    4: "#FF33A1",  # Pink
    5: "#FFA833",  # Orange
    6: "#A833FF",  # Purple
    7: "#33FFF3",  # Aqua
    8: "#FF5733",  # Coral
    9: "#57FF33",  # Bright Green
    10: "#5733FF", # Electric Blue
    11: "#FF33FF", # Magenta
    12: "#FFAA33", # Amber
    13: "#33A8FF", # Sky Blue
    14: "#8D33FF", # Violet
    15: "#FF3388", # Deep Pink
    16: "#33FFAA", # Mint Green
}

host_colors = {
    'Bovine': '#75A9CB',
    'Badger': '#FF5733',
    'Ovine': '#FFA833',
    'Deer': '#33FF57',
    'Cat': '#8D33FF',
    'Alpaca': '#6D3075',
    'Llama': '#FFCC00'
}

counties_gdf = gpd.read_file(os.path.join(data_path,'counties.shp'))
counties_gdf['geometry'] = counties_gdf.to_crs("EPSG:3857").geometry.simplify(300)
#hex grid of ireland
#iregrid = tools.get_irish_grid()

def git_version() -> str:
    """Get get version"""
    return subprocess.check_output(['git','describe','--tags']).decode('ascii').strip()

