#!/usr/bin/env python

"""
    btbwgstool core.
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

import sys,os,subprocess,glob,re
import time, datetime
import platform

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
config_path = os.path.join(home, '.config','tracebtb')

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
    "Carlow": "#FF5733",
    "Cavan": "#3366FF",
    "Clare": "#FF3333",
    "Cork": "#33FFCC",
    "Donegal": "#6633FF",
    "Dublin": "#FFD700",
    "Galway": "#4792BD",
    "Kerry": "#FF66CC",
    "Kildare": "#663300",
    "Kilkenny": "#6600CC",
    "Laois": "#FF6600",
    "Leitrim": "#00CC66",
    "Limerick": "#FF3399",
    "Longford": "#0000FF",
    "Louth": "#00FFFF",
    "Mayo": "#CC0000",
    "Meath": "#800000",
    "Monaghan": "#85C552",
    "Offaly": "#6600FF",
    "Roscommon": "#00FF33",
    "Sligo": "#339966",
    "Tipperary": "#FF9900",
    "Waterford": "#99FF33",
    "Westmeath": "#BDBB47",
    "Wexford": "#993399",
    "Wicklow": "#FF0000",
    "Londonderry": "white",
    "Antrim": "white",
    "Armagh": "white",
    "Tyrone": "white",
    "Fermanagh": "white",
    "Down": "white"
}