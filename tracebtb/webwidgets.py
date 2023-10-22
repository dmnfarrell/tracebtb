# -*- coding: utf-8 -*-

"""
    Qt widgets module.
    Created Sep 2023
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

import sys, os, io, platform
import numpy as np
import pandas as pd
import pylab as plt
import matplotlib as mpl
import string
from .qt import *
from . import core, widgets, plotting
import geopandas as gpd

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

module_path = os.path.dirname(os.path.abspath(__file__))
iconpath = os.path.join(module_path, 'icons')
data_path = os.path.join(module_path,'data')
counties = os.path.join(data_path,'counties.shp')

class FoliumViewer(QWidget):
    """folium plot widget"""
    def __init__(self, parent=None):
        super(FoliumViewer, self).__init__(parent)
        self.main = QWebEngineView()
        l = QVBoxLayout()
        self.setLayout(l)
        l.addWidget(self.main)
        #l.addWidget(QLabel('test'))
        #self.test()
        return

    def test(self):
        code = '<html> HELLO </html>'
        self.main.setHtml(code)
        return

    def show(self):
        map = self.create_map()
        data = io.BytesIO()
        map.save(data, close_file=False)
        #print (data.getvalue().decode())
        self.main.setHtml(data.getvalue().decode())
        return

    def create_map(self, location=[54.1, -7.0]):
        """Make base map with tiles"""

        import folium
        from folium import plugins
        #print (location)

        map = folium.Map(location=location, crs='EPSG3857',tiles='openstreetmap',
                                    width=1500, height=1200)
        tile = folium.TileLayer(
            tiles = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
            attr = 'Esri',
            name = 'Esri Satellite',
            #overlay = False,
            control = True
            ).add_to(map)
        tile = folium.TileLayer(
            tiles = 'https://mt1.google.com/vt/lyrs=m&x={x}&y={y}&z={z}',
            attr = 'Google',
            name = 'Google Maps',           
            control = True
            ).add_to(map)
        tile = folium.TileLayer(
            tiles = 'https://mt1.google.com/vt/lyrs=p&x={x}&y={y}&z={z}',
            attr = 'Google',
            name = 'Google Terrain',          
            control = True
            ).add_to(map)
        tile = folium.TileLayer(
            tiles = 'https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}',
            attr = 'Google',
            name = 'Google Satellite',       
            control = True
            ).add_to(map)             
        folium.LayerControl().add_to(map)
        #GPS
        plugins.LocateControl().add_to(map)
        #c = gpd.read_file(counties).to_crs("EPSG:29902")
        #folium.GeoJson(c.to_json()).add_to(map)
        return map
    
    def plot(self, sub, lpis, borders=None, colorcol=None):

        import folium
        from branca.element import Figure
        fig = Figure(width=600, height=600)
        map = self.create_map()   

        df = sub.to_crs('EPSG:4326')
        print (df)
        c  = df.dissolve().centroid.geometry
        bounds = df.bounds            
        style1 = {'fillColor': 'blue', 'color': 'gray','weight':1}
        p = folium.GeoJson(lpis.to_crs('EPSG:4326'),style_function=lambda x:style1)    
        #map.add_child(p)
        
        #col='snp3'
        labels = df[colorcol].unique()
        #colors=plotting.gen_colors(cmap="nipy_spectral",n=len(labels))
        #colors = plotting.random_colors(n=len(labels),seed=20)
        df['color'] = df.Species.map({'Bovine':'blue','Badger':'orange'})
        #lut = dict(zip(labels, colors))
        #df['color'] = df[colorcol].map(lut)
        #c['descr'] = c.apply(lambda x: x[col]+' ',1)
        
        for i,r in df.iterrows():
            if r.geometry.is_empty: continue
            x=r.geometry.x
            y=r.geometry.y
            w=0.005
            pts = ((y-w/1.5,x-w),(y+w/1.5,x+w))        
            folium.CircleMarker(location=(y,x), radius=10, 
                            color=False,fill=True,fill_opacity=0.6,
                            fill_color=r.color,tooltip=r.Species).add_to(map)
    
            icon=folium.Icon(color='blue', icon_color='white', icon='info-sign') 

        fig.add_child(map)
        data = io.BytesIO()
        fig.save(data, close_file=False)
        self.main.setHtml(data.getvalue().decode())
        return
