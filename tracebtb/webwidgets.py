# -*- coding: utf-8 -*-

"""
    Qt web widgets module.
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
from . import core, tools, widgets, plotting
import geopandas as gpd
import shapely

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

module_path = os.path.dirname(os.path.abspath(__file__))
iconpath = os.path.join(module_path, 'icons')
data_path = os.path.join(module_path,'data')
borders = gpd.read_file(os.path.join(data_path,'counties.shp'))

def get_bounds(gdf):
    """Get bounding coords for points in gdf"""

    bounding_box = gdf.bounds
    minx = bounding_box['minx'].min()
    miny = bounding_box['miny'].min()
    maxx = bounding_box['maxx'].max()
    maxy = bounding_box['maxy'].max()
    return minx, miny, maxx, maxy

def create_map(location=[-51, 8]):
    """Make a map"""
    map = folium.Map(location=[c.y, c.x], crs='EPSG3857',tiles='openstreetmap',
                          width=600, height=600 ,max_bounds=True)
    return map

def add_tiles(map):
    """Add tile layers to base map"""

    import folium
    from folium import plugins

    tile = folium.TileLayer(
        tiles = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
        attr = 'Esri',
        name = 'Esri Satellite',
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
    #style2 = {'fillColor': '#00000000', 'color': 'gray','weight':1}
    return map

def plot_moves_folium(moves, lpis_cent, map):
    """Plot moves in folium"""

    import folium
    if moves is None:
        return
    #fg = folium.FeatureGroup("Moves")
    moves = moves[moves.geometry.notnull()].to_crs('EPSG:4326')

    i=0
    for tag,t in moves.groupby('Animal_ID'):
        if t is not None:
            moved = lpis_cent[lpis_cent.SPH_HERD_N.isin(t.move_to)].to_crs('EPSG:4326')
            coords = tools.get_coords_data(t)
            #reverse coords to long-lat
            coords = coords.map(lambda l: shapely.ops.transform(lambda x, y: (y, x), l))

            if len(coords)==0:
                continue
            for c in coords:
                loc = list(c.coords)
                tip=tag
                l=folium.PolyLine(locations=loc, color='black', weight=1, tooltip=tip).add_to(map)

            for i,r in moved.iterrows():
                if r.geometry is None or r.geometry.is_empty: continue
                x=r.geometry.x
                y=r.geometry.y
                tip =  """{}<br>""".format(r.SPH_HERD_N)
                folium.CircleMarker(location=(y,x), radius=5, weight=1,
                                  color='black',fill=False,
                                  tooltip=tip).add_to(map)

            i+=1

    #fg.add_to(map)
    return

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

    def clear(self):
        """Clear"""

        code = '<html>  </html>'
        self.main.setHtml(code)
        return

    def test_map(self):
        import folium
        map = self.create_map()
        #b = folium.GeoJson(borders[:5].to_crs('EPSG:4326'),name='borders',
        #                  smooth_factor=3)
        #map.add_child(b)
        data = io.BytesIO()
        map.save(data, close_file=False)
        self.main.setHtml(data.getvalue().decode())

    def show(self):
        """Show base map"""

        map = create_map()
        data = io.BytesIO()
        map.save(data, close_file=False)
        #print (data.getvalue().decode())
        self.main.setHtml(data.getvalue().decode())
        return

    def plot(self, sub, parcels, neighbours=None, moves=None, lpis_cent=None,
             colorcol=None, parcelscol=None, cmap='Set1'):
        """Plot selected"""

        import folium
        from branca.element import Figure

        df = sub.to_crs('EPSG:4326')
        minx, miny, maxx, maxy = get_bounds(df)
        pad=.05
        bbox = [(miny-pad,minx-pad),(maxy+pad,maxx+pad)]
        c = df.dissolve().centroid.geometry
        bounds = df.bounds
        #fig = Figure()
        map = folium.Map(location=[c.y, c.x], crs='EPSG3857',tiles='openstreetmap',
                            width=1500, height=1200 ,max_bounds=True, control_scale = True)

        basestyle = {'fillColor': 'blue', 'fillOpacity': 0.4, 'color': 'gray','weight':1}
        def style_func(feature):
            return {
                'fillColor': feature['properties']['color'],
                'weight': 1,
                'color': 'gray',
                'fillOpacity': 0.4,
            }
        if parcels is not None:
            if parcelscol not in [None, '']:
                labels = parcels[parcelscol].unique()
                colors = plotting.gen_colors(cmap=cmap,n=len(labels))
                lut = dict(zip(labels, colors))
                parcels['color'] = parcels[parcelscol].map(lut)
                style = lambda x:style_func(x)
            else:
                style = lambda x:basestyle
            tooltip = folium.features.GeoJsonTooltip(fields=['SPH_HERD_N'], aliases=['Herd No.'])
            p = folium.GeoJson(parcels.to_crs('EPSG:4326'),style_function=style,
                            tooltip=tooltip, name='parcels')
            p.add_to(map)

        if neighbours is not None:
            try:
                tooltip = folium.features.GeoJsonTooltip(fields=['SPH_HERD_N'], aliases=['Herd No.'])
                p = folium.GeoJson(neighbours.to_crs('EPSG:4326'),style_function=style,
                                tooltip=tooltip, name='neighbours')
                p.add_to(map)
            except Exception as e:
                print (e)
        #colors = plotting.random_colors(n=len(labels),seed=20)
        if colorcol == None or colorcol == '':
            df['color'] = df.Species.map({'Bovine':'blue','Badger':'orange'})
        else:
            labels = df[colorcol].unique()
            colors = plotting.gen_colors(cmap=cmap,n=len(labels))
            lut = dict(zip(labels, colors))
            df['color'] = df[colorcol].map(lut)

        for i,r in df.iterrows():
            if r.geometry == None or r.geometry.is_empty: continue
            x=r.geometry.x
            y=r.geometry.y
            w=0.005
            pts = ((y-w/1.5,x-w),(y+w/1.5,x+w))
            #tip =  """{}<br>{}<br>{}<br>{}<br>snp3={}<br>st={}""".format(
            #    r.Animal_ID,r.Species,r.Aliquot,r.HERD_NO,r.snp3,r.strain_name)
            tip = self.get_tooltip(r)
            folium.CircleMarker(location=(y,x), radius=10,
                            color=False,fill=True,fill_opacity=0.6,
                            fill_color=r.color,tooltip=tip).add_to(map)

            icon=folium.Icon(color='blue', icon_color='white', icon='info-sign')

        folium.FitBounds(bbox).add_to(map)
        if moves is not None:
            plot_moves_folium(moves, lpis_cent, map)
        add_tiles(map)
        data = io.BytesIO()
        map.save(data, close_file=False)
        self.main.setHtml(data.getvalue().decode())
        return map

    def get_tooltip(self, x):
        """Tooltip"""

        cols = ['Animal_ID','Species','Aliquot','HERD_NO','snp3','strain_name']
        tip = ''
        for i,val in x.items():
            if i in cols:
                tip+='{}: {}<br>'.format(i,val)
        return tip

    def get_parcel_tooltip(self, x):

        cols = ['SPH_HERD_N']
        tip = ''
        for i,val in x.items():
            if i in cols:
                tip+='{}: {}<br>'.format(i,val)
        return tip

    def screen_capture(self, filename=None):
        """Capture map"""

        if filename is None:
            filename = 'capture.png'
        size = self.main.contentsRect()
        img = QPixmap(size.width(), size.height())
        self.main.render(img)
        img.save(filename)
        #browser.close()
        return filename