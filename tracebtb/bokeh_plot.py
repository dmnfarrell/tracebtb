"""
    Bokeh plotting methods
    Created May 2024
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

import os, sys, random
import string
import numpy as np
import pandas as pd
import matplotlib as mpl

try:
    from bokeh.io import show
    from bokeh.plotting import figure
    from bokeh.models import (ColumnDataSource, GeoJSONDataSource, GMapOptions, GMapPlot, TileSource,
                            HoverTool, BoxZoomTool,
                            Arrow, NormalHead, OpenHead, VeeHead)
    from bokeh.models.glyphs import Patches, Circle
    from bokeh.layouts import layout
    from bokeh.plotting import figure, output_file, save
except:
    print ('bokeh not installed')
from shapely.geometry import Polygon, Point
#from bokeh.tile_providers import Vendors
from . import tools

providers = [
    "CartoDB Positron",
    "CartoDB Positron retina",
    "Stamen Terrain",
    "Stamen Terrain retina",
    "Stamen Toner",
    "Stamen Toner Background",
    "Stamen Toner Labels",
    "OpenStreetMap Mapnik",
    "Esri World Imagery"
]

def test():
    """Test plot"""

    p = figure(tools=['pan,wheel_zoom,reset,save'])
    x = list(range(20))
    y = np.random.random(20)
    p = figure(sizing_mode="stretch_width", max_width=800, height=550)
    p.scatter(x, y, fill_color="red", size=10)
    return p

def save_figure(p):
    output_file(filename="test.html", title="Static HTML file")
    save(p)

def plot_selection(gdf, parcels=None, col=None, provider='CartoDB Positron'):
    """Plot dataframe selections with bokeh"""

    geojson = gdf.to_crs('EPSG:3857').to_json()
    parcelsjson = parcels.to_crs('EPSG:3857').to_json()

    # Create GeoJSONDataSources
    geo_source = GeoJSONDataSource(geojson=geojson)
    poly_source = GeoJSONDataSource(geojson=parcelsjson)

    # Create figure
    box_zoom = BoxZoomTool(match_aspect=True)
    p = figure(tools=['pan,wheel_zoom,reset,save',box_zoom], active_scroll="wheel_zoom",
               x_axis_type="mercator", y_axis_type="mercator",
               match_aspect=True)
    p.add_tile(provider, retina=True)

    # Add polygons and points to the plot
    r1 = p.patches('xs', 'ys', source=poly_source, fill_color='color',
                   fill_alpha=0.5, line_width=1, line_color='black')
    r2 = p.scatter('x', 'y', source=geo_source, color='color',
                   line_color='black', marker="marker", fill_alpha=0.5, size=10)#, legend_label=col)
    #hover tool
    h1 = HoverTool(renderers=[r1], tooltips=([("Herd", "@SPH_HERD_N")
                                           ]), mode='vline')
    h2 = HoverTool(renderers=[r2], tooltips=([("Sample", "@sample"),
                                            ("Animal_id", "@Animal_ID"),
                                            ("Herd", "@HERD_NO"),
                                            ("Homebred","@Homebred"),
                                            ("Clade", "@IE_clade")
                                           ]), mode='vline')
    #p.add_tools(h1)
    p.add_tools(h2)
    p.axis.visible = False
    p.toolbar.logo = None
    return p

def plot_moves(p, moves, lpis_cent):
    """Plot moves with bokeh)"""

    nh = VeeHead(size=12, fill_color='blue', fill_alpha=0.5, line_color='black')
    moves = moves[moves.geometry.notnull()].to_crs('EPSG:3857')
    #print (moves)
    for tag,t in moves.groupby('tag'):
        if t is not None:
            #print (t)
            moved = lpis_cent[lpis_cent.SPH_HERD_N.isin(t.move_to)].to_crs('EPSG:3857')
            coords = tools.get_coords_data(t)
            if len(coords)>0:
                mlines = gpd.GeoDataFrame(geometry=coords)
                #print (t)
                for i,l in mlines.iterrows():
                    #print (list(l.geometry.coords))
                    p1 =  l.geometry.coords[0]
                    p2 =  l.geometry.coords[1]
                    p.add_layout(Arrow(end=nh, line_color='black', line_dash=[10, 5],
                               x_start=p1[0], y_start=p1[1], x_end=p2[0], y_end=p2[1]))

    return p
