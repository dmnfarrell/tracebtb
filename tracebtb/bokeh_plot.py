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

import os, sys, random, math
import string
import numpy as np
import pandas as pd
import matplotlib as mpl
import geopandas as gpd

try:
    from bokeh.io import show
    from bokeh.plotting import figure
    from bokeh.models import (ColumnDataSource, GeoJSONDataSource, GMapOptions, GMapPlot,
                               TileSource, FactorRange,
                                HoverTool, BoxZoomTool,
                                Legend, LegendItem, GlyphRenderer,
                                Arrow, NormalHead, OpenHead, VeeHead)
    from bokeh.transform import jitter, factor_cmap
    from bokeh.layouts import layout
    from bokeh.plotting import figure, output_file, save
except:
    print ('bokeh not installed')
from shapely.geometry import Polygon, Point
#from bokeh.tile_providers import Vendors
from . import tools

providers = [
    "CartoDB Positron",
    "OpenStreetMap Mapnik",
    "Esri World Imagery"
]
module_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(module_path,'data')
counties_gdf = gpd.read_file(os.path.join(data_path,'counties.shp')).to_crs("EPSG:3857")

def calculate_grid_dimensions(n):
    """Calculate the number of rows and columns that best approximate a square layout"""

    sqrt_n = math.ceil(n**0.5)
    rows = math.ceil(n / sqrt_n)
    columns = math.ceil(n / rows)
    return rows, columns

def test():
    """Test plot"""

    p = figure(tools=['pan,wheel_zoom,reset,save'])
    x = list(range(30))
    y = np.random.random(30)
    p = figure(sizing_mode="stretch_width", max_width=800, height=550)
    p.scatter(x, y, fill_color="red", size=10)
    return p

def save_figure(p):
    output_file(filename="test.html", title="Static HTML file")
    save(p)

def init_figure(title=None, provider=None):
    """Create base figure"""

    box_zoom = BoxZoomTool(match_aspect=True)
    p = figure(tools=['pan,wheel_zoom,reset,save',box_zoom], active_scroll="wheel_zoom",
            x_axis_type="mercator", y_axis_type="mercator",
            title=title,
            width=600, height=600,
            match_aspect=True)
    if provider in providers:
        p.add_tile(provider, retina=True)
    p.sizing_mode = 'stretch_both'
    return p

def plot_selection(gdf, parcels=None, provider='CartoDB Positron', col=None,
                   legend=False, title=None, ms=10):
    """
    Plot geodataframe selections with bokeh
    Args:
        gdf: locations of samples
        parcels: land parcels geodataframe
        provider: context map provider
        ms: marker size
    """

    #create figure
    p = init_figure(title, provider)

    gdf = gdf[~gdf.geometry.is_empty]
    if len(gdf) == 0:
        return p
    geojson = gdf.to_crs('EPSG:3857').to_json()
    geo_source = GeoJSONDataSource(geojson=geojson)

    #add parcel polygons if provided
    if parcels is not None and len(parcels) > 0:
        parcelsjson = parcels.to_crs('EPSG:3857').to_json()
        poly_source = GeoJSONDataSource(geojson=parcelsjson)
        r1 = p.patches('xs', 'ys', source=poly_source, fill_color='color',
                       fill_alpha=0.5, line_width=1, line_color='black')
        #hover tool
        h1 = HoverTool(renderers=[r1], tooltips=([("Herd", "@SPH_HERD_N")
                                               ]), mode='vline')
        #p.add_tools(h1)

    #draw points
    r2 = p.scatter('x', 'y', source=geo_source, color='color',
                   line_color='black', marker="marker", fill_alpha=0.5, size=ms)
    h2 = HoverTool(renderers=[r2], tooltips=([("Sample", "@sample"),
                                            ("Animal_id", "@Animal_ID"),
                                            ("Herd", "@HERD_NO"),
                                            ("Homebred","@Homebred"),
                                            ("Clade", "@IE_clade")
                                           ]))
    p.add_tools(h2)

    if legend == True and col != None:
        color_map = dict(zip(gdf[col],gdf.color))
        legend_items = []
        x = (p.x_range.end-p.x_range.start)/2
        y = (p.y_range.end-p.y_range.start)/2
        for c, color in color_map.items():
            r = p.scatter(x=[x], y=[y], color=color, size=5)
            legend_items.append(LegendItem(label=c,renderers=[r]))
            r.visible=False
        legend = Legend(items=legend_items, location="top_left", title=col)
        p.add_layout(legend, 'right')

    p.axis.visible = False
    p.toolbar.logo = None
    return p

def plot_counties(p):

    geojson = counties_gdf.to_json()
    source = GeoJSONDataSource(geojson=geojson)
    r = p.patches('xs', 'ys', source=source,
                  fill_color=None, line_width=1, line_color='gray')
    return

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

def split_view(gdf, col, parcels=None, provider=None, limit=8):
    """Plot selection split by a column"""

    from bokeh.layouts import gridplot

    common = gdf[col].value_counts().index[:limit]
    l = len(common)
    if len(common) < 2: return
    nr, nc = calculate_grid_dimensions(l)
    i=0
    figures=[]
    for c, sub in gdf.groupby(col):
        if c in common:
            title  = f'{col}={c} len={len(sub)}'
            if parcels is not None:
                pcl = parcels[parcels.SPH_HERD_N.isin(sub.HERD_NO)].copy()
            else:
                pcl = None
            s = plot_selection(sub, pcl, provider=provider, title=title)
            figures.append(s)
            i+=1

    grid = gridplot(figures, ncols=nc)
    grid.sizing_mode = 'stretch_both'
    return grid

def plot_moves_timeline(mov, herdcolors):
    """Plot movement timeline"""

    if mov is None:
        return figure()
    cols = ['move_to','move_date','end_date','data_type','duration']
    new = []
    for tag,t in mov.groupby('tag'):
        t=t.sort_values('move_date')
        t['end_date'] = t.move_date.shift(-1)
        t['duration'] = t.end_date-t.move_date
        #print (t[cols])
        new.append(t[cols])
    df = pd.concat(new)
    df['color'] = df.move_to.map(herdcolors).fillna('grey')
    groups = df.groupby('tag')
    source = ColumnDataSource(df)
    source.add(df.duration.astype(str), 'length')
    p = figure(y_range=groups, width=600, height=300,
               title="timeline", tools='pan,wheel_zoom,reset,save', x_axis_type="datetime")
    r = p.hbar(source=source, y="tag", left='move_date', right='end_date', height=0.8,
               line_width=0, fill_color='color')#, legend_field="move_to")
    h = HoverTool(renderers=[r], tooltips=([("Herd", "@move_to"),
                                            ("tag", "@tag"),
                                            ("time", "@length")]),
                                           )
    p.add_tools(h)
    p.toolbar.logo = None
    p.xaxis.axis_label = "Time"
    if len(groups) > 30:
        p.yaxis.visible = False
    return p

def cat_plot(df, row, col, colorcol=None):
    """Categorical scatter plot"""

    from bokeh.palettes import Spectral7
    if row == None or col == None:
        return
    df = df.drop(columns='geometry').astype(str)
    source = ColumnDataSource(df)
    xrange = df.groupby(col)
    yrange = df.groupby(row)
    if colorcol:
        unique_factors = df[colorcol].unique().tolist()
        color_mapper = factor_cmap(field_name=colorcol, palette=Spectral7, factors=unique_factors)
        fill_color = color_mapper
    else:
        fill_color = "blue"

    p = figure(width=600, height=300, x_range=xrange, y_range=yrange,
               title="Category Plot")
    r = p.scatter(x=jitter(col, width=0.2, range=p.x_range), y=jitter(row, width=0.6, range=p.y_range),
                  source=source, alpha=0.8, color=fill_color)#, legend_field=colorcol)
    h = HoverTool(renderers=[r], tooltips=([("Sample", "@sample"),
                                            ("Animal_id", "@Animal_ID"),
                                            ("Herd", "@HERD_NO"),
                                            ("Homebred","@Homebred")
                                           ]))
    #p.add_tools(h)
    p.xaxis.axis_label = col
    p.yaxis.axis_label = row
    p.xaxis.major_label_orientation = "vertical"
    p.toolbar.logo = None
    return p
