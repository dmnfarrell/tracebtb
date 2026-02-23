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
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import geopandas as gpd
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
from datetime import datetime
from collections import OrderedDict
from . import core

try:
    from bokeh.io import show
    from bokeh.models import (ColumnDataSource, GeoJSONDataSource, GMapOptions, GMapPlot, TileSource, FactorRange, Range1d,
                                HoverTool, BoxZoomTool, TapTool, CustomJS,
                                Legend, LegendItem, GlyphRenderer, ColorBar, LinearColorMapper, LabelSet, Label,
                                Arrow, NormalHead, OpenHead, VeeHead)
    from bokeh.models.glyphs import Patches, Circle
    from bokeh.layouts import layout, column
    from shapely.geometry import Polygon, Point
    from bokeh.transform import jitter, factor_cmap
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
speciesmarkers = {'Bovine':'circle','Badger':'square',
                  'Deer':'triangle','Ovine':'diamond',None:'x'}

module_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(module_path,'data')
counties_gdf = core.counties_gdf

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

def random_color():
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))

def random_circles(n=20):
    """Plot random circles"""

    x = np.random.rand(n) * 10
    y = np.random.rand(n) * 10
    radius = np.random.rand(n) * 0.5 + 0.1
    #color = np.random.choice(['red', 'blue', 'green', 'purple', 'orange'], n)
    color = [random_color() for _ in range(n)]
    source = ColumnDataSource(data=dict(x=x, y=y, radius=radius, color=color))
    p = figure(title="Random Circles", x_range=(0, 10), y_range=(0, 10))
    p.circle(x='x', y='y', radius='radius', fill_color='color', fill_alpha=0.8,
              line_color="color", source=source)
    p.sizing_mode = 'stretch_both'
    p.axis.visible = False
    p.toolbar.logo = None
    #p.background_fill_color = 'gray'
    return p

# Function to generate a polygon with random number of sides
def generate_polygon(center_x, center_y, min_radius, max_radius, num_sides):
    angles = np.linspace(0, 2 * np.pi, num_sides, endpoint=False)
    radii = np.random.uniform(min_radius, max_radius, num_sides)
    x_vertices = center_x + radii * np.cos(angles)
    y_vertices = center_y + radii * np.sin(angles)
    return x_vertices, y_vertices

def random_polygons(n=10):
    """Plot random polygons"""

    p = figure(title="Random Polygons", width=800,
               x_range=(-10, 10), y_range=(-10, 10),
                match_aspect=True)
    # Generate random polygons and plot them
    w=10;h=10
    for _ in range(n):
        num_sides = random.randint(4, 10)
        center_x = random.uniform(-w, w)
        center_y = random.uniform(-h, h)
        min_radius = random.uniform(0.5, 1.0)
        max_radius = random.uniform(1.0, 4.0)
        x, y = generate_polygon(center_x, center_y, min_radius, max_radius, num_sides)
        # Plot the polygon
        color = random_color()
        p.patch(x, y, alpha=0.6, line_width=2, color=color)
    p.sizing_mode = 'stretch_height'
    p.axis.visible = False
    p.toolbar.logo = None
    return p

def save_figure(p):
    output_file(filename="test.html", title="Static HTML file")
    save(p)

def init_figure(title=None, provider=None, width=600, height=600, sizing_mode='stretch_both'):
    """Create base figure with map"""

    box_zoom = BoxZoomTool(match_aspect=True)
    p = figure(tools=['pan,wheel_zoom,reset,save',box_zoom], active_scroll="wheel_zoom",
            x_axis_type="mercator", y_axis_type="mercator",
            title=title,
            width=width, height=height,
            match_aspect=True)
    if provider in providers:
        p.add_tile(provider, retina=True)
    p.sizing_mode = sizing_mode
    if title != None:
        p.title.text_font_size = '12pt'
    return p

def add_scalebar(p):
    """Add a scale bar to a map plot"""

    scale_source = ColumnDataSource(data={'x': [], 'y': []})
    scale_bar = p.line('x', 'y', source=scale_source, line_width=5, color='black')
    scale_label = Label(x=0, y=0, text=' ', text_font_size='11pt', text_font_style='bold')
    p.add_layout(scale_label)

    # CustomJS callback for updating scale bar on zoom/pan
    callback = CustomJS(args=dict(scale_source=scale_source, scale_label=scale_label, plot=p), code="""
        // Get the current extents of the x_range and y_range
        var x_range = plot.x_range;
        var y_range = plot.y_range;

        // Compute the scale bar length as approximately 1/4 of the x range
        var scale_length = (x_range.end - x_range.start) / 4;

        // Round the scale length to a nice number (significant digits)
        var nice_scale = Math.pow(10, Math.floor(Math.log10(scale_length)));
        scale_length = Math.round(scale_length / nice_scale) * nice_scale;

        // Compute the bottom-left corner position for the scale bar
        var x_start = x_range.start + (x_range.end - x_range.start) * 0.05;
        var x_end = x_start + scale_length;
        var y_pos = y_range.start + (y_range.end - y_range.start) * 0.05;

        // Update the scale bar data source
        scale_source.data = {
            'x': [x_start, x_end],
            'y': [y_pos, y_pos]
        };

        // Update the label to show the scale length
        scale_label.x = x_start;
        scale_label.y = y_pos + (y_range.end - y_range.start) * 0.01;
        scale_label.text = scale_length.toFixed(2)/1000 + 'km';

        scale_source.change.emit();
        scale_label.change.emit();
    """)

    # Attach the callback to the plot's x_range and y_range
    p.x_range.js_on_change('start', callback)
    p.x_range.js_on_change('end', callback)
    p.y_range.js_on_change('start', callback)
    p.y_range.js_on_change('end', callback)
    return

def plot_selection(gdf, parcels=None, provider='CartoDB Positron', col=None,
                   legend=False, legend_fontsize=12,  label_fontsize=14,
                   title=None, ms=10, lw=1.2, labels=False, scalebar=True,
                   p=None):
    """
    Plot geodataframe selections with bokeh
    Args:
        gdf: locations of samples
        parcels: land parcels geodataframe
        provider: context map provider
        ms: marker size
        p: existing figure if needed
    """

    if p == None:
        p = init_figure(title, provider)
    gdf = gdf[~gdf.geometry.is_empty]
    gdf['last_move'] = gdf['last_move'].astype(str)
    if col in gdf.columns:
        gdf = gdf.sort_values(col)
    if len(gdf) == 0:
        p = error_message(msg='no samples with locations')
        return p
    if not 'color' in gdf.columns:
        gdf['color'] = 'blue'
    gdf['marker'] = gdf.Species.map(speciesmarkers).fillna('circle')
    geojson = gdf.to_crs('EPSG:3857').to_json()
    geo_source = GeoJSONDataSource(geojson=geojson)

    #add parcel polygons if provided
    if parcels is not None and len(parcels) > 0:
        parcelsjson = parcels.to_crs('EPSG:3857').to_json()
        poly_source = GeoJSONDataSource(geojson=parcelsjson)
        r1 = p.patches('xs', 'ys', source=poly_source, fill_color='color',
                       fill_alpha=0.7, line_width=1, line_color='black', name='parcels')
        #r1.selection_glyph = p.circle(size=15, fill_color="firebrick",
        #                                 line_color="black", alpha=.6)
        #hover tool
        h1 = HoverTool(renderers=[r1],
                       tooltips=([("Herd", "@SPH_HERD_N")]))
        p.add_tools(h1)

    #sett outlines
    #plot_setts(gdf, p)
    #draw points
    r2 = p.scatter('x', 'y', source=geo_source, color='color', line_width=lw,
                   line_color='black', marker="marker", fill_alpha=0.8, size=ms, name='points')
    h2 = HoverTool(renderers=[r2], tooltips=([("Sample", "@sample"),
                                            ("Animal_id", "@Animal_ID"),
                                            ("Herd/Sett", "@HERD_NO"),
                                            ("Year", "@Year"),
                                            ("Last move", "@last_move"),
                                            ("Homebred","@Homebred"),
                                            ("Lineage", "@lineage"),
                                            ("Strain", "@short_name"),
                                            ('snp5',"@snp5"),
                                            ('snp12',"@snp12")
                                           ]))
    p.add_tools(h2)

    if legend == True and col != None:
        add_legend(p, gdf, col, legend_fontsize)

    if labels == True and parcels is not None:
        cent = tools.calculate_parcel_centroids(parcels).to_crs('EPSG:3857')
        cent['color'] = parcels.color
        labels_source = GeoJSONDataSource(geojson=cent.to_json())
        labels = LabelSet(x='x', y='y', text='SPH_HERD_N', source=labels_source,
                          #x_offset=5, y_offset=5,
                          text_align='right', background_fill_color='white', background_fill_alpha=0.8,
                          text_font_size = f"{label_fontsize}px")
        p.add_layout(labels)

    #if only one sample and no parcels set ranges
    if len(gdf)==1 and parcels is None:
        #print (len(gdf))
        geo_data = json.loads(geojson)
        coords = geo_data["features"][0]["geometry"]["coordinates"]
        x_point, y_point = coords[0], coords[1]
        radius = 2000
        p.x_range = Range1d(x_point - radius, x_point + radius)
        p.y_range = Range1d(y_point - radius, y_point + radius)

    if scalebar == True:
        add_scalebar(p)

    p.axis.visible = False
    p.toolbar.logo = None
    p.add_tools(TapTool())
    return p

def add_legend(p, gdf, col, fontsize=14):
    """
    Add legend to figure given gdf and col color.
    If the column is numeric we add a colorbar.
    """

    if gdf[col].dtype == object:
        vals = gdf[col].astype(str)
        color_map = OrderedDict(zip(vals,gdf.color))
        legend_items = []
        x = (p.x_range.end-p.x_range.start)/2
        y = (p.y_range.end-p.y_range.start)/2
        for c, color in color_map.items():
            r = p.scatter(x=[x], y=[y], color=color, size=5)
            legend_items.append(LegendItem(label=c,renderers=[r]))
            r.visible=False
        legend = Legend(items=legend_items, location="top_left", title=col)
        p.add_layout(legend, 'right')
        p.legend.label_text_font_size = f'{fontsize}pt'
    else:
        vals = gdf[col]
        clrs = gdf.sort_values(col)['color']
        from bokeh.palettes import Viridis256
        mapper = LinearColorMapper(palette=Viridis256, low=vals.min(), high=vals.max())
        color_bar = ColorBar(color_mapper=mapper, label_standoff=12, location=(0,0))
        p.add_layout(color_bar, 'right')
    return

def plot_setts(gdf, p):
    """Add circles for setts where there are badger samples"""

    b = gdf[gdf.Species=='Badger']
    grouped = b.groupby('HERD_NO')['geometry'].apply(lambda x: x.unary_union.centroid)
    centroids_gdf = gpd.GeoDataFrame(grouped, geometry='geometry', crs=gdf.crs).reset_index()
    geojson = centroids_gdf.to_crs('EPSG:3857').to_json()
    source = GeoJSONDataSource(geojson=geojson)
    r = p.scatter('x', 'y', source=source, line_width=1, color=None,
                   line_color='black', marker="circle", size=40, name='setts')
    return

def plot_gdf(gdf, p, **kwargs):
    geojson = gdf.to_json()
    source = GeoJSONDataSource(geojson=geojson)
    r = p.patches('xs', 'ys', source=source,
                  fill_alpha=0, **kwargs)
    return

def plot_counties(p):
    """Plot county borders"""

    geojson = counties_gdf.to_json()
    source = GeoJSONDataSource(geojson=geojson)
    r = p.patches('xs', 'ys', source=source,
                  line_width=2, line_color='red', fill_alpha=0)
    return p

def plot_lpis(gdf, p=None, provider='CartoDB Positron',
               labels=False, label_fontsize=14, **kwargs):
    """Plot LPIS land parcels"""

    if p is None:
        p = init_figure(None, provider)
    if gdf is None:
        return

    parcelsjson = gdf.to_crs('EPSG:3857').to_json()
    source = GeoJSONDataSource(geojson=parcelsjson)
    r = p.patches('xs', 'ys', source=source, fill_color='color',
                           line_color='black', **kwargs)
    h = HoverTool(renderers=[r], tooltips=([("herd", "@SPH_HERD_N"),
                                            ('herd size', "@COUNT"),
                                            ('size', "@Size"),
                                            ('base year', "@base_year"),
                                            ('total std reactors', "@sr_total")]))
    p.add_tools(h)
    if labels == True:
        cent = tools.calculate_parcel_centroids(gdf).to_crs('EPSG:3857')
        cent['color'] = gdf.color
        labels_source = GeoJSONDataSource(geojson=cent.to_json())
        labels = LabelSet(x='x', y='y', text='SPH_HERD_N', source=labels_source,
                          #x_offset=5, y_offset=5,
                          text_align='right', background_fill_color='white', background_fill_alpha=0.8,
                          text_font_size = f"{label_fontsize}px")
        p.add_layout(labels)
    return p

def plot_moves(p, moves, lpis_cent, limit=300, name='moves'):
    """Plot moves with bokeh - replace with network plot"""

    nh = VeeHead(size=12, fill_color='blue', fill_alpha=0.5, line_color='black')
    moves = moves[moves.geometry.notnull()].to_crs('EPSG:3857')
    groups = moves.groupby('tag')
    if len(groups) > limit:
        print ('too many moves')
        return
    for tag,t in groups:
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
                               x_start=p1[0], y_start=p1[1], x_end=p2[0], y_end=p2[1],
                               name=name))
    return p

def plot_group_symbols(gdf, p, lw=4, ms=50):
    """Plot simplified symbols for herds and setts"""

    g = gdf.groupby('HERD_NO').first().set_crs(gdf.crs)
    g['marker'] = g.Species.map(speciesmarkers)
    geojson = g.to_crs('EPSG:3857').to_json()
    geo_source = GeoJSONDataSource(geojson=geojson)
    r = p.scatter('x', 'y', source=geo_source, color='white', line_width=lw,
                   line_color='black', marker="marker", fill_alpha=0.2, size=ms)
    return p

def plot_radius(geom, p, radius_km=4, source_crs="EPSG:29902", line_width=1.5):
    """
    Takes a single row from a GeoDataFrame and plots a 4km radius.
    The ColumnDataSource will contain all columns from that row.
    """

    if geom is None or geom.is_empty:
        return p

    import pyproj
    from shapely.ops import transform
    # 1. Project the single point to meters (3857)
    transformer = pyproj.Transformer.from_crs(source_crs, "EPSG:3857", always_xy=True).transform
    point_m = transform(transformer, geom)
    data = {}
    # Update the coordinates to the projected ones for Bokeh
    data['x'] = [point_m.x]
    data['y'] = [point_m.y]
    source = ColumnDataSource(data=data)
    p.circle(
        x='x',
        y='y',
        radius=radius_km * 1000,
        source=source,
        fill_color=None,
        line_color="black",
        line_dash="dashed",
        line_width=line_width
    )
    return

def error_message(msg=''):
    """Return plot with message"""

    p = figure(x_range=(-1, 1), y_range=(-1, 1), match_aspect=True)
    x_center = (p.x_range.start + p.x_range.end) / 2
    y_center = (p.y_range.start + p.y_range.end) / 2
    label = Label(x=x_center, y=y_center, text=msg, text_align='center',
                  text_baseline='middle', text_font_size='12pt', text_font_style="bold")
    p.add_layout(label)
    p.toolbar.logo = None
    p.xaxis.visible = False
    p.yaxis.visible = False
    p.sizing_mode='stretch_both'
    return p

def split_view(gdf, col, parcels=None, provider=None, limit=9, kde=False, **kwargs):
    """Plot selection split by a column"""

    from bokeh.layouts import gridplot
    groups = gdf.groupby(col)
    if len(groups)>limit or len(groups)<2:
        p = error_message('too many groups to display')
        return p
    common = gdf[col].value_counts().index[:limit]
    l = len(common)
    if len(common) < 2: return
    nr, nc = calculate_grid_dimensions(l)
    i=0
    figures=[]
    for c, sub in groups:
        if c in common:
            title  = f'{col}={c} len={len(sub)}'
            if parcels is not None:
                pcl = parcels[parcels.SPH_HERD_N.isin(sub.HERD_NO)].copy()
            else:
                pcl = None
            p = plot_selection(sub, pcl, provider=provider, title=title, **kwargs)
            if kde == True:
                kde_plot_groups(sub, p, col, 6)
            figures.append(p)
            i+=1

    grid = gridplot(figures, ncols=nc)
    grid.sizing_mode = 'stretch_both'
    return grid

def get_timeline_data(mov, meta, limit=300):
    """Get data for timeline plot"""

    if mov is None:
        return
    cols = ['move_to','move_date','end_date','data_type','duration','sample']
    mcols = ['sample']#,'snp5','snp7']
    cols = cols+mcols
    new = []
    #default end date if no death
    end = datetime(2022, 12, 31)
    groups = mov.groupby('tag')
    if len(groups)>limit:
        return
    for tag,t in groups:
        if len(t)==1:
            t['end_date'] = end
        else:
            t = t.sort_values('move_date')
            t['end_date'] = t.move_date.shift(-1)
        t['duration'] = t.end_date-t.move_date
        #combine meta data for sample
        row = meta[meta.Animal_ID==tag].iloc[0]
        t[mcols] = row[mcols]
        #print (tag)
        #print (t[cols])
        new.append(t[cols])
    df = pd.concat(new).reset_index()
    return df

def plot_moves_timeline(df, height=300):
    """Plot movement timeline.
        df is derived from get_timeline_data
    """

    if df is None:
        return figure()
    groups = df.groupby('tag')
    source = ColumnDataSource(df)
    source.add(df.duration.astype(str), 'length')
    p = figure(y_range=groups, width=600, height=height,
               title="timeline", tools='pan,wheel_zoom,reset,save', x_axis_type="datetime")
    r = p.hbar(source=source, y="tag", left='move_date', right='end_date', height=0.8,
               line_width=0, fill_color='color')#, legend_field="move_to")
    h = HoverTool(renderers=[r], tooltips=([("Herd", "@move_to"),
                                            ("tag", "@tag"),
                                            ("sample", "@sample"),
                                            ("snp7", "@snp7"),
                                            ("time", "@length")]),
                                           )
    p.add_tools(h)
    p.toolbar.logo = None
    p.xaxis.axis_label = "Time"
    if len(groups) > 30:
        p.yaxis.visible = False
    return p

def cat_plot(df, row, col, colorcol=None, ms=5, marker='circle', width=500, height=400):
    """Categorical scatter plot"""

    from bokeh.palettes import Spectral7
    if row == None or col == None:
        return
    df = df.drop(columns='geometry').astype(str)
    if 'color' not in df.columns:
        if colorcol:
            unique_factors = df[colorcol].unique().tolist()
            color_mapper = factor_cmap(field_name=colorcol, palette=Spectral7, factors=unique_factors)
            df['color'] = color_mapper
        else:
            df['color'] = 'blue'

    source = ColumnDataSource(df)
    xrange = df.groupby(col)
    yrange = df.groupby(row)

    p = figure(x_range=xrange, y_range=yrange,
               title="Category Plot", width=width, height=height)
    r = p.scatter(x=jitter(col, width=0.2, range=p.x_range), y=jitter(row, width=0.6, range=p.y_range),
                  source=source, alpha=0.7, size=ms, color='color', marker=marker)
    if len(df)<100:
        h = HoverTool(renderers=[r], tooltips=([("Sample", "@sample"),
                                            ("Animal_id", "@Animal_ID"),
                                            ("Herd", "@HERD_NO"),
                                            ("Homebred","@Homebred"),
                                            ("Clade", "@IE_clade")
                                           ]))
        p.add_tools(h)
    p.xaxis.axis_label = col
    p.yaxis.axis_label = row
    p.xaxis.major_label_orientation = "vertical"
    p.toolbar.logo = None
    return p

def heatmap(df):
    """Dataframe heatmap in bokeh"""

    from bokeh.palettes import Blues256
    from bokeh.transform import transform
    d = df.stack()
    d.index.names = ['row','column']
    d = d.reset_index()
    d.columns = ['row','column','value']
    source = ColumnDataSource(d)
    names = list(df.index)
    mapper = LinearColorMapper(palette=Blues256, low=df.min().min(), high=df.max().max())

    p = figure(x_axis_location="above", tools="hover,save",
               x_range=names, y_range=names,
               tooltips = [('samples', '@row, @column'), ('value', '@value')])

    p.width = 600
    p.height = 600
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "9px"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = np.pi/3
    p.toolbar.logo = None

    p.rect(x='row', y='column', width=1, height=1,
           source=source,
           fill_color=transform('value', mapper),
           line_color=None)

    # Add a color bar
    color_bar = ColorBar(color_mapper=mapper, location=(0, 0))
    p.add_layout(color_bar, 'right')
    return p

def adjust_brightness(color, factor):
    # Convert the color from hex to RGB
    rgb = np.array(mcolors.hex2color(color))

    # Adjust the brightness by the given factor
    # Clipping is used to ensure values stay within [0, 1]
    adjusted_rgb = np.clip(rgb * factor, 0, 1)

    # Convert back to hex
    return mcolors.to_hex(adjusted_rgb)

def create_linear_palette(base_color, n, light_factor=1.6, dark_factor=0.7):
    """Create color palette around a base color"""

    lighter_color = adjust_brightness(base_color, light_factor)
    darker_color = adjust_brightness(base_color, dark_factor)
    # Create a colormap from the lighter to the darker color
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", [lighter_color, darker_color])
    # Generate `n` evenly spaced colors from the colormap
    colors = [mcolors.to_hex(cmap(i/n)) for i in range(n)]
    return colors

def create_custom_palette(base_hex, num_colors=9):
    """
    Create a custom palette by generating darker variations of the base color while maintaining intensity.
    :param base_hex: A string representing the base color in hex format (e.g., "#FF6347").
    :param num_colors: The number of color variations to generate.
    :return: A list of hex color codes.
    """
    # Convert hex to RGB and then to HSL
    base_rgb = mcolors.hex2color(base_hex)
    base_hsv = mcolors.rgb_to_hsv(base_rgb)

    # Create a palette by darkening the base color while maintaining saturation
    palette = []
    hue = base_hsv[0]
    for i in range(num_colors):
        # Decrease lightness to make the color darker
        factor = 1 - (i / num_colors)  # Linearly reduce brightness
        # Maintain saturation to keep the color intense
        s = base_hsv[1]
        v = base_hsv[2] * factor
        #if s>.9: s=.9
        darker_hsv = (hue, s, v)
        darker_rgb = mcolors.hsv_to_rgb(darker_hsv)
        #print (darker_hsv)
        palette.append(mcolors.rgb2hex(darker_rgb))
    return palette

def kde(gdf, N):
    """Get kde points of geodataframe"""

    from scipy.stats import gaussian_kde
    # Extract x and y coordinates from the GeoDataFrame
    g = gdf.to_crs('EPSG:3857')
    x = g.geometry.x
    y = g.geometry.y
    # Compute min and max for the grid
    offset=x.max()-x.min()
    xmin, xmax = x.min()-offset, x.max()+offset
    ymin, ymax = y.min()-offset, y.max()+offset
    # Create a grid of points where KDE will be evaluated
    X, Y = np.mgrid[xmin:xmax:N*1j, ymin:ymax:N*1j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([x, y])
    # Perform KDE
    kernel = gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z

def kde_plot(gdf, p, color='#507CBD', levels=10, alpha=0.5):
    """kde plot of points in map"""

    x, y, z = kde(gdf, 100)
    z_range = np.max(z) - np.min(z)
    #levels = int((x.max()-x.min())/20000)
    #print (levels)
    p.grid.level = "overlay"
    palette = create_custom_palette(color, num_colors=levels)
    lvl = np.linspace(np.min(z), np.max(z), levels)
    p.contour(x, y, z, lvl[1:], fill_color=palette, line_color=palette, fill_alpha=alpha)
    return

def kde_plot_groups(gdf, p, col='strain_name', min_samples=5, alpha=0.5):
    """Kde plot of separate groups"""

    for c,sub in gdf.groupby(col):
        sub = sub[~sub.geometry.is_empty]
        #remove redundant points in same herd to avoid skewed plot?
        sub = sub.drop_duplicates('HERD_NO')
        #sub = tools.remove_outliers_zscore(sub,2)
        sub = tools.remove_outliers_mahalanobis(sub, 2)
        if len(sub)<min_samples:
            continue
        clr = sub.iloc[0].color
        kde_plot(sub, p, color=clr, levels=15, alpha=alpha)
    return

def hexbin(gdf, n_bins=10, p=None):
    """Hex bin plot of point counts"""

    from bokeh.transform import linear_cmap
    from bokeh.palettes import Reds256, OrRd9
    from bokeh.util.hex import hexbin

    if p is None:
        p = init_figure(provider='CartoDB Positron')

    gdf = gdf.to_crs('EPSG:3857')
    gdf['x'] = gdf.geometry.x
    gdf['y'] = gdf.geometry.y
    x = gdf['x'].values
    y = gdf['y'].values

    width = x.max() - x.min()
    hex_size = width/n_bins
    bins = hexbin(x, y, hex_size)
    source = ColumnDataSource(data=dict(
        q=bins.q, r=bins.r, counts=bins.counts
    ))
    #color mapper for coloring hexagons based on point counts
    color_mapper = linear_cmap(field_name='counts', palette=Reds256[::-1], low=min(bins.counts), high=max(bins.counts))
    tiles = p.hex_tile(q="q", r="r", size=hex_size, line_color='black', source=source,
               fill_alpha=0.7, fill_color=color_mapper)
    hover = HoverTool(
        tooltips=[("Count", "@counts")],
        mode="mouse",
        renderers=[tiles]
    )
    p.add_tools(hover)
    return p

def scatter_pie(gdf, groupby='HERD_NO', col='snp12', colormap=None,
                radius=2000, legend=True, legend_fontsize=12, scalebar=False,
                provider='CartoDB Positron', p=None):
    """Draw wedges colored by proportion of points in a group"""

    from sklearn.cluster import KMeans
    from bokeh.transform import cumsum
    from bokeh.models import Wedge

    gdf = gdf.to_crs('EPSG:3857')
    if p == None:
        p = init_figure('', provider)
    #if groupby == 'cl':
    #    gdf = tools.spatial_cluster(gdf,2000)
    groups = gdf.groupby(groupby)

    rend = []
    for herd, group in groups:
        if len(group) <= 1:
            continue
        geo_source = GeoJSONDataSource(geojson=group.to_json())
        pt = group.union_all().centroid
        if pt.is_empty: continue
        x, y = pt.x,pt.y

        #summarize group
        summary = group.groupby(col).size()
        summary = summary.reset_index(name='value').sort_values(by=col)
        summary['angle'] = summary['value'] / summary['value'].sum() * 2 * np.pi
        summary['color'] = summary[col].map(colormap)
        summary['radius'] = len(group)*200+200
        summary['size'] = len(group)
        row = group.iloc[0]
        summary['herd'] = row.HERD_NO

        #print (summary)
        source = ColumnDataSource(summary)

        r = p.wedge(x=x, y=y, radius='radius', source=source,
                start_angle=cumsum('angle', include_zero=True), end_angle=cumsum('angle'),
                color='color', line_color="black", fill_alpha=0.5)
        rend.append(r)

    h = HoverTool(renderers=rend, tooltips=([(col, f"@{col}"),
                                            ('herd',"@herd"),
                                            ('size',"@size")]))
    p.add_tools(h)
    if legend == True and col != None:
        add_legend(p, gdf, col, legend_fontsize)
    if scalebar == True:
        add_scalebar(p)
    p.axis.visible = False
    p.toolbar.logo = None
    return p

def plot_phylogeny(tree, df, tip_size=10, lw=1, font_size='10pt',
                   tip_labels=True, labelcol='name'):
    """
    Plots a phylogenetic tree with tips colored according to metadata.

    Parameters:
    - tree (Bio.Phylo.BaseTree.Tree): Phylogenetic tree from Biopython.
    - df (pd.DataFrame): DataFrame containing metadata for each tip with a column for coloring.
    """

    def calculate_positions(tree):
        x_positions = {}
        y_positions = {}

        def assign_y_positions(clade, y_base=0):
            if clade.is_terminal():
                y_positions[clade] = y_base
                return y_base + 1
            else:
                child_positions = []
                for child in clade.clades:
                    y_base = assign_y_positions(child, y_base)
                    child_positions.append(y_positions[child])
                y_positions[clade] = sum(child_positions) / len(child_positions)
                return y_base

        def assign_x_positions(clade, x_position=0):
            x_positions[clade] = x_position
            for child in clade.clades:
                branch_length = child.branch_length if child.branch_length else 0
                assign_x_positions(child, x_position + branch_length)

        assign_y_positions(tree.root)
        assign_x_positions(tree.root)
        return x_positions, y_positions

    tip_names = [clade.name for clade in tree.get_terminals()]
    if len(tip_names)<2:
        return error_message('too few samples')
    # Get positions
    x_positions, y_positions = calculate_positions(tree)

    # Data for Bokeh
    x_h0, y_h0, x_h1, y_h1 = [], [], [], []
    x_v0, y_v0, x_v1, y_v1 = [], [], [], []

    for clade in tree.find_clades(order='level'):
        if clade.clades:
            children_y = [y_positions[child] for child in clade.clades]
            x_v0.append(x_positions[clade])
            x_v1.append(x_positions[clade])
            y_v0.append(min(children_y))
            y_v1.append(max(children_y))

            for child in clade.clades:
                x_h0.append(x_positions[clade])
                x_h1.append(x_positions[child])
                y_h0.append(y_positions[child])
                y_h1.append(y_positions[child])

    # Define ranges for plot
    x_min, x_max = min(x_positions.values()), max(x_positions.values())
    y_min, y_max = min(y_positions.values()), max(y_positions.values())
    x_padding = (x_max - x_min) * 0.1
    y_padding = (y_max - y_min) * 0.05
    x_range = (x_min - x_padding, x_max + x_padding)
    y_range = (y_min - y_padding, y_max + y_padding)
    if len(tip_names)>80:
        tip_labels = False
    if tip_labels == True:
        x_range = (x_range[0], x_range[1]*1.2)

    p = figure(tools="pan,ywheel_zoom,box_zoom,reset,save",
               x_axis_label="Branch length", y_axis_label="Clade",
               x_range=x_range,
               y_range=y_range,
               sizing_mode='stretch_both')

    # Add segments for branches
    p.segment(x_h0, y_h0, x_h1, y_h1, line_width=lw, line_color="black")
    p.segment(x_v0, y_v0, x_v1, y_v1, line_width=lw, line_color="black")

    # Add terminal nodes with color mapping
    tip_coords = [(x_positions[clade], y_positions[clade]) for clade in tree.get_terminals()]
    x_tip, y_tip = zip(*tip_coords)

    metadata = pd.DataFrame(zip(tip_names,x_tip,y_tip),columns=['name','x','y']).fillna('')
    #print (metadata)
    df = df.drop(columns='geometry').astype(str)
    metadata = metadata.merge(df, left_on='name', right_index=True)
    metadata['marker'] = metadata.Species.map(speciesmarkers).fillna('asterisk')

    source = ColumnDataSource(metadata)
    #draw tips
    r = p.scatter('x', 'y', source=source, size=tip_size, alpha=0.9, color='color',
                  marker='marker', name='tree_tips')
    # Add tip labels
    if tip_labels == True:
        labels = LabelSet(x='x', y='y', text=labelcol, level='glyph', text_color='color',
                      text_baseline='middle', x_offset=10, text_font_size=font_size,
                      source=source)
        p.add_layout(labels)
    h = HoverTool(renderers=[r], tooltips=([("name", "@name"),
                                           ("Animal_id", "@Animal_ID"),
                                            ("Herd/Sett", "@HERD_NO"),
                                            ("Year", "@Year"),
                                            ("Homebred","@Homebred"),
                                            ("Clade", "@IE_clade"),
                                            ("Strain", "@strain_name"),
                                            ('snp7',"@snp7")
                                           ]))
    p.add_tools(h)
    p.yaxis.visible = False
    p.xgrid.visible = False
    p.ygrid.visible = False
    p.toolbar.logo = None
    p.add_tools(TapTool())
    return p

'''def plot_network(G, df, pos=None, node_size=12, show_node_labels=False,
                 show_edge_labels=True):
    """
    Plot a networkx graph with edge weights displayed.
    """

    import networkx as nx
    from bokeh.plotting import from_networkx

    if pos is None:
        pos = nx.spring_layout(G)  # Default layout if not provided

    graph = from_networkx(G, pos, scale=1)

    df['marker'] = df.Species.map(speciesmarkers).fillna('asterisk')
    for key in ['sample','Animal_ID','HERD_NO','Year','snp7','marker','color']:
        if key in df.columns:
            graph.node_renderer.data_source.data[key] = [df.loc[node][key] for node in G.nodes()]

    graph.node_renderer.glyph.update(size=node_size, fill_color="color", marker='marker')

    # Extract edge weights
    edge_weights = [G[u][v].get('weight', 1) for u, v in G.edges()]  # Default weight = 1 if not set
    graph.edge_renderer.data_source.data['edge_weights'] = edge_weights

    p = figure(tools="pan,wheel_zoom,box_zoom,reset,save", tooltips=None,
               sizing_mode='stretch_both')
    p.renderers.append(graph)

    if show_edge_labels:

        # Add edge weight labels at edge midpoints
        edge_labels_source = ColumnDataSource(data={
            'x': [(pos[u][0] + pos[v][0]) / 2 for u, v in G.edges()],
            'y': [(pos[u][1] + pos[v][1]) / 2 for u, v in G.edges()],
            'weight': [str(round(w, 2)) for w in edge_weights]
        })

        edge_labels = LabelSet(x='x', y='y', text='weight', level='glyph',
                            text_color='black', text_font_size='8pt',
                            background_fill_color="white",
                            source=edge_labels_source)
        p.add_layout(edge_labels)

    # Add node labels if requested
    if show_node_labels:
        x, y = zip(*pos.values())
        node_labels = graph.node_renderer.data_source.data['sample']
        node_source = ColumnDataSource(data=dict(x=x, y=y, label=node_labels))
        node_labels = LabelSet(x='x', y='y', text='label', level='glyph',
                               text_color='black', text_baseline='middle', x_offset=0,
                               text_font_size='9pt', source=node_source)
        p.add_layout(node_labels)

    h = HoverTool(renderers=[graph], tooltips=[
        ("index", "@index"),
        ("Animal_id", "@Animal_ID"),
        ("Herd/Sett", "@HERD_NO"),
        ("Year", "@Year"),
        ('snp7', "@snp7"),
        ('Edge Weight', "@edge_weights")  # Show edge weight in hover
    ])
    p.add_tools(h)

    p.xaxis.visible = False
    p.yaxis.visible = False
    p.grid.grid_line_color = None
    p.toolbar.logo = None
    return p'''

def transform_coords_to_mercator(pos, from_epsg=29902):
    """Transform graph positions to mercartor before plotting"""

    from pyproj import Transformer
    transformer = Transformer.from_crs(from_epsg, 3857, always_xy=True)
    mercator_pos = {}
    for herd_id, (x, y) in pos.items():
        m_x, m_y = transformer.transform(x, y)
        mercator_pos[herd_id] = (m_x, m_y)
    return mercator_pos

def plot_herd_network(G, pos, p=None, line_width=1, line_color="#746969"):
    """Plot herd network with bokeh and map overlay"""

    from bokeh.models import (BoxZoomTool, HoverTool, MultiLine, Plot, Range1d,
                            ResetTool, Circle, NodesAndLinkedEdges, TapTool)
    from bokeh.plotting import from_networkx
    from bokeh.transform import factor_cmap

    pos = transform_coords_to_mercator(pos)
    if p is None:
        p = init_figure(provider='CartoDB Positron')
    # Create the Bokeh graph renderer from NetworkX
    graph_renderer = from_networkx(G, pos, scale=1, center=(0, 0))
    mapper = factor_cmap('link', palette=['red', 'gray'], factors=['direct', 'indirect'])
    graph_renderer.node_renderer.glyph = Circle(radius=5, radius_units='screen',
                                                fill_color=mapper, fill_alpha=0.6, line_color='black')
    graph_renderer.node_renderer.hover_glyph = Circle(radius=8, fill_alpha=0.6, radius_units='screen')
    graph_renderer.edge_renderer.glyph = MultiLine(line_color=line_color, line_alpha=0.6, line_width=line_width)
    graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color='green', line_width=3)
    # Add Interactive Hover Policy
    graph_renderer.inspection_policy = NodesAndLinkedEdges()

    #hover_tool = HoverTool(tooltips=[("Herd", "@herd")])
    #p.add_tools(hover_tool)
    p.renderers.append(graph_renderer)
    p.xaxis.visible = False
    p.yaxis.visible = False
    p.toolbar.logo = None
    return p

def plot_mst(dm, df, **kwargs):
    """
    Plot minimum spanning tree from dist matrix with Bokeh.
    Requires networkx and graphviz.
    """

    from networkx.drawing.nx_agraph import graphviz_layout
    import networkx as nx

    G = nx.Graph()
    for i, row in dm.iterrows():
        for j, weight in row.items():
            G.add_edge(i, j, weight=weight)

    T = nx.minimum_spanning_tree(G, algorithm='kruskal')
    pos = graphviz_layout(T)
    p = plot_network(T, df, pos, **kwargs)
    return p
