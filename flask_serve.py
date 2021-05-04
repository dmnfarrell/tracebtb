#!/usr/bin/env python

"""
    BTBgenie prototype web app.
    Created Mar 2021
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

import sys, os, io
import numpy as np
import pandas as pd
import string
#from . import tools, widgets, tables
from flask import Flask, render_template, request
from wtforms import Form, TextField, validators, StringField, SelectField, FloatField
import sqlite3
import toytree, toyplot

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__))
DATABASE = '../notebooks/'
species_colors = {'Cow':'green', 'Badger':'blue', 'Deer':'red'}
sb_colors = {'SB0054':'blue','SB0041':'green'}
clade_colors = {2:'yellow',3:'green'}
cmaps = {'species': species_colors,'spoligotype':sb_colors,'clade':clade_colors}
providers = ['CARTODBPOSITRON','STAMEN_TERRAIN','STAMEN_TONER','OSM','ESRI_IMAGERY']

tree_style = {
    "layout":'r',
    "edge_type": 'p',
    "edge_style": {
        "stroke": 'black',
        "stroke-width": 2,
    },
    "tip_labels": True,
    "tip_labels_align": True,
    "tip_labels_colors": 'black',
    "tip_labels_style": {
        "font-size": "14px"
    },
    "node_labels": False,
    "node_sizes": 10,
    "node_colors": toytree.colors[2],
    "node_markers":"c",
    "use_edge_lengths":True,
}

webapp = Flask(__name__)

class ControlsForm(Form):
    colorby = SelectField('Color By:', choices=['clade','spoligotype','species'])
    tile = SelectField('Tile:', choices=providers)
    database = TextField('database', default='localhost')
    animal = SelectField('Animal:', choices=[])
    lat = FloatField('Lat',default=-6.2)
    long = FloatField('Lat',default=53.5)

def get_db():
    db = getattr(g, '_database', None)
    if db is None:
        db = g._database = sqlite3.connect(DATABASE)
    return db

def get_dataframe_db(conn):

    df = pd.read_sql_query('SELECT * FROM ', conn)
    return df

def help_msg():
    msg = '<a>BTBgenie sample web interface</a><br>'
    #msg += '<a href="%s"> see help page</a>' %wikipage
    return msg

def show_dataframe(df, plot):
    """Show points from a dataframe"""

    for i,r in df.iterrows():
        popup = get_popup(r)
        cm.add_to(m)
    return

def wgs84_to_web_mercator(df, lon="LON", lat="LAT"):

      k = 6378137
      df["x"] = df[lon] * (k * np.pi/180.0)
      df["y"] = np.log(np.tan((90 + df[lat]) * np.pi/360.0)) * k
      return df

def bokeh_map(df=None, long=None, lat=None,
              tile_provider='CARTODBPOSITRON', colorby='species'):
    """Bokeh map"""

    from bokeh.plotting import figure, show
    from bokeh.embed import components
    from bokeh.tile_providers import get_provider
    from bokeh.models import (DataTable, GeoJSONDataSource, ColumnDataSource, HoverTool,
            CustomJS, MultiSelect, Dropdown, Div)
    from bokeh.layouts import column, row
    import json
    from collections import OrderedDict

    tile_provider = get_provider(tile_provider)
    tools = "pan,wheel_zoom,hover,tap,lasso_select,reset,save"
    sizing_mode='stretch_both'

    # range bounds supplied in web mercator coordinates
    k = 6378137
    if lat == None:
        lat = 53.2
    if long == None:
        long = -6.2

    #get coords
    x = long * (k * np.pi/180.0)
    y = np.log(np.tan((90 + lat) * np.pi/360.0)) * k
    df = wgs84_to_web_mercator(df, lon="LONG", lat="LAT")
    colormap = cmaps[colorby]
    df['color'] = [colormap[i] if i in colormap else 'gray' for i in df[colorby]]

    source = ColumnDataSource(df)
    #draw figure
    p = figure(x_range=(x-100000, x+100000), y_range=(y-100000, y+100000),
               x_axis_type="mercator", y_axis_type="mercator", tools=tools,
               plot_height=500, sizing_mode=sizing_mode)
    p.add_tile(tile_provider)
    p.circle(x='x', y='y', size=15, alpha=0.7, color='color', source=source)
    #data_table = DataTable(source=source, columns=columns, width=400, height=280)
    p.toolbar.logo = None

    infodiv = Div(width=400, height=p.plot_height, height_policy="fixed")
    infodiv.text = 'info'

    #callbacks
    #tap event - get mapr coords from plot and send to form so
    #we can refresh lat and long??
    tap_event = CustomJS(args=dict(s=source,plot=p,d=infodiv),code="""

        var inds = cb_obj.indices;
        var line = "<span style=%r><b>";
        var text = d.text.concat(line+'\\n');
        text = text.concat(String(cb_obj.x));
        d.text = text;

    """)
    from bokeh import events
    point_events = [ events.Tap, events.DoubleTap ]
    p.js_on_event(events.Tap, tap_event)

    source.selected.js_on_change('indices',
        CustomJS(args=dict(s=source,div=infodiv), code="""
        var inds = cb_obj.indices;
        var d = s.data;
        for (var i = 0; i < inds.length; i++) {
            console.log(i);
        }
        
    """))

    #dropdown
    menu = [('carttodbpositron','CARTODBPOSITRON')]
    dropdown = Dropdown(label="Tile", menu=menu)
    callback = CustomJS(args=dict(source=source,plot=p),code="""
        var r = plot.renderers;
        console.log(r);
        var b = cb_obj.value;
        var data = source.data;
        var x = data['x']
        source.change.emit();
    """)
    dropdown.js_on_event("menu_item_click", callback)

    names = list(zip(df.name,df.name))
    callback = CustomJS(args=dict(s=source,p=p,d=infodiv),code="""
        d.text = "hello";
    """)
    multiselect = MultiSelect(value=["s"], options=names)
    multiselect.js_on_change("value", callback)

    hover = p.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([
        ("name", "@name"),
        ("species", "@species"),
        ("spoligotype", "@spoligotype"),
        ("clade", "@clade"),
        ("nearest", "@nearest"),
    ])

    l = row(column(dropdown,multiselect,p),column(infodiv))
    script, div = components(l)
    return script, div

def draw_tree():
    """Draw tree based on sample distances"""


    return

def sample_tree(n=20):

    import toytree
    tre = toytree.rtree.coaltree(n)
    ## assign random edge lengths and supports to each node
    for node in tre.treenode.traverse():
        node.dist = np.random.exponential(1)
        node.support = int(np.random.uniform(50, 100))

    canvas,axes,mark = tre.draw(
                    width=400,
                    height=500,
                    scalebar=True, **tree_style)
    toyplot.html.render(canvas, "templates/tree.html")
    return

@webapp.route('/')
def index():
    """main index page"""

    from bokeh.resources import INLINE
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()

    colorby = request.args.get("colorby")
    tile = request.args.get("tile")
    lat = request.args.get('lat')
    long = request.args.get('long')

    if tile == None:
        tile= 'CARTODBPOSITRON'
    if colorby == None:
        colorby = 'species'

    msg = help_msg()
    df = pd.read_csv('wicklow_test.csv')
    form = ControlsForm(**request.args)

    script, div = bokeh_map(df, tile_provider=tile, colorby=colorby, lat=lat)

    sample_tree()
    return render_template("index.html",form=form,script=script,plotdiv=div,
            long=long,lat=lat,
            js_resources=js_resources, css_resources=css_resources, msg=msg)

@webapp.route('/links')
def links():
    return render_template("links.html")

def main():
    webapp.run(port=5000, debug=True)

if __name__ == '__main__':
	main()
