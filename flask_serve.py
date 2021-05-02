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
    name = SelectField('name', choices=[])
    path = TextField('path', default='results')
    n = TextField('n', default='2')
    #kinds = [(i,i) for i in plotkinds]
    #kind = SelectField('plot kind', choices=kinds)
    #submit = SubmitField()

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

def bokeh_map(df=None):
    """Bokeh map"""

    from bokeh.plotting import figure, show
    from bokeh.embed import components
    from bokeh.tile_providers import CARTODBPOSITRON,STAMEN_TERRAIN, get_provider
    from bokeh.models import DataTable, GeoJSONDataSource, ColumnDataSource, HoverTool, CustomJS, Dropdown
    from bokeh.layouts import column, row
    import json
    from collections import OrderedDict
    tile_provider = get_provider(CARTODBPOSITRON)
    tile_provider2 = get_provider(STAMEN_TERRAIN)
    tools = "pan,wheel_zoom,hover,tap,reset,save"
    sizing_mode='stretch_width'

    # range bounds supplied in web mercator coordinates
    k = 6378137
    long = -6.2
    lat = 53.5
    x = long * (k * np.pi/180.0)
    y = np.log(np.tan((90 + lat) * np.pi/360.0)) * k

    p = figure(x_range=(x-100000, x+100000), y_range=(y-100000, y+100000),
               x_axis_type="mercator", y_axis_type="mercator", tools=tools,
               plot_height=500, sizing_mode=sizing_mode)
    p.add_tile(tile_provider)
    #p.add_tile(tile_provider2)

    #geo_source = GeoJSONDataSource(df.to_json())
    df = wgs84_to_web_mercator(df, lon="LONG", lat="LAT")
    df['color'] = [species_colors[i] if i in species_colors else 'gray' for i in df.species]
    #print (df.iloc[4])
    source = ColumnDataSource(df)
    p.circle(x='x', y='y', size=15, alpha=0.7, color='color', source=source)
    #data_table = DataTable(source=source, columns=columns, width=400, height=280)

    menu = [("Item 1", "CARTODBPOSITRON"), ("Item 2", "STAMEN_TERRAIN")]

    dropdown = Dropdown(label="Layer", menu=menu)
    dropdown.js_on_event("menu_item_click",
            CustomJS(code="console.log('dropdown: ' + this.item, this.toString())"))

    hover = p.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([
        ("name", "@name"),
        ("species", "@species"),
        ("spo type", "@SB"),
        ("clade", "@clade"),
        ("nearest", "@nearest"),
    ])
    p.toolbar.logo = None
    l = column(p, dropdown)
    script, div = components(l)
    return script, div

def get_tree():
    """show a tree"""

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

    path = request.args.get("path")
    msg = help_msg()
    df = pd.read_csv('wicklow_test.csv')
    form = ControlsForm()
    #map = base_map()
    script, div = bokeh_map(df)

    #show_dataframe(df, map)
    #map = folium_map._repr_html_()
    #map.save('templates/map.html')
    sample_tree()
    return render_template("index.html",form=form,script=script,plotdiv=div,
            js_resources=js_resources, css_resources=css_resources, msg=msg)

@webapp.route('/links')
def links():
    return render_template("links.html")

def main():
    webapp.run(port=5000, debug=True)

if __name__ == '__main__':
	main()
