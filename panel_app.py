#!/usr/bin/env python3

"""
    BTBgenie prototype web app with Panel.
    Created April 2021
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
import sqlite3
from collections import OrderedDict
import toytree, toyplot
from snipgenie import snp_typing, trees, tools

from collections import OrderedDict
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, Slider
from bokeh.plotting import figure
from bokeh.themes import Theme
from bokeh.io import show, output_notebook
from bokeh.models import (DataTable, GeoJSONDataSource, ColumnDataSource, HoverTool, renderers,
                          Label, LabelSet, CustomJS, MultiSelect, Dropdown, Div)
from bokeh.tile_providers import CARTODBPOSITRON, get_provider

import panel as pn
import panel.widgets as pnw

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__))
DATABASE = '../notebooks/'
df = pd.read_csv('ireland_test_data.csv')
df = df.fillna('')

species_colors = {'Cow':'green', 'Badger':'blue', 'Deer':'red','Dog':'orange'}
sb_colors = trees.colors_from_labels(df,'name','spoligotype')
clade_colors = {'Monaghan-1':'brown','Monaghan-2':'orange','Monaghan-3':'red',
                'Wicklow-1':'green','NI-1':'blue','Unknown':'gray','':'gray'}
county_colors = trees.colors_from_labels(df,'name','county')
cmaps = {'species': species_colors,'spoligotype':sb_colors,'clade':clade_colors,'county':county_colors}
providers = ['CARTODBPOSITRON','STAMEN_TERRAIN','OSM','ESRI_IMAGERY']


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

template = """
{% extends base %}

<!-- goes in body -->
{% block postamble %}
<link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css">
{% endblock %}

<!-- goes in body -->
{% block contents %}
{{ app_title }}
<p>This is a Panel app with a custom template allowing us to compose multiple Panel objects into a single HTML document.</p>
<br>
<div class="container">
  <div class="row">
    <div class="col-sm">
      {{ embed(roots.A) }}
    </div>
    <div class="col-sm">
      {{ embed(roots.B) }}
    </div>
  </div>
</div>
{% endblock %}
"""

style1 = {'background':'#f5f5d6','padding':'5px','font':'monospace'}

def sample_tree(n=10):

    tre = toytree.rtree.coaltree(n)
    ## assign random edge lengths and supports to each node
    for node in tre.treenode.traverse():
        node.dist = np.random.exponential(1)
        node.support = int(np.random.uniform(50, 100))

    canvas,axes,mark = tre.draw(
                    width=350,
                    height=500,
                    scalebar=True, **tree_style)
    toyplot.html.render(canvas, "tree.html")
    return

def get_tree(df):
    """Get a tree from a selection of samples
       uses encoded snp data from dataframe/db to make a distance matrix
    """

    #decode snps
    snpmat = df.snps.apply(snp_typing.decode_snps)
    snpmat.index = df.name
    #print (snpmat[:4])
    tre = snp_typing.tree_from_snps(snpmat.T)
    return tre

def draw_tree(tre, df, colorby=None, layout='r', font_size=10, root=None):
    """draw the tree with given options"""

    if root not in ['',None]:
        tre = tre.root(root)
    tipnames = tre.get_tip_labels()
    node_colors = None
    node_sizes = None
    if colorby not in ['',None]:
        mapping = dict(zip(df.name,df[colorby]))
        colormap =  cmaps[colorby]
        tip_colors = [colormap[mapping[i]] if (i in mapping and i!='') else 'gray' for i in tipnames]
        if len(tipnames)>40:
            node_sizes=[0 if i else 6 for i in tre.get_node_values(None, 1, 0)]
            node_colors = [colormap[mapping[n]] if n in mapping else 'gray' for n in tre.get_node_values('name', True, True)]
            tipnames = ['' for i in tipnames]
    else:
        tip_colors = None
    tip_labels_style={
            "font-size": "%spx" %font_size,
            "-toyplot-anchor-shift": "13px",
        }
    #render to html
    canvas, axes, mark = tre.draw(tip_labels=tipnames, tip_labels_colors=tip_colors,tip_labels_style=tip_labels_style,
                                  layout=layout,node_colors=node_colors,node_sizes=node_sizes,
                                  scalebar=True, width=400, height=500);
    toyplot.html.render(canvas, "tree.html");
    return canvas

def wgs84_to_web_mercator(df, lon="LON", lat="LAT"):

      k = 6378137
      df["x"] = df[lon] * (k * np.pi/180.0)
      df["y"] = np.log(np.tan((90 + df[lat]) * np.pi/360.0)) * k
      return df

def test_map():
    tile_provider = get_provider(CARTODBPOSITRON)
    # range bounds supplied in web mercator coordinates
    p = figure(x_range=(-2000000, 6000000), y_range=(-1000000, 7000000),
               x_axis_type="mercator", y_axis_type="mercator")
    p.add_tile(tile_provider)
    return p

def bokeh_map(df=None, long=None, lat=None,
              tile_provider='CARTODBPOSITRON', colorby='species',
              labels=None):
    """Bokeh map"""

    tile_provider = get_provider(tile_provider)
    tools = "pan,wheel_zoom,box_zoom,hover,tap,lasso_select,reset,save"
    sizing_mode='stretch_both'

    # range bounds supplied in web mercator coordinates
    k = 6378137
    if lat == None:
        lat = 53.5
    if long == None:
        long = -7

    #get coords
    x = long * (k * np.pi/180.0)
    y = np.log(np.tan((90 + lat) * np.pi/360.0)) * k
    df = wgs84_to_web_mercator(df, lon="LONG", lat="LAT")
    colormap = cmaps[colorby]
    df['color'] = [colormap[i] if i in colormap else 'gray' for i in df[colorby]]
    df['label'] = ''
    source = ColumnDataSource(df)
    #draw figure
    p = figure(x_range=(x-200000, x+200000), y_range=(y-200000, y+200000),
               x_axis_type="mercator", y_axis_type="mercator", tools=tools,
               plot_width=500, plot_height=500, sizing_mode=sizing_mode)
    p.add_tile(tile_provider)
    p.circle(x='x', y='y', size=15, alpha=0.7, color='color', source=source)#, legend_group=colorby)

    labels = LabelSet(x='x', y='y', text='label',text_font_size='10pt',
                     x_offset=5, y_offset=5, source=source, render_mode='canvas')
    p.add_layout(labels)
    p.toolbar.logo = None
    p.match_aspect = True
    hover = p.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([
        ("name", "@name"),
        ("species", "@species"),
        ("spoligotype", "@spoligotype"),
        ("clade", "@clade"),
        ("nearest", "@nearest"),
        ("county", "@county"),
    ])
    return p

def map_dash():
    """Map dashboard"""

    names = sorted(list(df.name.unique()))
    tre = None
    cols = df.columns[:6]
    cats=['species','clade','spoligotype','county']
    labels=['','name','clade']
    map_pane = pn.pane.Bokeh(width=400)
    tree_pane = pn.pane.HTML(width=300)
    #tree_slider = pnw.FloatSlider(start=0,end=1,width=200)
    tree_layout_select = pnw.Select(name='tree layout',options=['r','c','d'],width=200)
    root_select = pnw.Select(name='root on',options=[''],width=200)
    tile_select = pnw.Select(name='tile layer',options=providers,width=200)
    colorby_select = pnw.Select(name='color by',options=cats,width=200)
    label_select = pnw.Select(name='label',options=labels,width=200)
    name_select = pnw.MultiSelect(name='name',options=names,size=8,width=200)
    btn = pnw.Button(name='Reset', button_type='primary',width=200)
    info_pane = pn.pane.HTML(style=style1,width=200, height=200,sizing_mode='stretch_both')
    df_pane = pn.pane.DataFrame(df[cols],width=500,height=200,sizing_mode='scale_both',max_rows=20,index=False)
    empty_pane = pn.pane.HTML(width=300,style=style1,sizing_mode='scale_height')
    empty_pane.object = 'lskdklasdlkjsad'
    loading = pn.indicators.LoadingSpinner(value=False, width=100, height=100)

    def update_info(attr,new,old):
        #print(new,old)
        info_pane.object = '<p>%s,%s</p>' %(int(new),int(old))

    def zoom_to_points(sel, p, pad=50000):
        #zoom with aspect conserved
        x1=sel.x.min()-pad
        x2=sel.x.max()+pad
        y1=sel.y.min()-pad
        y2=sel.y.max()+pad
        xr = x2-x1
        yr = y2-y1
        if yr<xr:
            y2 = y1+xr-pad*1.3
        else:
            x2 = x1+yr-pad*1.3

        p.x_range.update(start=x1,end=x2)
        p.y_range.update(start=y1,end=y2)
        return

    def items_selected(event):

        global tre
        items = name_select.value
        root_select.options = ['']+items
        colorby = colorby_select.value
        info_pane.object = '\n'.join(items)
        p = map_pane.object
        source = p.renderers[1].data_source
        sel = df[df.name.isin(items)]
        df_pane.object = sel[cols]
        #show these points only on map
        source.data = dict(sel)

        #zoom to points selected
        zoom_to_points(sel, p)

        #get a tree
        if len(sel)>=3:
            loading.value = True
            tre = get_tree(sel)
            draw_tree(tre, sel, colorby, layout=tree_layout_select.value)
            tree_pane.object = open('tree.html','r').read()
            loading.value = False
        else:
            tree_pane.object = ''
        return

    def points_selected(attr,new,old):
        """bokeh callback for lasso"""

        global tre
        colorby = colorby_select.value
        ind =[int(n) for n in new]
        sel = df.loc[ind]
        df_pane.object = sel[cols]
        if len(sel)>=3:
            tre = get_tree(sel)
            draw_tree(tre, sel, colorby, layout=tree_layout_select.value)
            tree_pane.object = open('tree.html','r').read()
        return

    def tap_callback(event):
        #tap tool callback

        p = map_pane.object
        source = p.renderers[1].data_source
        ind = source.selected.indices
        sel = df.loc[ind]
        df_pane.object = sel[cols]
        #get nearest
        if len(sel)>0:
            s = sel.iloc[0].nearest
            near = s.split()
            info_pane.object = s
        return

    def draw_map(event):

        p = map_pane.object = bokeh_map(df)
        p.x_range.on_change('start', update_info)
        source = p.renderers[1].data_source
        source.selected.on_change('indices', points_selected)
        p.on_event('tap',tap_callback)
        tree_pane.object = ''
        return

    def update_map(event):

        p = map_pane.object
        source = p.renderers[1].data_source
        p.renderers = [x for x in p.renderers if not str(x).startswith('TileRenderer')]
        rend = renderers.TileRenderer(tile_source= get_provider(tile_select.value))
        p.renderers.insert(0, rend)
        colorby = colorby_select.value
        colormap = cmaps[colorby]
        df['color'] = [colormap[i] if i in colormap else 'gray' for i in df[colorby]]

        info_pane.object = '<p>%s,%s</p>' %(p.x_range.start,p.x_range.end)
        if label_select.value != '':
            df['label'] = df[label_select.value]
        else:
            df['label'] = ''
        source.data = dict(df)
        update_tree(event)
        return

    def update_tree(event):
        global tre
        #use subset of samples if selected
        items = name_select.value
        if tre != None:
            sel = df[df.name.isin(items)]
            draw_tree(tre, sel, colorby_select.value, layout=tree_layout_select.value, root=root_select.value)
            tree_pane.object = open('tree.html','r').read()

    draw_map(None)
    btn.on_click(draw_map)
    #label_box = pnw.Checkbox(name='Show labels')
    tile_select.param.watch(update_map,'value')
    colorby_select.param.watch(update_map,'value')
    label_select.param.watch(update_map,'value')
    name_select.param.watch(items_selected,'value')
    tree_layout_select.param.watch(update_tree,'value')
    root_select.param.watch(update_tree,'value')

    #layout dashboard
    app = pn.Column(pn.Row(pn.Column(tile_select,colorby_select,label_select,name_select,tree_layout_select,root_select,btn,
                                     background='whitesmoke',sizing_mode='stretch_height'),
                           pn.Column(map_pane,width=600),tree_pane,loading),pn.Row(pn.Row(info_pane,scroll=True),df_pane,height=200,scroll=True))
    return app

bootstrap = pn.template.BootstrapTemplate(title='BTBGenie Sample App',
            favicon='static/logo.png',logo='static/logo.png',header_color='blue')
pn.config.sizing_mode = 'stretch_width'
app = map_dash()
bootstrap.main.append(app)
bootstrap.servable()

if __name__ == '__main__':
	pn.serve(bootstrap, port=5000)
