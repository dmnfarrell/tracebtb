#!/usr/bin/env python3

"""
    TraceBTB panel dashboard
    Created June 2024
    Copyright (C) Damien Farrell
"""

import sys,os,time,re
import pickle
import glob,io
import json
import math
import pylab as plt
import pandas as pd

from bokeh.io import show
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, GeoJSONDataSource, GMapOptions, GMapPlot, TileSource, HoverTool, BoxZoomTool
from bokeh.models.glyphs import Patches, Circle
from bokeh.layouts import layout, column
from shapely.geometry import Polygon, Point
from bokeh.tile_providers import get_provider, Vendors
from bokeh.models import Arrow, NormalHead, OpenHead, VeeHead

import panel as pn
import panel.widgets as pnw
pn.extension('tabulator', css_files=[pn.io.resources.CSS_URLS['font-awesome']])

from tracebtb import gui, tools, plotting, trees, bokeh_plot

def get_icon(name):
    filename = f'../tracebtb/icons/{name}.svg'
    with open(filename, "r") as f:
        data = f.read()
    svg_tag = re.search(r'(<svg[^>]*>.*?</svg>)', data, re.DOTALL).group(1)
    return svg_tag

def dashboard(gdf, parcels, moves=None):

    def update(event=None, sub=None):
        """Update selection"""

        provider = provider_input.value
        cmap = cmap_input.value
        col = colorby_input.value

        global selected

        if sub is None:
            sub = selected
        else:
            selected = sub
        if len(sub)>150:
            return

        mov = gui.get_moves_bytag(sub, moves, lpis_cent)

        sub['color'],c = tools.get_color_mapping(sub, col, cmap)
        sub['marker'] = sub.Species.map({'Bovine':'circle','Badger':'square'})
        if parcelsbtn.value == True:
            sp = parcels[parcels.SPH_HERD_N.isin(sub.HERD_NO)].copy()
            sp['color'],c = tools.get_color_mapping(sp, 'SPH_HERD_N','Paired')
        else:
            sp = None
        p = bokeh_plot.plot_selection(sub, sp, provider=provider)
        if movesbtn.value == True:
            bokeh_plot.plot_moves(p, mov, lpis_cent)
        plot_pane.object = p
        plot_pane.param.trigger('object')
        #change selection table
        #selected_pane.value = sub
        return

    def show_moves(event=None):

        return

    def update_groups(event=None):
        """change group choices"""

        groupby = groupby_input.value
        vals = pd.DataFrame(meta[groupby].value_counts())#.reset_index()
        groups_table.value = vals
        return

    def select_group(event=None):

        df = groups_table.value
        rows = groups_table.selection
        groups = list(df.iloc[rows].index)
        #print (groups)
        key = groupby_input.value
        sub = meta[meta[key].isin(groups)].copy()
        update(sub=sub)
        return

    def export_selection(event):
        """export selection to html file"""

        return

    def split_view(event=None):

        provider = provider_input.value
        cmap = cmap_input.value
        col = colorby_input.value
        global selected
        f = bokeh_plot.split_view(selected, 'snp3', parcels, provider)
        print (f)
        split_pane.object = f
        split_pane.param.trigger('object')
        return

    def get_tree(idx):
        #replace this with just selecting tips on main tree?
        treefile = 'tree.newick'
        M = snpdist.loc[idx,idx]
        trees.tree_from_distmatrix(M, treefile)
        return treefile

    def update_tree(event):
        global selected
        treefile = get_tree(selected.index)
        col = colorby_input.value
        html = draw_tree(treefile, selected, col)
        #html = '<h1>hello</h1>'
        tree_pane.object = html

    def set_stop(event):
        global stop
        stop = True
        print ('STOP')

    #plot_pane = pn.pane.Matplotlib(plt.Figure(),height=300)
    #p = bokeh_plot.test()
    plot_pane = pn.pane.Bokeh()
    split_pane = pn.pane.Bokeh()
    #log_pane = pnw.TextAreaInput(disabled=True,height=600,width=400)
    meta_pane = pnw.Tabulator(disabled=True,height=600,theme='modern',
                              header_filters=True,frozen_columns=['sample'])
    selected_pane = pnw.Tabulator(show_index=False,disabled=True,width=200,height=600)
    tree_pane = pn.pane.HTML()

    w=140
    progress = pn.indicators.Progress(name='Progress', value=0, width=600, bar_color='primary')
    cols = list(gdf.columns)
    #tables = pn.Column(pn.Tabs(('metadata',meta_pane),('selected',selected_pane)))

    groupby_input = pnw.Select(name='group by',options=cols,value='snp7',width=w)
    groups_table = pnw.Tabulator(disabled=True, widths={'index': 70}, layout='fit_data',pagination=None, height=250, width=w)
    colorby_input = pnw.Select(name='color by',options=cols,width=w)
    cmap_input = pnw.Select(name='colormap',options=gui.colormaps,width=w)
    provider_input = pnw.Select(name='provider',options=['']+bokeh_plot.providers,width=w)
    widgets = pn.Column(pn.WidgetBox(groupby_input,groups_table,colorby_input,cmap_input,provider_input), width=w+30)

    splitbtn = pnw.Button(icon=get_icon('plot-grid'), description='split view', icon_size='1.8em')
    treebtn = pnw.Button(icon=get_icon('tree'), icon_size='1.8em')
    parcelsbtn = pnw.Toggle(icon=get_icon('parcels'), icon_size='1.8em')
    movesbtn = pnw.Toggle(icon=get_icon('plot-moves'), icon_size='1.8em')
    toolbar = pn.Column(pn.WidgetBox(splitbtn,treebtn,parcelsbtn,movesbtn),width=70)

    groupby_input.param.watch(update_groups, 'value')
    groups_table.param.watch(select_group, 'selection')
    provider_input.param.watch(update, 'value')
    colorby_input.param.watch(update, 'value')
    cmap_input.param.watch(update, 'value')
    parcelsbtn.param.watch(update, 'value')
    movesbtn.param.watch(update, 'value')

    pn.bind(update_tree, treebtn, watch=True)
    pn.bind(split_view, splitbtn, watch=True)

    def execute(event):
        #run the model with widget
        pass

    app = pn.Column(
                pn.Row(widgets,
                       toolbar,
                pn.Column(
                    pn.Tabs(('map view',pn.Column(plot_pane, height=800)),
                            ('split view',pn.Column(split_pane, height=800)),
                                  ('main table', meta_pane),
                                  ('selected', selected_pane),
                             sizing_mode='stretch_both'),
                          ),
                pn.Column(tree_pane, width=300)),
        #sizing_mode='stretch_both', max_width=1300,
    )
    #app.sizing_mode='stretch_both'
    meta_pane.value = gdf
    update_groups()
    selected = meta[meta.snp7.isin(['24'])].copy()
    update(sub=selected)
    return app


bootstrap = pn.template.BootstrapTemplate(title='TraceBTB Web')
            #favicon='static/logo.png',logo='static/logo.png',header_color='blue')
pn.config.sizing_mode = 'stretch_width'

data = pickle.load(open('../test.tracebtb','rb'))
meta = data['meta']
moves = data['moves']
lpis_cent = data['lpis_cent']
parcels = data['parcels']
snpdist = data['snpdist']

app = dashboard(meta, parcels, moves)
bootstrap.main.append(app)
bootstrap.servable()

if __name__ == '__main__':
    pn.serve(bootstrap, port=5010)
