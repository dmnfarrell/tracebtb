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

module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
data_path = os.path.join(module_path,'data')
logoimg = os.path.join(module_path, 'logo.png')
iconpath = os.path.join(module_path, 'icons')

def get_icon(name):
    """Get svg icon"""

    filename = os.path.join(iconpath, f'{name}.svg')
    with open(filename, "r") as f:
        data = f.read()
    svg_tag = re.search(r'(<svg[^>]*>.*?</svg>)', data, re.DOTALL).group(1)
    return svg_tag

def draw_toytree(treefile, df, col, **kwargs):
    import toyplot
    from tracebtb import trees
    canvas = trees.draw_tree(treefile, df, col=col, **kwargs)
    toyplot.html.render(canvas, "temp.html")
    with open('temp.html', 'r') as f:
        html = f.read()
    return html

def draw_tree(treefile, df, col):
    """Draw newick tree with phylocanvas"""

    with open(treefile, 'r') as file:
        newick_data = file.read()
    #color_mapping = df[col].to_dict()
    #color_json = json.dumps(color_mapping)

    html = """

        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Phylocanvas Tree Example</title>
            <script src="https://unpkg.com/@phylocanvas/phylocanvas.gl@latest/dist/bundle.min.js"></script>

        </head>
        <body>
            <h1> TEST </h1>
            <div id="demo" style="border: 1px solid lightgray"></div>

            <script>
            const phylotree = new phylocanvas.PhylocanvasGL(
            document.querySelector("#demo"),
            {{
                showLabels: true,
                showLeafLabels: true,
                size: {{ width: 400, height: 500 }},
                source: `{n}`,
            }},
            );
            </script>
        </body>
        </html>
    """.format(n=newick_data)
    return html

def dashboard(meta, parcels, moves=None, lpis_cent=None, snpdist=None):
    """
    Dashboard app with panel for tracebtb.
    Args:
        meta: geodataframe of meta with locations of samples
        parcels: parcels for relevant samples from lpis
        moves: dataframe of moves
        lpis_cent: centroids from LPIS
    """

    def update(event=None, sub=None):
        """Update selection"""

        provider = provider_input.value
        cmap = cmap_input.value
        col = colorby_input.value
        ms = markersize_input.value

        global selected

        if sub is None:
            sub = selected
        else:
            selected = sub
        if len(sub)>2000:
            return

        mov = gui.get_moves_bytag(sub, moves, lpis_cent)

        sub['color'],c = tools.get_color_mapping(sub, col, cmap)
        sub['marker'] = sub.Species.map({'Bovine':'circle','Badger':'square'})

        sp = parcels[parcels.SPH_HERD_N.isin(sub.HERD_NO)].copy()
        sp['color'],c = tools.get_color_mapping(sp, 'SPH_HERD_N','Paired')

        if parcelsbtn.value == True:
            pcls=sp
        else:
            pcls=None
        p = bokeh_plot.plot_selection(sub, pcls, provider=provider, ms=ms)
        #if borders == True:
        #    bokeh_plot.plot_counties(p)
        if movesbtn.value == True:
            bokeh_plot.plot_moves(p, mov, lpis_cent)
            moves_pane.value = mov
            #moves_pane.param.trigger('value')
        plot_pane.object = p
        plot_pane.param.trigger('object')
        #change selection table
        #selected_pane.value = sub
        #timeline
        herdcolors = dict(zip(sp.SPH_HERD_N,sp.color))
        p = bokeh_plot.plot_moves_timeline(mov, herdcolors)
        timeline_pane.object = p
        timeline_pane.param.trigger('object')

        if treebtn.value ==True:
            update_tree()
        return

    def set_provider(event=None):
        p = plot_pane.object
        #remove old tile
        p.renderers = [x for x in p.renderers if not str(x).startswith('TileRenderer')]
        provider = provider_input.value
        p.add_tile(provider, retina=True)

    def update_groups(event=None):
        """change group choices"""

        groupby = groupby_input.value
        vals = pd.DataFrame(meta[groupby].value_counts())#.reset_index()
        groups_table.value = vals
        return

    def select_group(event=None):
        """Select groups from table"""

        df = groups_table.value
        rows = groups_table.selection
        groups = list(df.iloc[rows].index)
        #print (groups)
        key = groupby_input.value
        sub = meta[meta[key].isin(groups)].copy()
        update(sub=sub)
        return

    def select_region(event=None):
        """Select samples in region"""

        p = plot_pane.object
        xmin, xmax = p.x_range.start, p.x_range.end
        ymin, ymax = p.y_range.start, p.y_range.end
        sub = meta.cx[xmin:xmax, ymin:ymax]
        update(sub=sub)
        return

    def export_selection(event):
        """export selection to html file"""

        return

    def split_view(event=None):
        """split view"""

        provider = provider_input.value
        cmap = cmap_input.value
        col = colorby_input.value
        global selected
        f = bokeh_plot.split_view(selected, 'snp3', parcels, provider)
        #print (f)
        split_pane.object = f
        split_pane.param.trigger('object')
        return

    def get_tree(idx):
        #replace this with just selecting tips on main tree?
        treefile = 'tree.newick'
        M = snpdist.loc[idx,idx]
        trees.tree_from_distmatrix(M, treefile)
        return treefile

    def update_tree(event=None):
        global selected
        treefile = get_tree(selected.index)
        col = colorby_input.value
        cmap = cmap_input.value
        #html = draw_tree(treefile, selected, col)
        html = draw_toytree(treefile,  selected, col, cmap=cmap, tiplabelcol='sample', height=400)
        tree_pane.object = html
        return

    def set_stop(event):
        global stop
        stop = True
        print ('STOP')

    cmaps = ['Blues','Reds','Greens','RdBu','coolwarm','summer','winter','icefire','hot','viridis']
    #plot_pane = pn.pane.Matplotlib(plt.Figure(),height=300)
    #p = bokeh_plot.test()
    plot_pane = pn.pane.Bokeh()
    split_pane = pn.pane.Bokeh()
    timeline_pane = pn.pane.Bokeh()
    #log_pane = pnw.TextAreaInput(disabled=True,height=600,width=400)
    meta_pane = pnw.Tabulator(disabled=True,height=600,#theme='modern',
                              header_filters=True,frozen_columns=['sample'])
    selected_pane = pnw.Tabulator(show_index=False,disabled=True,
                                  frozen_columns=['sample'],height=600)
    moves_pane = pnw.Tabulator(show_index=False,disabled=True,
                                  frozen_columns=['tag'],height=600)
    tree_pane = pn.pane.HTML(height=400)

    w=140
    progress = pn.indicators.Progress(name='Progress', value=0, width=600, bar_color='primary')
    cols = gui.get_ordinal_columns(meta)
    #tables = pn.Column(pn.Tabs(('metadata',meta_pane),('selected',selected_pane)))

    groupby_input = pnw.Select(name='group by',options=cols,value='snp7',width=w)
    groups_table = pnw.Tabulator(disabled=True, widths={'index': 70}, layout='fit_data',pagination=None, height=250, width=w)
    colorby_input = pnw.Select(name='color by',options=cols,value='snp7',width=w)
    cmap_input = pnw.Select(name='colormap',options=gui.colormaps,width=w)
    provider_input = pnw.Select(name='provider',options=['']+bokeh_plot.providers,width=w)
    markersize_input = pnw.FloatInput(name='marker size', value=10, step=1, start=2, end=100,width=w)
    #showborders_input = pnw.Checkbox(name='Show counties', value=False)
    widgets = pn.Column(pn.WidgetBox(groupby_input,groups_table,colorby_input,cmap_input,
                                     provider_input,markersize_input), width=w+30)

    splitbtn = pnw.Button(icon=get_icon('plot-grid'), description='split view', icon_size='1.8em')
    selectregionbtn = pnw.Button(icon=get_icon('plot-region'), description='select in region', icon_size='1.8em')
    treebtn = pnw.Toggle(icon=get_icon('tree'), icon_size='1.8em')
    parcelsbtn = pnw.Toggle(icon=get_icon('parcels'), icon_size='1.8em')
    movesbtn = pnw.Toggle(icon=get_icon('plot-moves'), icon_size='1.8em')
    legendbtn = pnw.Toggle(icon=get_icon('legend'), icon_size='1.8em')
    reportbtn = pnw.Button(icon=get_icon('report'), description='report', icon_size='1.8em')
    toolbar = pn.Column(pn.WidgetBox(splitbtn,selectregionbtn,treebtn,
                                     parcelsbtn,movesbtn,legendbtn,reportbtn),width=70)

    groupby_input.param.watch(update_groups, 'value')
    groups_table.param.watch(select_group, 'selection')
    provider_input.param.watch(set_provider, 'value')
    colorby_input.param.watch(update, 'value')
    cmap_input.param.watch(update, 'value')
    treebtn.param.watch(update, 'value')
    parcelsbtn.param.watch(update, 'value')
    movesbtn.param.watch(update, 'value')

    #pn.bind(update_tree, treebtn, watch=True)
    pn.bind(split_view, splitbtn, watch=True)
    pn.bind(select_region, selectregionbtn, watch=True)

    def execute(event):
        #run the model with widget
        pass

    app = pn.Column(
                pn.Row(widgets,
                       toolbar,
                pn.Column(
                    pn.Tabs(('map view',pn.Column(plot_pane)),
                            ('split view',pn.Column(split_pane)),
                              ('main table', meta_pane),
                              ('selected', selected_pane),
                              ('moves', moves_pane),
                             sizing_mode='stretch_both'),
                          ),
                pn.Column(tree_pane,timeline_pane, width=500))
    )
    app.sizing_mode='stretch_both'
    meta_pane.value = meta
    update_groups()
    selected = meta[meta.snp7.isin(['24'])].copy()
    update(sub=selected)
    return app

def main():
    "Run the application"
    from argparse import ArgumentParser
    parser = ArgumentParser(description='TracebTB')
    parser.add_argument("-p", "--proj", dest="project",default=None,
                            help="load project file", metavar="FILE")
    args = parser.parse_args()
    #load data
    if args.project == None:
        print ('provide a project file')
        exit()

    data = pickle.load(open(args.project,'rb'))
    meta = data['meta']
    moves = data['moves']
    lpis_cent = data['lpis_cent']
    parcels = data['parcels']
    snpdist = data['snpdist']

    #create template
    bootstrap = pn.template.BootstrapTemplate(title='TracebTB Web',
                        favicon=logoimg,logo=logoimg,header_color='green')
    #pn.config.sizing_mode = 'stretch_width'
    app = dashboard(meta, parcels, moves, lpis_cent, snpdist)
    bootstrap.main.append(app)
    bootstrap.servable()
    pn.serve(bootstrap, port=5010,
             basic_auth={'guest':"mbovis"}, cookie_secret='cookie_secret')

if __name__ == '__main__':
    main()
