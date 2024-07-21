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

    view_history = {}
    def update(event=None, sub=None):
        """Update selection"""

        provider = provider_input.value
        pcmap = cmap = cmap_input.value
        col = colorby_input.value
        ms = markersize_input.value
        legend = legendbtn.value

        global selected

        if sub is None:
            sub = selected
        else:
            selected = sub
        #if len(sub)>2000:
        #    return
        #filters
        sub = apply_filters(sub)
        if len(sub[col].unique())>20:
            cmap=None
        sub['color'],c = tools.get_color_mapping(sub, col, cmap)
        sub['marker'] = sub.Species.map({'Bovine':'circle','Badger':'square'})

        sp = parcels[parcels.SPH_HERD_N.isin(sub.HERD_NO)].copy()
        if len(sp.SPH_HERD_N.unique())>20:
            pcmap=None
        sp['color'],c = tools.get_color_mapping(sp, 'SPH_HERD_N',pcmap)

        mov = gui.get_moves_bytag(sub, moves, lpis_cent)

        if parcelsbtn.value == True:
            pcls=sp
        else:
            pcls=None
        p = bokeh_plot.plot_selection(sub, pcls, provider=provider, ms=ms, col=col, legend=legend)

        if movesbtn.value == True:
            bokeh_plot.plot_moves(p, mov, lpis_cent)
            #moves_pane.value = mov
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
        info_pane.object = f'{len(sub)} samples'
        if treebtn.value ==True:
            update_tree(sub=sub)
        return

    def apply_filters(df):
        minsize = clustersizeslider.value
        key = groupby_input.value
        groups = df[key].value_counts()
        groups = groups[groups>minsize]
        df = df[df[key].isin(groups.index)].copy()
        #print (df)
        start, end = timeslider.value
        df = df[(df.Year>=start) & (df.Year<=end)].copy()
        if homebredbox.value == True:
            df = df[df.Homebred=='yes'].copy()
        return df

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
        f = bokeh_plot.split_view(selected, col, parcels, provider)
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

    def update_tree(event=None, sub=None):

        treefile = get_tree(sub.index)
        col = colorby_input.value
        #html = draw_tree(treefile, selected, col)
        html = draw_toytree(treefile, sub, col, cmap=cmap_input.value,
                            tiplabelcol=tiplabel_input.value, markercol='Species', height=400)
        tree_pane.object = html
        return

    def set_stop(event):
        global stop
        stop = True
        print ('STOP')

    #p = bokeh_plot.test()
    plot_pane = pn.pane.Bokeh()
    split_pane = pn.pane.Bokeh()
    timeline_pane = pn.pane.Bokeh()
    info_pane = pn.pane.Str('', styles={'font-size': '10pt'})
    #log_pane = pnw.TextAreaInput(disabled=True,height=600,width=400)
    meta_pane = pnw.Tabulator(disabled=True,height=600,page_size=50,
                              header_filters=True,frozen_columns=['sample'])
    selected_pane = pnw.Tabulator(show_index=False,disabled=True,
                                  frozen_columns=['sample'],height=600)
    moves_pane = pnw.Tabulator(show_index=False,disabled=True,
                                  frozen_columns=['tag'],height=600)
    tree_pane = pn.pane.HTML(height=400)
    network_pane = pn.pane.Bokeh()
    snpdist_pane = pn.pane.HTML(height=400)

    w=140
    progress = pn.indicators.Progress(name='Progress', value=0, width=600, bar_color='primary')
    cols = [None]+gui.get_ordinal_columns(meta)
    #tables = pn.Column(pn.Tabs(('metadata',meta_pane),('selected',selected_pane)))

    groupby_input = pnw.Select(name='group by',options=cols,value='snp7',width=w)
    groups_table = pnw.Tabulator(disabled=True, widths={'index': 70}, layout='fit_data',pagination=None, height=250, width=w)
    colorby_input = pnw.Select(name='color by',options=cols,value='snp7',width=w)
    cmap_input = pnw.Select(name='colormap',options=gui.colormaps,value='Set1',width=w)
    provider_input = pnw.Select(name='provider',options=['']+bokeh_plot.providers,value='CartoDB Positron',width=w)
    markersize_input = pnw.FloatInput(name='marker size', value=10, step=1, start=2, end=100,width=w)
    tiplabel_input = pnw.Select(name='tip label',options=list(meta.columns),value='sample',width=w)
    #showborders_input = pnw.Checkbox(name='Show counties', value=False)
    widgets = pn.Column(pn.WidgetBox(groupby_input,groups_table,colorby_input,cmap_input,tiplabel_input,
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

    timeslider = pnw.IntRangeSlider(name='Time',width=150,
                    start=2000, end=2024, value=(2000, 2024), step=1)
    clustersizeslider = pnw.IntSlider(name='Min. Cluster Size',width=150,
                    start=1, end=10, value=1, step=1)
    homebredbox = pnw.Checkbox(name='Homebred',value=False)
    filters = pn.Row(timeslider,clustersizeslider,homebredbox,info_pane)

    groupby_input.param.watch(update_groups, 'value')
    groups_table.param.watch(select_group, 'selection')
    provider_input.param.watch(set_provider, 'value')
    colorby_input.param.watch(update, 'value')
    cmap_input.param.watch(update, 'value')
    tiplabel_input.param.watch(update, 'value')
    treebtn.param.watch(update, 'value')
    parcelsbtn.param.watch(update, 'value')
    movesbtn.param.watch(update, 'value')
    legendbtn.param.watch(update, 'value')
    timeslider.param.watch(update, 'value')
    clustersizeslider.param.watch(update, 'value')
    homebredbox.param.watch(update, 'value')

    #pn.bind(update_tree, treebtn, watch=True)
    pn.bind(split_view, splitbtn, watch=True)
    pn.bind(select_region, selectregionbtn, watch=True)
    #pn.bind(filter_selection, timeslider, watch=True)

    #cat scatter plot
    scatter_pane = pn.pane.Bokeh(height=400)
    row_input = pnw.Select(name='row',options=cols,width=w)
    col_input = pnw.Select(name='col',options=cols,width=w)
    scolor_input = pnw.Select(name='color',options=cols,width=w)
    analysis_pane = pn.Column(scatter_pane,pn.Row(row_input,col_input,scolor_input))
    def update_scatter(event=None):
        global selected
        p = bokeh_plot.cat_plot(selected, row_input.value, col_input.value, scolor_input.value)
        scatter_pane.object = p
    row_input.param.watch(update_scatter, 'value')
    col_input.param.watch(update_scatter, 'value')
    scolor_input.param.watch(update_scatter, 'value')

    app = pn.Column(
                pn.Row(widgets,
                       toolbar,
                pn.Column(
                    pn.Tabs(('map view',pn.Column(plot_pane,filters)),
                            ('split view',pn.Column(split_pane)),
                            ('main table', meta_pane),
                            ('selected', selected_pane),
                            ('moves', moves_pane),
                             sizing_mode='stretch_both'),
                          ),
                pn.Column(pn.Tabs(('tree',tree_pane),('network',network_pane),
                                  ('snpdist',snpdist_pane),('scatter',analysis_pane)),timeline_pane, width=500))
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
        print ('please provide a project file')
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
