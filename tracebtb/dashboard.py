#!/usr/bin/env python3

"""
    TraceBTB panel dashboard
    Created June 2024
    Copyright (C) Damien Farrell
"""

import sys,os,time,re
import platform
from datetime import datetime
import pickle
import glob,io
import json
import math
import pylab as plt
import pandas as pd
import geopandas as gpd

from bokeh.io import show
from bokeh.plotting import figure
from bokeh.models import Range1d, CustomJS, TapTool
import panel as pn
import panel.widgets as pnw

from tracebtb import tools, plotting, trees, bokeh_plot

module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
data_path = os.path.join(module_path,'data')
logoimg = os.path.join(module_path, 'logo.png')
iconpath = os.path.join(module_path, 'icons')
home = os.path.expanduser("~")
if platform.system() == 'Windows':
    configpath = os.path.join(os.environ['APPDATA'], 'tracebtb')
else:
    configpath = os.path.join(home, '.config','tracebtb')
if not os.path.exists(configpath):
    os.makedirs(configpath, exist_ok=True)
configfile = os.path.join(configpath, 'settings.json')
report_file = 'report.html'
selections_file = os.path.join(configpath,'selections.json')
layers_file = os.path.join(configpath,'layers.gpkg')

speciescolors = {'Bovine':'blue','Badger':'red','Ovine':'green'}
card_style = {
    'background': '#f9f9f9',
    'border': '1px solid #bcbcbc',
    'padding': '5px',
    'margin': '5px',
    'box-shadow': '4px 4px 4px #bcbcbc'
}
colormaps = ['Paired', 'Dark2', 'Set1', 'Set2', 'Set3',
            'tab10', 'tab20', 'tab20b', 'tab20c']
stylesheet = """
.tabulator-cell {
    font-size: 12px;
}
.tabulator-col-title {
    font-size: 11px;
}
"""
defaults = {'dashboard':{'lpis_master_file':''}}

def get_icon(name):
    """Get svg icon"""

    filename = os.path.join(iconpath, f'{name}.svg')
    with open(filename, "r") as f:
        data = f.read()
    svg_tag = re.search(r'(<svg[^>]*>.*?</svg>)', data, re.DOTALL).group(1)
    return svg_tag

def get_tree(snpdist, idx):
    """Get tree from snpdist"""

    treefile = 'tree.newick'
    M = snpdist.loc[idx,idx]
    trees.tree_from_distmatrix(M, treefile)
    return treefile

def keep_tips(tree, tips):
    """Prune  tips not in a list"""

    import copy
    new = copy.deepcopy(tree)
    for tip in new.get_terminals():
        if tip.name not in tips:
            new.prune(tip)
    return new

def draw_toytree(treefile, df, **kwargs):
    """Toytree plot"""

    import toyplot
    from tracebtb import trees
    with open(treefile, 'r') as f:
        s = f.read()
    if len(df)<1800:
        canvas = trees.draw_tree(s, df, **kwargs)
        toyplot.html.render(canvas, "temp.html")
    else:
        return '<h2>too many tips</b>'
    with open('temp.html', 'r') as f:
        html = f.read()
    return html

def get_figure_coords(p):
    """Coords of plot figure"""

    xmin, xmax = p.x_range.start, p.x_range.end
    ymin, ymax = p.y_range.start, p.y_range.end
    return xmin, xmax, ymin, ymax

def find_neighbours(gdf, dist, lpis_cent, lpis):
    """Find neighbours"""

    found = []
    for x in gdf.geometry:
        dists = lpis_cent.distance(x)
        points = lpis_cent[(dists<=dist) & (dists>10)]
        found.append(points)

    found = pd.concat(found).drop_duplicates()
    x = lpis[lpis.SPH_HERD_N.isin(found.SPH_HERD_N)]
    #exclude those in source gdf
    x = x[~x.SPH_HERD_N.isin(gdf.HERD_NO)]
    return x

def shared_borders(parcels, lpis):
    """Find herds with shared borders"""

    found = []
    for i, r in parcels.iterrows():
        #farms with shared borders
        #polygon = lpis[lpis.SPH_HERD_N==r.HERD_NO].geometry
        polygon = r.geometry
        x = lpis[lpis.touches(polygon)]
        found.append(x)
    if len(found) == 0:
        return
    found = pd.concat(found)
    found = found[~found.SPH_HERD_N.isin(parcels.SPH_HERD_N)]
    return found

def plot_overview(df):
    """Plot overview map"""

    fig,ax=plt.subplots(1,1)
    bokeh_plot.counties_gdf.plot(color='#E7F7C3', edgecolor='gray',
                        lw=1,alpha=0.7,
                        ax=ax)
    if len(df) != 0:
        df.to_crs('EPSG:3857').plot(ax=ax,color=df.color)
    ax.axis('off')
    plt.tight_layout()
    return fig

def report(sub, parcels, moves, col, lpis_cent, snpdist, cmap='Set1'):
    """
    Save a report with custom panes.
    """

    from datetime import date
    styles={ "margin-top": "10px", "margin-bottom": "20px"}
    CARD_STYLE = {
      "box-shadow": "rgba(50, 50, 93, 0.25) 0px 6px 12px -2px, rgba(0, 0, 0, 0.3) 0px 3px 7px -3px",
      "padding": "10px",
      "border-radius": "5px"
    }
    title = f'# TracebTB report. {len(sub)} samples. {date.today()}'
    header = pn.pane.Markdown(title, styles=styles, margin=(5, 20))

    #prepare data
    if len(sub[col].unique())>20:
        cmap=None

    sub['color'],c = tools.get_color_mapping(sub, col, cmap)
    herds = list(sub.HERD_NO)
    mov = tools.get_moves_bytag(sub, moves, lpis_cent)
    if mov is not None:
        herds.extend(mov.move_to)
    sp = parcels[parcels.SPH_HERD_N.isin(herds)].copy()
    sp['color'],c = tools.get_color_mapping(sp, 'SPH_HERD_N')

    #plot map
    plot = bokeh_plot.plot_selection(sub, sp, col=col, legend=True)
    bokeh_plot.plot_moves(plot, mov, lpis_cent)
    #overview
    fig = plot_overview(sub)
    overview = pn.pane.Matplotlib(fig, width=250)
    plt.close(fig)
    #tree
    treefile = get_tree(snpdist, sub.index)
    html = draw_toytree(treefile, sub, tiplabelcol='sample', markercol='Species')
    treepane = pn.pane.HTML(html)
    #table
    cols = ['sample','HERD_NO','Animal_ID','Species','County','Homebred','IE_clade']
    table = pn.pane.DataFrame(sub[cols], styles=CARD_STYLE)
    herdcolors = dict(zip(sp.SPH_HERD_N,sp.color))
    #timeline
    timeline = pn.pane.Bokeh(bokeh_plot.plot_moves_timeline(mov, herdcolors, height=500))

    #main = pn.Column(pn.Row(plot, treepane, overview, height=600),pn.Row(table,timeline))
    main = pn.FlexBox(pn.Row(overview, plot, treepane, height=600), table, timeline)
    report = pn.Column(header, main)
    return report

def layers_from_file(layers_file):
    """load layers from file"""

    import fiona
    layers = {}
    names = fiona.listlayers(layers_file)
    for name in names:
        gdf = gpd.read_file(layers_file, layer=name)
        layers[name] = gdf
    return layers

class Dashboard:

    def __init__(self, meta, parcels, moves=None, lpis_cent=None,
                snpdist=None, lpis_master_file=None, treefile=None,
                selections={}, layers={}):
        """
        Dashboard app with panel for tracebtb.
        Args:
            meta: geodataframe of meta with locations of samples
            parcels: parcels for relevant samples from lpis
            moves: dataframe of moves
            lpis_cent: centroids from LPIS
        """

        self.meta = meta
        self.moves = moves
        self.lpis_cent = lpis_cent
        self.parcels = parcels
        self.snpdist = snpdist
        self.lpis = None
        self.lpis_master_file = lpis_master_file
        self.selections = selections
        self.layers = layers
        if treefile is not None:
            from Bio import Phylo
            self.tree = Phylo.read(treefile, "newick")
            self.treefile = treefile
        self.view_history = []
        self.current_index = 0
        self.layout = self.setup_widgets()
        return

    def show(self):
        return self.layout

    def setup_widgets(self):

        w=140
        pn.config.throttled = True
        #main panes
        self.plot_pane = pn.pane.Bokeh()
        self.overview_pane = pn.pane.Matplotlib(height=300)
        self.split_pane = pn.Column(sizing_mode='stretch_both')
        #pane = pn.Column()
        tccols = ['','snp7','snp5','snp12']
        small_style = """
            .bk-root .bk-select {
                font-size: 10px;
            }
            }"""
        self.timlinecolor_input = pnw.Select(name='Color by', value='',options=tccols,
                                             width=100,stylesheets=[small_style])
        self.timeline_pane = pn.pane.Bokeh(sizing_mode='stretch_height')
        self.info_pane = pn.pane.Markdown('', styles={'color': "red"})

        #main table
        self.meta_pane = pnw.Tabulator(disabled=True,page_size=100,
                                frozen_columns=['sample'],stylesheets=[stylesheet],
                                sizing_mode='stretch_both')
        self.showselected_btn = pnw.Button(name='Select Samples', button_type='primary', align="end")
        pn.bind(self.select_from_table, self.showselected_btn, watch=True)
        self.selectrelated_btn = pnw.Button(name='Find Related', button_type='primary', align="end")
        pn.bind(self.select_related, self.selectrelated_btn, watch=True)
        self.threshold_input = pnw.IntInput(name='Threshold', value=7, step=1, start=2, end=20,width=60)
        #search
        scols = ['sample','Year','HERD_NO','Animal_ID','Species','County','IE_clade','snp7','snp12','snp20']
        self.search_input = pnw.TextInput(name="Search", value='',sizing_mode='stretch_width')
        self.searchcol_select = pnw.Select(name='Column',value='HERD_NO',options=scols,width=100)
        self.search_btn = pnw.Button(icon=get_icon('search'), icon_size='1.8em', align="end")
        pn.bind(self.do_search, self.search_btn, watch=True)

        self.reset_btn = pnw.Button(icon=get_icon('refresh'), icon_size='1.8em', align="end")
        pn.bind(self.reset_table, self.reset_btn, watch=True)

        self.table_widgets = pn.Row(self.showselected_btn, self.selectrelated_btn, self.threshold_input,
                            self.search_input, self.searchcol_select, self.search_btn, self.reset_btn,
                            sizing_mode='stretch_width')
        self.table_pane = pn.Column(self.meta_pane,self.table_widgets,sizing_mode='stretch_both')

        #selected table
        self.selected_pane = pnw.Tabulator(show_index=False,disabled=True,page_size=100,
                                    frozen_columns=['sample'],sizing_mode='stretch_both')
        #moves table
        self.moves_pane = pnw.Tabulator(show_index=False,disabled=True,page_size=100,
                                    frozen_columns=['tag'],sizing_mode='stretch_both')
        #herds table
        self.herds_table = pnw.Tabulator(show_index=False,disabled=True,page_size=100,
                                frozen_columns=['HERD_NO'],sizing_mode='stretch_both')
        self.selectherds_btn = pnw.Button(name='Select Herds', button_type='primary', align="end")
        pn.bind(self.select_herds, self.selectherds_btn, watch=True)
        downloadherds_btn = pnw.FileDownload(callback=pn.bind(self.herds_file),
                                            filename='herds.csv', button_type='success')
        self.herds_pane = pn.Column(self.herds_table, pn.Row(self.selectherds_btn,downloadherds_btn))

        self.clusters_table = pnw.Tabulator(show_index=False,disabled=True,page_size=100,
                                    frozen_columns=['cluster'],sizing_mode='stretch_both')
        self.selectclusters_btn = pnw.Button(name='Select Groups', button_type='primary', align="end")
        pn.bind(self.select_clusters, self.selectclusters_btn, watch=True)
        self.clusters_pane = pn.Column(self.clusters_table, pn.Row(self.selectclusters_btn))

        self.tree_pane = pn.pane.HTML(height=400)
        #self.network_pane = pn.pane.Bokeh()
        #self.snpdist_pane = pn.pane.Bokeh(height=400)
        #details
        self.details_pane = pn.pane.DataFrame(sizing_mode='stretch_both')

        cols = [None]+tools.get_ordinal_columns(self.meta)
        #selections pane
        self.selections_input = pnw.Select(name='Selections',options=list(self.selections.keys()),value='',width=w)
        self.loadselection_btn = pnw.Button(name='Load Selection', button_type='primary',width=w)
        pn.bind(self.load_selection, self.loadselection_btn, watch=True)
        self.deleteselection_btn = pnw.Button(name='Delete Selection', button_type='primary',width=w)
        pn.bind(self.delete_selection, self.deleteselection_btn, watch=True)
        self.saveselection_btn = pnw.Button(name='Save Selection', button_type='primary',width=w)
        pn.bind(self.save_selection, self.saveselection_btn, watch=True)
        self.selectionname_input = pnw.TextInput(value='',width=w)
        upload_label = pn.pane.Markdown("### Import Selections")
        self.upload_selection_input = pnw.FileInput(accept='.json')
        self.upload_selection_input.param.watch(self.import_selections, 'value')
        self.exportselections_btn = pnw.FileDownload(name='Export Selections', callback=pn.bind(self.export_selections),
                                            filename='selections.json', button_type='success')

        card1 = pn.Column('## Selections', pn.Row(pn.Column(self.selections_input, self.selectionname_input),
                        pn.Column(self.loadselection_btn, self.deleteselection_btn, self.saveselection_btn)),
                                    pn.Column(upload_label,self.upload_selection_input,
                                            self.exportselections_btn,styles=dict(background='#ECECEC'),
                                            sizing_mode='stretch_width'),
                        width=340, styles=card_style)
        #layers
        label = pn.pane.Markdown("### Import Layer")
        self.importlayer_input = pnw.FileInput(accept='.zip,.geojson')
        self.importlayer_input.param.watch(self.import_layer, 'value')
        self.layers_select = pn.widgets.MultiSelect(name='Layers to show:', value=[],
                                            options=list(self.layers.keys()), height=100, size=3)
        card2 = pn.Column('## Layers', label, self.importlayer_input, self.layers_select,
                        width=340, styles=card_style)

        #lpis
        self.loadlpis_btn = pnw.Button(name='Load LPIS', button_type='primary')
        pn.bind(self.load_lpis, self.loadlpis_btn, watch=True)
        card3 = pn.Column('## LPIS', self.loadlpis_btn, width=340, styles=card_style)

        #settings
        self.markersize_input = pnw.IntSlider(name='marker size', value=10, start=0, end=80,width=w)
        self.edgewidth_input = pnw.FloatSlider(name='marker edge width', value=0.8, start=0, end=4, step=.1,width=w)
        self.labelsize_input = pnw.IntSlider(name='label size', value=18, start=6, end=80,width=w)
        self.legendsize_input = pnw.IntSlider(name='legend size', value=12, start=6, end=40,width=w)
        self.hexbins_input = pnw.IntSlider(name='hex bins', value=10, start=5, end=100,width=w)
        self.scalebar_input = pnw.Checkbox(name='show scalebar',value=False)
        self.randseed_input = pnw.IntSlider(name='random color seed', value=647, start=1, end=5000,width=w)
        savesettings_btn = pnw.Button(name='Save Settings',button_type='primary')
        pn.bind(self.save_settings, savesettings_btn, watch=True)
        card4 = pn.Column('## Settings', self.markersize_input,self.edgewidth_input,self.labelsize_input,
                          self.legendsize_input, self.hexbins_input, self.scalebar_input, self.randseed_input,
                          savesettings_btn, width=340, styles=card_style)

        #reporting
        self.doreport_btn = pnw.Button(name='Generate', button_type='primary')
        self.savereport_btn = pnw.FileDownload(file=report_file, button_type='success', auto=False,
                                        embed=False, name="Download Report")
        card5 = pn.Column('## Reporting', self.doreport_btn, self.savereport_btn, width=340, styles=card_style)
        pn.bind(self.create_report, self.doreport_btn, watch=True)

        utils_pane = pn.FlexBox(*[card1,card2, card3, card4, card5], flex_direction='column', min_height=400,
                                styles={'margin': '10px'}, sizing_mode='stretch_both')

        #nav buttons
        self.prev_btn = pnw.Button(icon=get_icon('arrow-left'), description='back', icon_size='1.6em',width=60)
        self.next_btn = pnw.Button(icon=get_icon('arrow-right'), description='forward', icon_size='1.6em',width=60)
        pn.bind(self.back, self.prev_btn, watch=True)
        pn.bind(self.forward, self.next_btn, watch=True)
        nav_pane = pn.Row(self.prev_btn,self.next_btn)

        #widgets
        self.groupby_input = pnw.Select(name='group by',options=cols,value='snp7',width=w)
        self.groups_table = pnw.Tabulator(disabled=True, widths={'count': 30}, layout='fit_columns',
                                          pagination=None, height=200, width=w,
                                          stylesheets=[stylesheet])
        self.colorby_input = pnw.Select(name='color by',options=cols,value='snp7',width=w)
        self.cmap_input = pnw.Select(name='colormap',options=colormaps,value='Set1',width=w)
        self.tiplabel_input = pnw.Select(name='tip label',options=list(self.meta.columns),value='sample',width=w)
        self.provider_input = pnw.Select(name='provider',options=['']+bokeh_plot.providers,value='CartoDB Positron',width=w)
        self.pointstyle_input = pnw.Select(name='points display',options=['default','pie'],value='default',width=w)
        widgets = pn.Column(pn.WidgetBox(nav_pane,self.groupby_input,self.groups_table,
                                         self.colorby_input,self.cmap_input,self.tiplabel_input,
                                        self.provider_input,self.pointstyle_input),self.info_pane,width=w+30)
        #button toolbar
        icsize = '1.9em'
        self.split_btn = pnw.Button(icon=get_icon('plot-grid'), description='split view', icon_size=icsize)
        self.selectregion_btn = pnw.Button(icon=get_icon('plot-region'), description='select in region', icon_size=icsize)
        self.selectradius_btn = pnw.Button(icon=get_icon('plot-centroid'), description='select within radius', icon_size=icsize)
        self.tree_btn = pnw.Toggle(icon=get_icon('tree'), icon_size=icsize)
        self.parcels_btn = pnw.Toggle(icon=get_icon('parcels'), icon_size=icsize)
        self.moves_btn = pnw.Toggle(icon=get_icon('plot-moves'), icon_size=icsize)
        self.legend_btn = pnw.Toggle(icon=get_icon('legend'), icon_size=icsize)
        self.neighbours_btn = pnw.Toggle(icon=get_icon('neighbours'), icon_size=icsize)
        self.parcellabel_btn = pnw.Toggle(icon=get_icon('parcel-label'), icon_size=icsize)
        self.showcounties_btn = pnw.Toggle(icon=get_icon('counties'), icon_size=icsize)
        self.kde_btn = pnw.Toggle(icon=get_icon('contour'), icon_size=icsize)
        self.hex_btn = pnw.Toggle(icon=get_icon('hexbin'), icon_size=icsize)
        #pie_btn = pnw.Toggle(icon=get_icon('scatter-pie'), icon_size='1.8em')
        #lockbtn = pnw.Toggle(icon=get_icon('lock'), icon_size='1.8em')
        toolbar = pn.Column(pn.WidgetBox(self.selectregion_btn,self.selectradius_btn,self.tree_btn,
                                        self.parcels_btn,self.parcellabel_btn,self.showcounties_btn,self.moves_btn,self.legend_btn,
                                        self.neighbours_btn,self.kde_btn,self.hex_btn,self.split_btn),width=70)

        #option below plot
        self.timeslider = pnw.IntRangeSlider(name='Time',width=150,
                        start=2000, end=2024, value=(2000, 2024), step=1)
        self.clustersizeslider = pnw.IntSlider(name='Min. Cluster Size',width=150,
                        start=1, end=20, value=1, step=1)
        self.homebredbox = pnw.Checkbox(name='Homebred',value=False)
        self.findrelated_btn = pnw.Button(name='Find Related', button_type='primary', align="end")
        pn.bind(self.find_related, self.findrelated_btn, watch=True)
        self.related_col_input = pnw.Select(options=cols,value='snp7',width=90)

        filters = pn.Row(self.findrelated_btn,self.related_col_input,self.timeslider,self.clustersizeslider,self.homebredbox)

        self.groupby_input.param.watch(self.update_groups, 'value')
        self.groups_table.param.watch(self.select_group, 'selection')
        self.provider_input.param.watch(self.set_provider, 'value')
        self.colorby_input.param.watch(self.update, 'value')
        self.cmap_input.param.watch(self.update, 'value')
        self.tiplabel_input.param.watch(self.update, 'value')
        self.pointstyle_input.param.watch(self.update, 'value')
        self.tree_btn.param.watch(self.update, 'value')
        self.parcels_btn.param.watch(self.update, 'value')
        self.parcellabel_btn.param.watch(self.update, 'value')
        self.showcounties_btn.param.watch(self.update, 'value')
        self.moves_btn.param.watch(self.update, 'value')
        self.legend_btn.param.watch(self.update, 'value')
        self.neighbours_btn.param.watch(self.update, 'value')
        self.kde_btn.param.watch(self.update, 'value')
        self.hex_btn.param.watch(self.update, 'value')
        self.timeslider.param.watch(self.update, 'value')
        #pn.bind(update, timeslider.param.value_throttled)
        self.clustersizeslider.param.watch(self.update, 'value')
        self.homebredbox.param.watch(self.update, 'value')

        #pn.bind(update_tree, treebtn, watch=True)
        pn.bind(self.split_view, self.split_btn, watch=True)
        pn.bind(self.select_region, self.selectregion_btn, watch=True)
        pn.bind(self.select_radius, self.selectradius_btn, watch=True)

        #categorical plot widget
        self.catplot_pane = pn.pane.Matplotlib(sizing_mode='stretch_both')
        self.catx_input = pnw.Select(name='x',options=cols,value='Species',width=w)
        self.caty_input = pnw.Select(name='y',options=cols,value=None,width=w)
        self.cathue_input = pnw.Select(name='hue',options=cols,value=None,width=w)
        kinds = ['count','bar','strip','swarm']
        self.catkind_input = pnw.Select(name='kind',options=kinds,
                                value='count',width=w)
        self.cathue_input = pnw.Select(name='hue',options=cols,value='County',width=w)
        catupdate_btn = pnw.Button(name='Update', button_type='primary',align='end')
        pn.bind(self.update_catplot, catupdate_btn, watch=True)
        self.analysis_pane1 = pn.Column(pn.Row(self.catx_input,self.caty_input,
                                       self.cathue_input,self.catkind_input,catupdate_btn),
                                self.catplot_pane)

        styles={ "margin-top": "10px", "font-size": "15px"}
        self.about_pane = pn.pane.Markdown('',styles=styles)
        self.about()

        app = pn.Column(
                    pn.Row(widgets,
                        toolbar,
                    pn.Column(
                        pn.Tabs(('Map',pn.Column(self.plot_pane,filters)),
                                ('Split View',self.split_pane),
                                ('Summary I',self.analysis_pane1),
                                ('Main Table', self.table_pane),
                                ('Selected', self.selected_pane),
                                ('Moves',self. moves_pane),
                                ('Herds',self.herds_pane),('Groups',self.clusters_pane),
                                ('Tools', utils_pane),
                                ('About',self.about_pane),
                                dynamic=True,
                                sizing_mode='stretch_both'),
                            ),
                    pn.Tabs(('Overview',pn.Column(self.overview_pane,
                                            pn.Column(self.timeline_pane,self.timlinecolor_input), width=500)),
                                    ('Tree',self.tree_pane),
                                    ('Details',self.details_pane),
                                    dynamic=True,width=500),
                max_width=2600,min_height=600
        ))
        app.sizing_mode='stretch_both'
        self.meta_pane.value = self.meta
        self.update_groups()
        self.selected = self.meta[self.meta.snp7.isin(['24'])].copy()
        self.update(sub=self.selected)
        return app

    def update(self, event=None, sub=None):
        """
        Update selection - does all the display.
        Args:
            sub: a selection may be provided and will set current 'selected' variable
        """

        provider = self.provider_input.value
        cmap = self.cmap_input.value
        self.col = col = self.colorby_input.value
        self.ms = ms = self.markersize_input.value
        lw = self.edgewidth_input.value
        legend = self.legend_btn.value

        if sub is None:
            sub = self.selected

        if len(sub[col].unique()) > 20:
            cmap = None
        elif len(sub[col].unique()) > 10:
            cmap = 'tab20'
        seed = self.randseed_input.value
        sub['color'], cm1 = tools.get_color_mapping(sub, col, cmap, seed=seed)
        self.selected = sub
        # Apply filters
        sub = self.apply_filters(sub)
        self.info_pane.object = f'**{len(sub)} samples**'
        self.update_overview(sub=sub)

        # Get the herds and moves
        herds = list(sub.HERD_NO)
        mov = tools.get_moves_bytag(sub, self.moves, self.lpis_cent)

        if mov is not None and self.moves_btn.value is True:
            herds.extend(mov.move_to)

        sp = self.parcels[self.parcels.SPH_HERD_N.isin(herds)].copy()
        pcmap = 'tab20'
        if len(sp.SPH_HERD_N.unique()) > 20:
            pcmap = None

        sp['color'], cm2 = tools.get_color_mapping(sp, 'SPH_HERD_N', pcmap)

        if self.parcels_btn.value is True:
            pcls = sp
            if mov is not None:
                herds.extend(mov.move_to)
        else:
            pcls = None

        labels = self.parcellabel_btn.value
        labelsize = self.labelsize_input.value
        legsize = self.legendsize_input.value

        if self.pointstyle_input.value == 'default':
            p = bokeh_plot.plot_selection(sub, pcls, provider=provider, ms=ms, lw=lw,
                                            col=col, legend=legend, labels=labels,
                                            legend_fontsize=legsize, label_fontsize=labelsize,
                                            scalebar=self.scalebar_input.value)

            from bokeh.events import Tap
            p.on_event(Tap, self.point_selection)
        else:
            p = bokeh_plot.scatter_pie(sub, 'HERD_NO', col, cm1, provider=provider,
                                            legend=legend, legend_fontsize=legsize,
                                            scalebar=self.scalebar_input.value)

        if self.showcounties_btn.value is True:
            bokeh_plot.plot_counties(p)
        self.show_layers(p)

        if self.hex_btn.value is True:
            bins = self.hexbins_input.value
            bokeh_plot.hexbin(sub, n_bins=bins, p=p)

        if self.kde_btn.value is True:
            bokeh_plot.kde_plot_groups(sub, p, col, 6)

        if self.neighbours_btn.value is True and self.lpis is not None:
            shb = shared_borders(sp, self.lpis)
            bokeh_plot.plot_lpis(shb, p)

        mov = tools.get_moves_bytag(sub, self.moves, self.lpis_cent)

        if self.moves_btn.value is True:
            bokeh_plot.plot_moves(p, mov, self.lpis_cent)

        if mov is not None:
            self.moves_pane.value = mov.reset_index().drop(columns=['geometry'])
        else:
            self.moves_pane.value = pd.DataFrame()

        self.plot_pane.object = p

        self.selected_pane.value = sub.drop(columns=['geometry'])

        def highlight(x):
            color = self.speciescolors[x]
            return 'background-color: %s' % color

        self.selected_pane.style.apply(highlight, subset=['Species'], axis=1)

        # Timeline plot
        tlc = self.timlinecolor_input.value
        tldf = bokeh_plot.get_timeline_data(mov, self.meta)
        #get colors for timeline data
        if tldf is not None:
            if tlc in ['',None]:
                #match parcel colors by default
                tlcolors = dict(zip(sp.SPH_HERD_N, sp.color))
                tldf['color'] = tldf.move_to.map(tlcolors).fillna('grey')
            else:
                #tlcolors = dict(zip(sub['sample'], sub.color))
                tldf['color'],c = tools.get_color_mapping(tldf, tlc, cmap)
                #tldf['color'] = tldf['sample'].map(tlcolors).fillna('grey')

            p = bokeh_plot.plot_moves_timeline(tldf)
            self.timeline_pane.object = p
            self.timeline_pane.param.trigger('object')

        if self.tree_btn.value is True:
            self.update_tree(sub=sub)

        # Update summaries
        self.update_herd_summary()
        self.update_cluster_summary()
        return

    def point_selection(self, event):
        """Point click callback"""

        p = event.model
        source = p.renderers[1].data_source
        indices = source.selected.indices
        if len(indices)==0:
            return
        idx = indices[0]
        geojson_dict = json.loads(source.geojson)
        selected_feature = geojson_dict['features'][idx]
        name = selected_feature['id']
        data = self.meta.loc[name]
        #print (data)
        self.details_pane.object = pd.DataFrame(data)
        return

    def update_herd_summary(self, event=None):
        """herd summary"""

        h = tools.herd_summary(self.selected, self.moves, self.snpdist)
        self.herds_table.value = h
        return

    def update_cluster_summary(self, event=None):

        col = self.colorby_input.value
        cl = tools.cluster_summary(self.selected, col, 5, self.snpdist)
        self.clusters_table.value = cl
        return

    def update_overview(self, event=None, sub=None):

        fig = plot_overview(self.selected)
        self.overview_pane.object = fig
        plt.close(fig)
        return

    def load_lpis(self, event=None):

        if self.lpis_master_file != None:
            self.lpis = gpd.read_file(self.lpis_master_file).set_crs('EPSG:29902')
        return

    def apply_filters(self, df):
        """Filter from widgets"""

        minsize = self.clustersizeslider.value
        key = self.colorby_input.value
        groups = df[key].value_counts()
        groups = groups[groups>=minsize]
        df = df[df[key].isin(groups.index)].copy()

        start, end = self.timeslider.value
        df = df[(df.Year>=start) & (df.Year<=end) | (df.Year.isnull())].copy()
        if self.homebredbox.value == True:
            df = df[df.Homebred=='yes'].copy()
        return df

    def set_provider(self, event=None):
        """Change map provider"""

        p = self.plot_pane.object
        #remove old tile
        p.renderers = [x for x in p.renderers if not str(x).startswith('TileRenderer')]
        provider = self.provider_input.value
        if provider in bokeh_plot.providers:
            p.add_tile(provider, retina=True)
        return

    def update_groups(self, event=None):
        """Change group choices"""

        groupby = self.groupby_input.value
        vals = pd.DataFrame(self.meta[groupby].value_counts())
        vals = vals[vals['count']>1]
        self.groups_table.value = vals
        return

    def select_group(self, event=None):
        """Select groups from table"""

        df = self.groups_table.value
        rows = self.groups_table.selection
        groups = list(df.iloc[rows].index)
        #print (groups)
        key = self.groupby_input.value
        sub = self.meta[self.meta[key].isin(groups)].copy()
        self.update(sub=sub)
        self.add_to_history()
        return

    def select_region(self, event=None):
        """Select samples in region"""

        p = self.plot_pane.object
        xmin, xmax = p.x_range.start, p.x_range.end
        ymin, ymax = p.y_range.start, p.y_range.end
        from pyproj import Transformer
        transformer_to_latlon = Transformer.from_crs("EPSG:3857", "EPSG:29902", always_xy=True)
        xmin, ymin = transformer_to_latlon.transform(xmin, ymin)
        xmax, ymax = transformer_to_latlon.transform(xmax, ymax)
        sub = self.meta.cx[xmin:xmax, ymin:ymax]
        self.update(sub=sub)
        self.add_to_history()
        return

    def select_radius(self, event=None):
        """Select within radius of a center"""

        point = self.selected.union_all().centroid
        distances = self.selected.geometry.apply(lambda x: x.distance(point))
        radius = distances.median()
        sub = self.meta[self.meta.geometry.distance(point) <= radius].copy()
        self.update(sub=sub)
        self.add_to_history()
        return

    def select_from_table(self, event=None):
        """Select samples from table"""

        df = self.meta_pane.selected_dataframe
        sub = self.meta.loc[df.index]
        self.update(sub=sub)
        self.add_to_history()
        return

    def select_related(self, event=None, df=None):
        """Find related samples to selected indexes i.e. within n snps"""

        df = self.meta_pane.selected_dataframe
        idx = list(df.index)
        dist = self.threshold_input.value
        names=[]
        for i in idx:
            found = tools.get_within_distance(self.snpdist, i, dist)
            names.extend(found)
        names = list(set(names))
        sub = self.meta.loc[names]
        self.update(sub=sub)
        self.add_to_history()
        return

    def find_related(self, event=None):
        """Find related samples"""

        col = self.related_col_input.value
        cl = list(self.selected[col].unique())
        sub = self.meta[(self.meta[col].isin(cl))].copy()
        self.update(sub=sub)
        self.add_to_history()
        return

    def select_herds(self, event=None):

        df = self.herds_table.selected_dataframe
        sub = self.meta[self.meta.HERD_NO.isin(df.HERD_NO)].copy()
        self.update(sub=sub)
        self.add_to_history()
        return

    def herds_file(self):
        """Return herds table for export"""

        df = self.herds_table.value
        sio = io.StringIO()
        df.to_csv(sio,index=False)
        sio.seek(0)
        return sio

    def select_clusters(self, event=None):

        df = self.clusters_table.selected_dataframe
        col = self.colorby_input.value
        sub = self.meta[(self.meta[col].isin(df.cluster))].copy()
        self.update(sub=sub)
        return

    def split_view(self, event=None):
        """split view"""

        provider = self.provider_input.value
        cmap = self.cmap_input.value
        col = self.colorby_input.value
        labels = self.parcellabel_btn.value
        labelsize = self.labelsize_input.value
        legsize = self.legendsize_input.value
        ms = self.markersize_input.value
        lw = self.edgewidth_input.value
        kde = self.kde_btn.value

        sub = self.apply_filters(self.selected)
        f = bokeh_plot.split_view(sub, col, self.parcels, provider, labels=labels,
                                    legend_fontsize=legsize, label_fontsize=labelsize,
                                    limit=9, ms=ms, lw=lw, kde=kde)
        self.split_pane.objects.clear()
        self.split_pane.append(pn.pane.Bokeh(f))
        return

    def update_tree(self, event=None, sub=None):
        """Update with toytree"""

        if len(sub)<=1 or len(sub)>1800:
            html = '<h1><2 or too many samples</h1>'
            self.tree_pane.object = html
            return
        if self.tree == None:
            return
        from Bio import Phylo
        stree = keep_tips(self.tree, list(sub.index))
        tempfile = 'temp.newick'
        Phylo.write(stree, tempfile, "newick")

        col = self.colorby_input.value
        html = draw_toytree(tempfile, sub,
                            tiplabelcol=self.tiplabel_input.value,
                            markercol='Species', node_hover=True,
                            height=700)
        #html = phylocanvas_tree(treefile, sub, col)
        self.tree_pane.object = html
        return

    def do_search(self, event=None):
        """Search main table"""

        query = self.search_input.value
        col = self.searchcol_select.value
        found = self.meta[self.meta[col]==query]
        self.meta_pane.value = found.drop(columns=['geometry'])
        return

    def reset_table(self, event=None):
        self.meta_pane.value = self.meta
        return

    def create_report(self, event=None):
        """Standard report"""

        sub = self.selected
        sp = self.parcels[self.parcels.SPH_HERD_N.isin(sub.HERD_NO)].copy()
        sp['color'],c = tools.get_color_mapping(sp, 'SPH_HERD_N')
        mov = tools.get_moves_bytag(sub, self.moves, self.lpis_cent)
        col = self.colorby_input.value
        result = report(self.selected, sp, mov, col, self.lpis_cent, self.snpdist)
        #report_file = f'report_{datetime.today().date()}.html'
        result.save(report_file)
        return

    def load_selection(self, event=None):
        """Load a selection"""

        name = self.selections_input.value
        idx = self.selections[name]['indexes']
        self.selected = self.meta.loc[idx]
        self.update(sub=self.selected)
        return

    def delete_selection(self, event=None):
        """Delete a selection"""

        name = self.selections_input.value
        del (self.selections[name])
        self.selections_input.options = list(self.selections.keys())
        with open(selections_file,'w') as f:
            f.write(json.dumps(self.selections))
        return

    def import_selections(self, event=None):
        """Import a selection"""

        data = self.upload_selection_input.value
        #print(filename)
        new = json.load(data)
        self.selections.update(new)
        return

    def export_selections(self, event=None):
        """Selections export"""

        sio = io.StringIO()
        #df.to_csv(sio,index=False)
        json.dump(self.selections, sio)
        sio.seek(0)
        return sio

    def save_selection(self, event=None):
        """Save a selection"""

        name = self.selectionname_input.value
        if name == '':
            print ('empty name')
            return
        self.selections[name] = {}
        self.selections[name]['indexes'] = list(self.selected.index)
        #add to menu
        self.selections_input.options = list(self.selections.keys())
        #save to file
        with open(selections_file,'w') as f:
            f.write(json.dumps(self.selections))
        return

    def import_layer(self, event=None):
        """Import shapefile"""

        fname = self.importlayer_input.filename
        ext = os.path.splitext(fname)[1]
        data = self.importlayer_input.value
        tempf='temp'+ext
        with open(tempf, 'wb') as f:
            f.write(data)
        name = os.path.basename(fname)
        name = os.path.splitext(name)[0]
        self.add_layer(tempf, name)
        return

    def add_layer(self, filename, name):
        """Add a map layer"""

        gdf = gpd.read_file(filename).to_crs('EPSG:3857')
        self.layers[name] = gdf
        self.layers_select.options = list(self.layers.keys())
        #save layer to file
        gdf.to_file(layers_file, layer=name, driver="GPKG", mode="a")
        return

    def show_layers(self, p):
        """Show visible layers"""

        show = self.layers_select.value
        colors = tools.random_colors(len(self.layers.keys()), 10)
        i=0
        for l in self.layers:
            if l in show:
                gdf = self.layers[l]
                bokeh_plot.plot_gdf(gdf, p, line_width=2, line_color=colors[i])
                i+=1
        return

    def add_to_history(self):
        """Add current selection to history"""

        view_history = self.view_history
        view_history.append(self.selected.index)
        if len(view_history) > 20:
            view_history.pop(0)
        # Move the current index to the new end
        self.current_index = len(view_history)-1
        return

    def back(self, event=None):
        """Go back"""

        view_history = self.view_history
        if len(view_history) == 0:
            return
        if self.current_index <= 0:
            return
        self.current_index -= 1
        idx = view_history[self.current_index]
        sub = self.meta.loc[idx]
        self.update(sub=sub)
        return

    def forward(self, event=None):

        view_history = self.view_history
        if len(view_history) == 0 or self.current_index >= len(view_history)-1:
            return
        self.current_index += 1
        idx = view_history[self.current_index]
        sub = self.meta.loc[idx]
        self.update(sub=sub)
        return

    def update_catplot(self, event=None, sub=None):
        """Update categorical plot"""

        if sub is None:
            sub = self.selected
        x = self.catx_input.value
        y = self.caty_input.value
        hue = self.cathue_input.value
        kind = self.catkind_input.value
        cmap = self.cmap_input.value
        row = None
        import seaborn as sns
        cg = sns.catplot(data=sub, x=x, y=y, hue=hue, kind=kind, aspect=2,
                        palette=cmap, dodge=True)
        self.catplot_pane.object = cg.fig
        plt.close(cg.fig)
        #print (fig)
        return

    def save_settings(self, event=None):
        """Save settings"""

        #print (self.settings)
        #d = {'colorby':self.col, 'ms':self.ms}
        #self.settings.update(d)
        #with open(configfile, "w") as outfile:
        #    json.dump(self.settings, outfile)
        return

    def current_data_info(self):
        """Currently loaded data summary"""

        m=f"{len(self.meta)} rows loaded\n"
        empty = len(self.meta[self.meta.geometry.is_empty])
        m+=f"{empty} rows with no geometry\n"
        return m

    def about(self):
        try:
            from . import core
            VERSION = core.git_version()
        except:
            from . import __version__ as VERSION
        pandasver = pd.__version__
        gpdver = gpd.__version__
        import bokeh
        bokehver = bokeh.__version__
        m="# TracebTB\n"
        m+="This software is developed as part of a DAFM PSSRC grant (2022PSS113)\n"
        m+="Licensed under the GPL v.3.0\n"
        m+="## Software\n"
        m+=f"* Version {VERSION}\n"
        m+=f"* pandas v{pandasver}\n"
        m+=f"* panel v{pn.__version__}\n"
        m+=f"* bokeh v{bokehver}\n"
        m+=f"* geopandas v{gpdver}\n"
        m+="## Links\n"
        m+="* [Homepage](https://github.com/dmnfarrell/tracebtb)\n"

        m+="## Current data\n"
        m+=self.current_data_info()
        self.about_pane.object = m
        return

def test_app():
    """test app"""
    bootstrap = pn.template.BootstrapTemplate(
        title='Testing'
    )

    def update_plot(event=None):
        #plot_pane.object = bokeh_plot.random_circles(n=20)
        plot_pane.object = bokeh_plot.random_polygons(n=ninput.value)

    plot_pane = pn.pane.Bokeh(bokeh_plot.random_circles(n=20))
    button = pnw.Button(name="TEST", button_type="primary")
    button.on_click(update_plot)
    ninput = pnw.IntInput(name='n', value=10, step=1, start=2, end=1000,width=60)
    app = pn.Column(plot_pane, pn.Row(button,ninput))

    bootstrap.main.append(app)
    bootstrap.servable()
    return bootstrap

def main():
    "Run the application"
    from argparse import ArgumentParser
    parser = ArgumentParser(description='TracebTB')
    parser.add_argument("-p", "--proj", dest="project",default=None,
                            help="load project file", metavar="FILE")
    parser.add_argument("-s", "--settings", dest="settings",default=configfile,
                            help="load a json settings file", metavar="FILE")
    parser.add_argument("-t", "--test", dest="test", action="store_true",
                        help="dummy test app")
    args = parser.parse_args()

    if args.test == True:
        pn.serve(test_app, port=5010, prefix='testapp',
                 websocket_origin=["localhost:5010"])
    elif args.project == None:
        print ('please provide a project file')
        exit()
    else:
        #load config file
        if not os.path.exists(args.settings):
            with open(configfile, "w") as outfile:
                json.dump(defaults, outfile)
            lpis_master_file = None
            treefile = None
        else:
            print('found settings file')
            with open(args.settings) as f:
                settings = json.load(f)['dashboard']
                print (settings)
            lpis_master_file = settings['lpis_master_file']
            if 'tree_file' in settings:
                treefile = settings['tree_file']
            else:
                treefile = None

        data = pickle.load(open(args.project,'rb'))
        meta = data['meta']
        moves = data['moves']
        lpis_cent = data['lpis_cent']
        parcels = data['parcels']
        snpdist = data['snpdist']

        if os.path.exists(selections_file):
            selections = json.load(open(selections_file,'r'))
        else:
            selections = {}
        if os.path.exists(layers_file):
            layers = layers_from_file(layers_file)
        else:
            layers = {}

        def create_app():
            # Create a session-specific app function
            bootstrap = pn.template.BootstrapTemplate(
                title='TracebTB',
                favicon=logoimg,
                logo=logoimg,
                header_color='white'
            )
            #Generate a new dashboard instance per session
            app = Dashboard(meta, parcels, moves, lpis_cent, snpdist, lpis_master_file,
                            treefile, selections, layers)
            app.project_file = args.project
            app.settings = settings
            app.treefile = treefile
            layout = app.show()
            bootstrap.main.append(layout)
            bootstrap.servable()
            return bootstrap

        if 'reverse_proxy' in settings:
            #using nginx proxy with basic auth
            s = settings['reverse_proxy']
            pn.serve(create_app, port=5010, prefix=s['prefix'],
                websocket_origin=s['origin'], #basic_auth='credentials.json',
                cookie_secret='cookie_secret')
        else:
            #default
            pn.serve(create_app, port=5010, websocket_origin=["localhost:5010"],
                     cookie_secret='cookie_secret')

if __name__ == '__main__':
    main()
