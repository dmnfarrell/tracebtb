#!/usr/bin/env python3

"""
    TraceBTB panel dashboard
    Created June 2024
    Copyright (C) Damien Farrell
"""

import sys,os,time,re
from datetime import date
import platform
from datetime import datetime
import pickle
import glob,io
import json
import math
import pylab as plt
import numpy as np
import pandas as pd
import geopandas as gpd
import seaborn as sns

from bokeh.io import show, curdoc
from bokeh.plotting import figure
from bokeh.models import Range1d, CustomJS, TapTool, Div
from bokeh.events import Tap
from bokeh.layouts import column
import panel as pn
import panel.widgets as pnw

pn.extension('tabulator')
pn.config.throttled = True

from tracebtb import tools, plotting, trees, bokeh_plot, source_attribution, movement

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

report_file = 'report.html'
selections_file = os.path.join(configpath,'selections.json')
layers_file = os.path.join(configpath,'layers.gpkg')

speciescolors = {'Bovine':'blue','Badger':'red','Ovine':'green'}
colormaps = ['Paired', 'Dark2', 'Set1', 'Set2', 'Set3',
            'tab10', 'tab20', 'tab20b', 'tab20c', 'RdBu']
card_style = {
    'background': '#f9f9f9',
    'border': '1px solid #bcbcbc',
    'padding': '5px',
    'margin': '5px',
    'box-shadow': '4px 4px 4px #bcbcbc'
}
df_stylesheet = """
tr {
    font-size: 14px;
}
"""
stylesheet = """
.tabulator-cell {
    font-size: 12px;
}
.tabulator-col-title {
    font-size: 11px;
}
"""
dm_css = """
.tabulator .tabulator-header .tabulator-col .tabulator-col-content .tabulator-col-title {
  transform: rotate(-90deg);
  white-space: normal;
  text-overflow: clip;
  font-size: 10px;
  width: 20px;
}
.tabulator-cell {
    font-size: 10px;
    width: 10px;
}
"""

icsize = '1.9em'
defaults = {'dashboard':{'lpis_master_file':'','tree_file':None}}
mainsearchcols = ['sample','Year','HERD_NO','Animal_ID','Species','County','Region','class_2021',
         'lineage','strain_name','short_name','last_move_type']

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

def get_figure_coords(p):
    """Coords of plot figure"""

    xmin, xmax = p.x_range.start, p.x_range.end
    ymin, ymax = p.y_range.start, p.y_range.end
    return xmin, xmax, ymin, ymax

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
    html = ''
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

def get_bokeh_widgets(panel_obj):
    """
    Recursively retrieves *all* underlying Bokeh widget models from a Panel layout.
    Returns a dict of {widget_name: bokeh_widget_model}.
    """

    from bokeh.models import Widget as BkWidget, Model
    from bokeh.document import Document
    from bokeh.core.property.descriptors import UnsetValueError

    doc = Document()
    root = panel_obj.get_root(doc)

    widgets = {}

    def walk(model):
        # Store the model if it's a Bokeh widget
        if isinstance(model, BkWidget):
            name = getattr(model, "name", None) or model.__class__.__name__
            widgets[name] = model

        if isinstance(model, Model):
            for attr_name in model.properties():
                try:
                    val = getattr(model, attr_name)
                except UnsetValueError:
                    continue  # Skip properties that aren't set yet
                except Exception:
                    continue  # Skip any weird attributes safely

                # Recurse into any nested Bokeh models or lists
                if isinstance(val, Model):
                    walk(val)
                elif isinstance(val, (list, tuple, set)):
                    for v in val:
                        if isinstance(v, Model):
                            walk(v)

    walk(root)
    return widgets

# --- JS to save all widget states into localStorage ---
save_js_code = """
// Collect all widget values into a settings object
console.log("Saving settings...");
const settings = {};
for (const [name, widget] of Object.entries(widgets)) {
    if (widget.type === 'CheckboxGroup') {
        settings[name] = widget.active;
    } else {
        settings[name] = widget.value;
    }
}
localStorage.setItem('tracebtb_settings', JSON.stringify(settings));
"""

# --- JS to load settings from localStorage ---
load_js_code = """
const saved = localStorage.getItem('tracebtb_settings');
if (saved) {
    const settings = JSON.parse(saved);
    for (const [name, widget] of Object.entries(widgets)) {
        if (settings.hasOwnProperty(name)) {
            if (widget.type === 'CheckboxGroup') {
                widget.active = settings[name];
            } else {
                widget.value = settings[name];
            }
        }
    }
}
"""

class Dashboard:

    def __init__(self, meta, parcels, moves=None, lpis_cent=None,
                snpdist=None, lpis=None, treefile=None,
                testing=None, feedlots=None, settings={},
                selections={}, layers={},
                parent=None,
                **kwargs):
        """
        Base class for dashboard app with panel for tracebtb.
        Args:
            meta: geodataframe of meta with locations of samples
            parcels: parcels for relevant samples from lpis
            moves: dataframe of moves
            lpis_cent: centroids from LPIS
            snpdist: snp distance matrix
            treefile: newick tree file
            testing: testing data
        """

        self.meta = meta
        self.badger = self.meta[self.meta.Species=='Badger']
        self.moves = moves
        self.lpis_cent = lpis_cent
        self.lpis = lpis
        self.parcels = parcels
        self.snpdist = snpdist
        self.selections = selections
        self.layers = layers
        self.hmoves = None
        self.settings = settings
        self.treefile = treefile
        self.parent = parent
        if treefile is not None and os.path.exists(treefile):
            from Bio import Phylo
            self.tree = Phylo.read(treefile, "newick")
            self.treefile = treefile
        else:
            print ('no tree found')
            self.tree = None
        te = self.testing = testing.set_index('HERD_NO')
        self.sr = te.filter(regex='^Sr')
        self.view_history = []
        self.current_index = 0
        self.cols = [None]+tools.get_ordinal_columns(self.meta)+['snp1','snp2','snp3','snp5','snp7','snp12']
        self.layout = self.setup_widgets()
        return

    def load_lpis(self, event=None):
        """Load LPIS"""

        if self.lpis_master_file is None:
            print ('no LPIS file')
        if self.lpis_master_file != None:
            self.lpis = gpd.read_file(self.lpis_master_file).set_crs('EPSG:29902')
            print ('loaded LPIS file')
        return

    def setup_widgets(self):
        """Create widgets - override this"""
        return

    def search_widgets(self, recent_items=8):
        """Quick search widgets"""

        w=140
        self.search_input = pnw.TextInput(name="Query", value='', width=w)
        self.search_btn = pnw.Button(name='Search', icon=get_icon('search'), icon_size='1.8em', width=w)
        pn.bind(self.quick_search, self.search_btn, watch=True)
        #self.search_input.bind('return', self.quick_search)
        self.recents_select = pnw.Select(name='Recent',value='',options=[],width=w,size=recent_items)
        self.recents_select.param.watch(self.recent_search,'value')
        herds_with_samples = list(self.meta.HERD_NO.dropna().sort_values().unique())
        widgets = pn.WidgetBox(self.search_input,self.search_btn,
                               self.recents_select,width=w+20)
        return widgets

    def group_widgets(self):
        """Groupby widgets"""

        w=140
        self.groupby_input = pnw.Select(name='group by',options=self.cols,value='short_name',width=w)
        self.groups_table = pnw.Tabulator(disabled=True, widths={'count': 30}, layout='fit_columns',
                                          pagination=None, height=200, width=w,
                                          stylesheets=[stylesheet])
        return pn.WidgetBox(self.groupby_input,self.groups_table)

    def recent_search(self, event=None):
        """Do recent search"""

        query = self.recents_select.value
        self.quick_search(query=query)
        return

    def quick_search(self, event=None, query=None):
        """Search herd query string"""

        if query == None:
            query = self.search_input.value
        if len(query)<=1:
            return

        #found = self.meta[self.meta.isin([query]).any(axis=1)]
        #found = self.meta[self.meta.map(lambda x: str(query).lower() in str(x).lower()).any(axis=1)]
        query_list = [q.strip().lower() for q in query.split(",")]
        found = self.meta[self.meta.map(lambda x: any(q in str(x).lower() for q in query_list)).any(axis=1)]
        if len(found)>100 or len(found)==0:
            return

        self.add_to_recent(query)
        self.update(sub=found)
        self.add_to_history()
        return

    def sampled_herd_search(self, event=None):
        """Search for herd selection"""

        query = str(self.sampled_herds_select.value)
        self.quick_search(query=query)
        return

    def set_provider(self, event=None):
        """Change map provider"""

        p = self.plot_pane.object
        #remove old tile
        p.renderers = [x for x in p.renderers if not str(x).startswith('TileRenderer')]
        provider = self.provider_input.value
        if provider in bokeh_plot.providers:
            p.add_tile(provider, retina=True)
        return

    def add_to_recent(self, query):
        """Add query to recents list"""

        x = list(self.recents_select.options)
        if query not in x:
            x = [query] + x
        self.recents_select.options = x
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

    def update_tree(self, event=None, sub=None, col='snp7',
                    tip_size=12, font_size='11pt', labelcol='name'):
        """Update tree"""

        if self.tree is None or len(sub)<=1 or len(sub)>4000:
            p = self.tree_pane.object = bokeh_plot.error_message('')
        else:
            from Bio import Phylo
            stree = keep_tips(self.tree, list(sub.index))
            tempfile = 'temp.newick'
            Phylo.write(stree, tempfile, "newick")
            p = bokeh_plot.plot_phylogeny(stree, sub,
                                          tip_size=tip_size, font_size=font_size,
                                          labelcol=labelcol)

        self.tree_pane.objects.clear()
        self.tree_pane.append(pn.pane.Bokeh(p))
        p.on_event(Tap, self.tip_selected)
        return

    def random_herd(self):
        return self.meta.sample(1).iloc[0].HERD_NO

    def tip_selected(self, event):
        """Point click callback"""
        return

    def update(self, event=None, sub=None):
        """Override this if needed"""

        return

    def plot_neighbours(self, df, parcels, column='SPH_HERD_N', pad=.3, legend=False):
        """Show herd parcels and its neighbours - used in multiple dashboards"""

        if len(parcels)>0:
            point = parcels.to_crs('EPSG:3857').iloc[0].geometry.centroid
            x1,y1,x2,y2 = df.union_all().bounds
            pad = (x2-x1)*pad
            if parcels[column].dtype == object:
                cmap=None
            else:
                cmap='viridis'
                legend = True #always need colorbar
            parcels['color'], cm = tools.get_color_mapping(parcels, column, cmap)
            p = bokeh_plot.plot_lpis(parcels, fill_alpha=0.5, line_width=0.2)
            if legend == True:
                bokeh_plot.add_legend(p, parcels, column, fontsize=12)
        else:
            p = bokeh_plot.init_figure()
        df['color'] = 'blue'
        bokeh_plot.plot_lpis(df, p, fill_alpha=0.7, line_width=3)
        return p

    def show(self):
        return self.layout

    def current_data_info(self):
        """Currently loaded data summary"""

        m=f"{len(self.meta)} rows loaded\n"
        empty = len(self.meta[self.meta.geometry.is_empty])
        m+=f"{empty} rows with no geometry\n"
        return m

    def about(self):
        """About details"""

        try:
            from . import core
            VERSION = core.git_version()
        except:
            from . import __version__ as VERSION
        pandasver = pd.__version__
        gpdver = gpd.__version__
        import polars as pl
        polarsver = pl.__version__
        import bokeh
        bokehver = bokeh.__version__

        m="# TracebTB\n"
        m+="This software has been developed as part of a DAFM PSSRC grant (2022PSS113)\n"
        m+="Licensed under the GPL v.3.0\n"
        m+="## Software\n"
        m+=f"* Version {VERSION}\n"
        m+=f"* pandas v{pandasver}\n"
        m+=f"* panel v{pn.__version__}\n"
        m+=f"* bokeh v{bokehver}\n"
        m+=f"* geopandas v{gpdver}\n"
        m+=f"* polars v{polarsver}\n"
        m+="## Links\n"
        m+="* [Homepage](https://github.com/dmnfarrell/tracebtb)\n"

        m+="## Current data\n"
        m+=self.current_data_info()
        #m+=f'![tracebtb]({logoimg})'
        self.about_pane.object = m
        return

class FullDashboard(Dashboard):
    """Full dashboard"""

    def setup_widgets(self):
        """Create widgets"""

        w=140
        #main panes
        self.plot_pane = pn.pane.Bokeh()
        self.overview_pane = pn.pane.Matplotlib(height=300)
        self.split_pane = pn.Column(sizing_mode='stretch_both')

        cols = self.cols
        tccols = ['','snp7','snp5']
        small_style = """
            .bk-root .bk-select {
                font-size: 10px;
            }
            }"""
        self.timelinecolor_input = pnw.Select(name='Color by', value='',options=tccols,
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
        self.threshold_input = pnw.IntInput(name='Threshold', value=5, step=1, start=1, end=20,width=60)
        #search  table
        self.search_input = pnw.TextInput(name="Search", value='',sizing_mode='stretch_width')
        self.searchcol_select = pnw.Select(name='Column',value='HERD_NO',options=mainsearchcols,width=100)
        self.search_btn = pnw.Button(icon=get_icon('search'), icon_size='1.8em', align="end")
        pn.bind(self.do_search, self.search_btn, watch=True)

        self.reset_btn = pnw.Button(icon=get_icon('refresh'), icon_size='1.8em', align="end")
        pn.bind(self.reset_table, self.reset_btn, watch=True)

        self.table_widgets = pn.Row(self.showselected_btn, self.selectrelated_btn, self.threshold_input,
                            self.search_input, self.searchcol_select, self.search_btn, self.reset_btn,
                            sizing_mode='stretch_width')
        self.table_pane = pn.Column(self.meta_pane,self.table_widgets,sizing_mode='stretch_both')

        #selected table
        self.selected_table = pnw.Tabulator(show_index=False,disabled=True,page_size=100,
                                    frozen_columns=['sample'],stylesheets=[stylesheet],
                                    sizing_mode='stretch_both')
        downloadselection_btn = pnw.FileDownload(callback=pn.bind(self.selection_file),
                                            filename='selection.csv', button_type='success')
        self.selected_pane = pn.Column(self.selected_table,pn.Row(downloadselection_btn))
        #moves table
        self.moves_pane = pnw.Tabulator(show_index=False,disabled=True,page_size=100,
                                    frozen_columns=['tag'],stylesheets=[stylesheet],
                                    sizing_mode='stretch_both')
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

        self.tree_pane = pn.Column(sizing_mode='stretch_both')
        self.mst_pane = pn.Column(sizing_mode='stretch_both')
        self.sdist_pane = pn.pane.DataFrame(sizing_mode='stretch_both')
        #self.sdist_pane = pnw.Tabulator(show_index=True,disabled=True,page_size=100,layout='fit_data_stretch',
        #                                stylesheets=[dm_css],
        #                                sizing_mode='stretch_both')
        #details
        self.details_pane = pn.pane.DataFrame(sizing_mode='stretch_both')

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
        #self.loadlpis_btn = pnw.Button(name='Load LPIS', button_type='primary')
        #pn.bind(self.load_lpis, self.loadlpis_btn, watch=True)
        #card3 = pn.Column('## LPIS', self.loadlpis_btn, width=340, styles=card_style)

        #settings
        self.markersize_input = pnw.IntSlider(name='marker size', value=10, start=0, end=80,width=w)
        self.edgewidth_input = pnw.FloatSlider(name='marker edge width', value=0.8, start=0, end=4, step=.1,width=w)
        self.labelsize_input = pnw.IntSlider(name='label size', value=18, start=6, end=80,width=w)
        self.legendsize_input = pnw.IntSlider(name='legend size', value=12, start=6, end=40,width=w)
        self.hexbins_input = pnw.IntSlider(name='hex bins', value=10, start=5, end=100,width=w)
        self.scalebar_input = pnw.Checkbox(name='show scalebar',value=True)
        self.randseed_input = pnw.IntSlider(name='random color seed', value=647, start=1, end=5000,width=w)
        self.tipsize_input = pnw.IntSlider(name='tree tip size', value=8, start=1, end=25,width=w)
        self.tiplabelsize_input = pnw.IntSlider(name='tip label font size', value=9, start=6, end=20,width=w)
        #self.showtiplabels_input = pnw.Checkbox(name='show tree tip labels',value=True)
        save_js = """
        var settings = {
            marker_size: document.querySelector('[id*="marker size"]').value,
            marker_edge_width: document.querySelector('[id*="marker edge width"]').value,
            label_size: document.querySelector('[id*="label size"]').value,
            legend_size: document.querySelector('[id*="legend size"]').value,
        };
        localStorage.setItem('tracebtb_settings', JSON.stringify(settings));
        """
        savesettings_btn = pnw.Button(name="Save Settings",button_type='primary')
        savesettings_btn.js_on_click(code=save_js)

        card4 = pn.Row(pn.Column('## Settings', self.markersize_input,self.edgewidth_input,self.labelsize_input,
                          self.legendsize_input, self.hexbins_input, self.scalebar_input, self.randseed_input),
                          pn.Column(self.tipsize_input,self.tiplabelsize_input, savesettings_btn),
                          width=340, styles=card_style)

        #reporting
        self.doreport_btn = pnw.Button(name='Generate', button_type='primary')
        self.savereport_btn = pnw.FileDownload(file=report_file, button_type='success', auto=False,
                                        embed=False, name="Download Report")
        card5 = pn.Column('## Reporting', self.doreport_btn, self.savereport_btn, width=340, styles=card_style)
        pn.bind(self.create_report, self.doreport_btn, watch=True)

        utils_pane = pn.FlexBox(*[card1,card2, card4, card5], flex_direction='column', min_height=400,
                                styles={'margin': '10px'}, sizing_mode='stretch_both')

        #nav buttons
        self.prev_btn = pnw.Button(icon=get_icon('arrow-left'), description='back', icon_size='1.6em',width=60)
        self.next_btn = pnw.Button(icon=get_icon('arrow-right'), description='forward', icon_size='1.6em',width=60)
        pn.bind(self.back, self.prev_btn, watch=True)
        pn.bind(self.forward, self.next_btn, watch=True)
        nav_pane = pn.Row(self.prev_btn,self.next_btn)

        #search and cluster widgets
        search_widgets = self.search_widgets()
        group_widgets = self.group_widgets()
        widgets1 = pn.Tabs(('search',search_widgets),('groups',group_widgets))

        #options
        self.colorby_input = pnw.Select(name='color by',options=cols,value='short_name',width=w)
        self.cmap_input = pnw.Select(name='colormap',options=colormaps,value='Set1',width=w)
        self.tiplabel_input = pnw.Select(name='tip label',options=list(self.meta.columns),value='sample',width=w)
        self.provider_input = pnw.Select(name='provider',options=['']+bokeh_plot.providers,value='CartoDB Positron',width=w)
        self.pointstyle_input = pnw.Select(name='points display',options=['default','pie'],value='default',width=w)
        widgets2 = pn.Column(pn.WidgetBox(nav_pane,
                                        self.colorby_input,self.cmap_input,self.tiplabel_input,
                                        self.provider_input,self.pointstyle_input),width=w+30)
        #button toolbar
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
        #self.related_col_input = pnw.Select(options=cols,value='snp7',width=90)

        filters = pn.Row(self.findrelated_btn,self.threshold_input,self.timeslider,self.clustersizeslider,self.info_pane)

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
        kinds = ['count','bar','strip']
        self.catkind_input = pnw.Select(name='kind',options=kinds,
                                value='count',width=w)
        self.cathue_input = pnw.Select(name='hue',options=cols,value='County',width=w)
        self.catcmap_input = pnw.Select(name='colormap',options=colormaps,value='Set1',width=w)
        catupdate_btn = pnw.Button(name='Update', button_type='primary',align='end')
        pn.bind(self.update_catplot, catupdate_btn, watch=True)
        self.analysis_pane1 = pn.Column(pn.Row(self.catx_input,self.caty_input,self.cathue_input,
                                               self.catkind_input,self.catcmap_input,catupdate_btn),
                                    self.catplot_pane)

        styles={ "margin-top": "10px", "font-size": "15px"}
        self.about_pane = pn.pane.Markdown('',styles=styles)
        self.about()

        app = pn.Column(
                    pn.Row(pn.Column(widgets1,widgets2),
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
                                            pn.Column(self.timeline_pane), width=500)),
                                    ('Tree',self.tree_pane),('MST',self.mst_pane),('SNPdist',self.sdist_pane),
                                    ('Details',self.details_pane),
                                    dynamic=True,width=500),
                max_width=2600,min_height=600
        ))

        #widgets = get_bokeh_widgets(app)
        #print (widgets)
        # Attach JS callbacks to save options
        '''save_callback = CustomJS(args=dict(widgets=widgets), code=save_js_code)
        load_callback = CustomJS(args=dict(widgets=widgets), code=load_js_code)

        for w in widgets.values():
            if hasattr(w, 'js_on_change'):
                print (w)
                if hasattr(w, "value"):
                    w.js_on_change("value", save_callback)
                elif hasattr(w, "active"):
                    w.js_on_change("active", save_callback)

        # ✅ Create a dummy Div to attach the startup event
        trigger = Div()
        trigger.js_on_event("document_ready", load_callback)
        #curdoc().add_root(pn.column(trigger, *widgets.values()))'''

        app.sizing_mode='stretch_both'
        self.meta_pane.value = self.meta.fillna('')
        self.update_groups()
        #select random cluster
        #self.selected = self.meta.sample(4).copy()
        strain = self.meta.sample(1).iloc[0].strain_name
        self.selected = self.meta[self.meta.strain_name==strain]
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
        labelcol = self.tiplabel_input.value

        if sub is None:
            sub = self.selected
        #print (sub.columns)
        if 'snp3' not in sub.columns:
            sub = self.update_clusters(sub)

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

        if mov is not None:# and self.moves_btn.value is True:
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
            p.on_event(Tap, self.point_selected)
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
            shb = tools.shared_borders(sp, self.lpis)
            shb['color'] = shb.apply(tools.random_grayscale_color, 1)
            bokeh_plot.plot_lpis(shb, p, line_width=1)

        #mov = tools.get_moves_bytag(sub, self.moves, self.lpis_cent)

        if self.moves_btn.value is True:
            bokeh_plot.plot_moves(p, mov, self.lpis_cent)

        if mov is not None:
            self.moves_pane.value = mov.reset_index().drop(columns=['geometry'])
        else:
            self.moves_pane.value = pd.DataFrame()

        self.plot_pane.object = p

        scols = ['sample','Animal_ID','Year','HERD_NO','snp3','snp5','snp7','lineage','short_name',
                 'County','Region','SB','last_move','last_move_type']
        self.selected_table.value = sub[scols].fillna('')

        def highlight(x):
            color = self.speciescolors[x]
            return 'background-color: %s' % color

        self.selected_table.style.apply(highlight, subset=['Species'], axis=1)

        # Timeline plot
        tlc = self.timelinecolor_input.value
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
            fs = f'{self.tiplabelsize_input.value}pt'
            ts = self.tipsize_input.value
            self.update_tree(sub=sub, col=col, tip_size=ts, font_size=fs, labelcol=labelcol)
            self.update_mst(sub=sub, node_size=ms)
            self.update_snpdist(sub=sub)

        # Update summaries
        self.update_herd_summary()
        self.update_group_summary()
        return

    def update_view(self, event):

        p = self.plot_pane.object
        cur_x = p.x_range.start, p.x_range.end
        cur_y = p.y_range.start, p.y_range.end
        self.update(event)
        p = self.plot_pane.object
        p.x_range.start, p.x_range.end = cur_x
        p.y_range.start, p.y_range.end = cur_y

        #p.x_range = Range1d(cur_x[0], cur_x[1])
        #p.y_range = Range1d(cur_y[0], cur_y[1])
        return

    def point_selected(self, event):
        """Map click callback"""

        p = event.model
        r = p.select({"name": "points"})[0]
        source = r.data_source
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

    def get_map_datasource(self):
        """Get point data source"""

        p = self.plot_pane.object
        r = p.select({"name": "points"})[0]
        return r.data_source

    def tip_selected(self, event):
        """Point click callback"""

        #get map source
        map_source = self.get_map_datasource()
        p = event.model
        print (p)
        r = p.select({"name": "tree_tips"})[0]
        source = r.data_source
        indices = source.selected.indices
        if len(indices)==0:
            return
        idx = indices[0]
        name = source.data['sample'][idx]
        #select point on map
        geojson_features = json.loads(map_source.geojson)['features']
        #print (geojson_features)
        name_to_index = {feature["properties"]['sample']: i for i, feature in enumerate(geojson_features)}
        #print (name_to_index)

        match_idx = name_to_index.get(name, None)
        map_source.selected.indices = [match_idx]

        #set details pane
        data = self.meta.loc[name]
        #print (data)
        self.details_pane.object = pd.DataFrame(data)
        return

    def update_herd_summary(self, event=None):
        """herd summary"""

        h = tools.herd_summary(self.selected, self.moves, self.snpdist)
        self.herds_table.value = h
        return

    def update_group_summary(self, event=None):

        col = self.colorby_input.value
        cl = tools.cluster_summary(self.selected, col, 5, self.snpdist)
        self.clusters_table.value = cl
        return

    def update_overview(self, event=None, sub=None):

        fig = plot_overview(self.selected)
        self.overview_pane.object = fig
        plt.close(fig)
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

        sub = self.selected
        idx = list(sub.index)
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

    def selection_file(self):
        """Return selection table for export"""

        df = self.selected_table.value
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

    def update_clusters(self, sub):
        """Get clusters within the selection - useful for within clade divisions"""

        return tools.add_clusters(sub, self.snpdist, linkage='average')

    def update_mst(self, event=None, sub=None, **kwargs):
        """Update mst"""

        idx = sub.index
        dm = self.snpdist.loc[idx,idx]
        p = bokeh_plot.plot_mst(dm, sub, **kwargs)
        #self.mst_pane.object = p
        self.mst_pane.objects.clear()
        self.mst_pane.append(pn.pane.Bokeh(p))
        return

    def update_snpdist(self, event=None, sub=None):
        """Update snpdist pane"""

        idx = sub.index
        dm = self.snpdist.loc[idx,idx]
        if len(idx)<60:
            self.sdist_pane.object = dm
        else:
            self.sdist_pane.object = None
        print (dm.to_numpy().flatten().max())
        return

    def do_search(self, event=None):
        """Search main table"""

        query = self.search_input.value
        col = self.searchcol_select.value
        #found = self.meta[self.meta[col]==query]
        found = self.meta[self.meta[col].astype(str).str.contains(query, case=False, na=False)]
        print (found)
        self.meta_pane.value = found.drop(columns=['geometry']).fillna('')
        return

    def reset_table(self, event=None):
        """Reset main table"""

        self.meta_pane.value = self.meta.fillna('')
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
        self.add_to_history()
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

    def update_catplot(self, event=None, sub=None):
        """Update categorical plot"""

        if sub is None:
            sub = self.selected
        x = self.catx_input.value
        y = self.caty_input.value
        hue = self.cathue_input.value
        kind = self.catkind_input.value
        cmap = self.catcmap_input.value
        row = None

        if kind== 'count':
            cg = sns.catplot(data=sub, x=x, y=y, hue=hue, kind=kind, aspect=2,
                            palette=cmap)
        else:
            cg = sns.catplot(data=sub, x=x, y=y, hue=hue, kind=kind, aspect=2,
                            palette=cmap, dodge=True)
        self.catplot_pane.object = cg.fig
        plt.close(cg.fig)
        #print (fig)
        return

    def save_settings(self, event=None):
        """Save user interface settings"""

        ignored_types = [pnw.Tabulator]
        settings_dict = {
            name: widget.value
            for name in dir(self)
            if isinstance((widget := getattr(self, name, None)), pnw.Widget)
            and hasattr(widget, "value")
            and not isinstance(widget, tuple(ignored_types))  # Exclude certain types
        }
        print (settings_dict)
        return

    def restore_settings(self, settings_dict):

        for name, value in settings_dict.items():
            widget = getattr(self, name, None)
            if isinstance(widget, pnw.Widget) and hasattr(widget, "value"):
                widget.value = value  # Restore the value

class QueryDashboard(FullDashboard):
    def __init__(self, **kwargs):
        super(QueryDashboard, self).__init__(**kwargs)
        return

    def setup_widgets(self):

        w=140
        icsize = '1.5em'
        cols = self.cols
        self.recents = []
        #search widgets
        search_widgets = self.search_widgets()
        #tools
        self.selectregion_btn = pnw.Button(icon=get_icon('plot-region'), description='select in region', icon_size=icsize)
        pn.bind(self.select_region, self.selectregion_btn, watch=True)
        self.parcels_btn = pnw.Toggle(icon=get_icon('parcels'), value=True, icon_size=icsize)
        #self.moves_btn = pnw.Toggle(icon=get_icon('plot-moves'), icon_size=icsize)
        #self.neighbours_btn = pnw.Toggle(icon=get_icon('neighbours'), icon_size=icsize)
        self.parcellabel_btn = pnw.Toggle(icon=get_icon('parcel-label'), icon_size=icsize)
        toolbar = pn.GridBox(self.selectregion_btn,self.parcels_btn,self.parcellabel_btn,
                             ncols=3,width=w)
        self.colorby_input = pnw.Select(name='color by',options=cols,value='snp5',width=w)
        self.colorby_input.param.watch(self.update, 'value')
        self.plot_pane = pn.pane.Bokeh()
        self.selected_table = pnw.Tabulator(disabled=True,page_size=100,
                                          pagination=None, sizing_mode='stretch_both',
                                          stylesheets=[stylesheet])
        downloadselection_btn = pnw.FileDownload(callback=pn.bind(self.selection_file),
                                            filename='selection.csv', button_type='success')
        self.selected_pane = pn.Column(self.selected_table,pn.Row(downloadselection_btn))
        self.tree_pane = pn.Column(sizing_mode='stretch_both')
        self.mst_pane = pn.Column(sizing_mode='stretch_both')
        self.details_pane = pn.pane.DataFrame(sizing_mode='stretch_both')
        self.info_pane = pn.pane.Markdown('', styles={'color': "red"})

        app = pn.Row(pn.Column(search_widgets,toolbar,self.colorby_input,self.info_pane),
                        pn.Tabs(('Map',pn.Column(self.plot_pane,sizing_mode='stretch_both')),
                                ('Selected',pn.Column(self.selected_pane,sizing_mode='stretch_both')),
                                dynamic=True,
                                sizing_mode='stretch_both'),
                     pn.Column(self.tree_pane, self.mst_pane, width=600))

        return app

    def update(self, event=None, sub=None):
        """
        Update selection - does all the display.
        """

        provider = "CartoDB Positron"
        cmap = 'Dark2'
        self.col = col = self.colorby_input.value
        self.ms = ms = 15
        lw = 2
        legend = True

        if sub is None:
            sub = self.selected

        if len(sub[col].unique()) > 20:
            cmap = None
        elif len(sub[col].unique()) > 10:
            cmap = 'tab20'
        seed = 50
        sub['color'], cm1 = tools.get_color_mapping(sub, col, cmap, seed=seed)
        self.selected = sub

        # Get the herds and moves
        herds = list(sub.HERD_NO)
        mov = tools.get_moves_bytag(sub, self.moves, self.lpis_cent)

        #if mov is not None and self.moves_btn.value is True:
        #    herds.extend(mov.move_to)

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
        labelsize = 12
        legsize = 15
        p = bokeh_plot.plot_selection(sub, pcls, provider=provider, ms=ms, lw=lw,
                                        col=col, legend=legend, labels=labels,
                                        legend_fontsize=legsize, label_fontsize=labelsize,
                                        scalebar=True)
        p.on_event(Tap, self.point_selected)
        self.plot_pane.object = p

        #tree
        ns = (120/len(sub))+6
        labelcol = self.tiplabel_input.value
        self.update_tree(sub=sub, col=col, tip_size=ns, labelcol=labelcol, font_size='9pt')
        #mst
        self.update_mst(sub=sub, node_size=ns, labels=False)
        self.info_pane.object = f'**{len(sub)} samples**'
        scols = ['sample','Animal_ID','Year','HERD_NO','snp3','snp5','snp7','short_name','County','SB']
        self.selected_table.value = sub[scols]
        return

    def quick_search(self, event=None, query=None):
        """Search query string, allows comma separated list"""

        if query == None:
            query = self.search_input.value
        if len(query)<=1:
            return
        query_list = [q.strip().lower() for q in query.split(",")]
        found = self.meta[self.meta.map(lambda x: any(q in str(x).lower() for q in query_list)).any(axis=1)]
        if len(found)>1000 or len(found)==0:
            return

        col='snp5'
        #related animals
        cl = list(found[col].unique())
        found = self.meta[(self.meta[col].isin(cl))].copy()
        #same herd?
        self.update(sub=found)
        self.add_to_history()
        self.add_to_recent(query)
        return

class HerdSelectionDashboard(Dashboard):

    def __init__(self, **kwargs):
        super(HerdSelectionDashboard, self).__init__(**kwargs)
        self.get_test_stats()
        self.feedlots = kwargs['feedlots']
        self.iregrid = kwargs['ireland_grid']
        self.herd = None
        movement.BASEDIR = self.settings['moves_lake']
        return

    def setup_widgets(self):
        """Add widgets"""

        w=140
        styles={"margin": "10px", "font-size": "15px"}
        self.plot_pane = pn.pane.Bokeh()
        self.seqselection_pane = pnw.Tabulator(show_index=False,disabled=True)
        self.seqselection_pane.value = pd.DataFrame()
        self.tree_pane = pn.Column(sizing_mode='stretch_both')
        self.fragments_pane = pn.pane.Matplotlib(sizing_mode='stretch_both')
        self.grid_pane = pn.pane.Matplotlib(sizing_mode='stretch_both')
        self.samples_pane = pnw.Tabulator(show_index=False,disabled=True,page_size=100,
                                          stylesheets=[stylesheet], sizing_mode='stretch_both')
        self.moves_pane = pnw.Tabulator(show_index=False,disabled=True,page_size=100,
                                    frozen_columns=['tag'],stylesheets=[stylesheet],
                                    sizing_mode='stretch_both')
        self.testingplot_pane = pn.pane.Bokeh(height=200)
        self.movesplot_pane = pn.pane.Bokeh(height=200)
        self.herdinfo_pane = pn.pane.DataFrame(stylesheets=[df_stylesheet], sizing_mode='stretch_both')
        self.indicator = pn.indicators.Number(
                    name='Priority',
                    value=0, format="{value}",
                    colors=[(0, "red"), (1, "green")])
        mstyles={"margin": "10px", "font-size": "17px", "color": 'red'}
        self.seq_priority_pane = pn.pane.Markdown('',width=480,height=20,styles=mstyles)
        self.pathways_pane = pn.pane.DataFrame(stylesheets=[df_stylesheet], sizing_mode='stretch_both')
        search_widgets = self.search_widgets(4)
        herds_with_samples = self.get_sampled_herds()
        self.sampled_herds_select = pnw.Select(name='Sampled Herds',value='',
                                               options=herds_with_samples,width=w,size=4)
        self.sampled_herds_select.param.watch(self.sampled_herd_search,'value')

        sim_btn = pnw.Button(name='Random Herd',button_type='primary',width=w)
        pn.bind(self.random_herd, sim_btn, watch=True)
        refresh_btn = pnw.Button(name='Refresh',width=w,button_type='success')
        pn.bind(self.update, refresh_btn, watch=True)
        #sendquery_btn = pnw.Button(name='Samples Query',width=w)
        #pn.bind(self.send_to_query_dashboard, sendquery_btn, watch=True)
        self.dist_input = pnw.IntSlider(name='Dist Threshold',width=w,value=1000,
                                        start=100,end=10000,step=100)
        opts = ['within any parcel','contiguous parcels']
        self.dist_method = pnw.Select(name='Dist Method',value='within any parcel',
                                      options=opts,width=w)
        self.provider_input = pnw.Select(name='Provider',options=['']+bokeh_plot.providers,
                                         value='CartoDB Positron',width=w)
        self.provider_input.param.watch(self.set_provider, 'value')
        colorby = ['SPH_HERD_N','Size','sr_total','COUNT','AREA','base_year']
        self.colorby_input = pnw.Select(name='Color Parcels by',value='SPH_HERD_N',options=colorby,width=w)

        widgets = pn.Column(search_widgets,self.sampled_herds_select,sim_btn,refresh_btn,
                            self.dist_method,self.dist_input,self.provider_input,self.colorby_input)
        self.about_pane = pn.pane.Markdown('',styles=styles)
        self.about()
        #self.rules_pane = pn.pane.Markdown('',styles=styles, sizing_mode='stretch_both')
        #doc_path = os.path.join(module_path,'docs','selection_framework.md')
        #rules_txt = open(doc_path,'r').read()
        #self.rules_pane.object = rules_txt
        #from . import agent
        #chat_pane = agent.ai_chat_widget(self.meta)

        app = pn.Row(widgets,
                  #pn.Column(self.plot_pane,sizing_mode='stretch_both'),
                  pn.Column(pn.Tabs(('Map',pn.Column(self.plot_pane,sizing_mode='stretch_both')),
                                    #('Sequence Selection',self.seqselection_pane),
                                    #('Rules',self.rules_pane),
                                    ('Samples',self.samples_pane),
                                    ('Moves',self.moves_pane),
                                    ('Pathways',self.pathways_pane),
                                    #('Chat',chat_pane),
                                    ('About',self.about_pane),
                             dynamic=True,
                             sizing_mode='stretch_both')),
                  pn.Column(pn.Tabs(('Herd info',self.herdinfo_pane),
                                    ('Metrics',pn.Column(self.testingplot_pane, self.movesplot_pane)),
                                    ('Tree',self.tree_pane),
                                    ('Fragments',self.fragments_pane),
                                    ('Grid',self.grid_pane)),
                            self.seq_priority_pane,width=500),
                    sizing_mode='stretch_both')
        app.sizing_mode='stretch_both'
        return app

    def get_sampled_herds(self):

        b = self.meta[self.meta.Species=='Bovine']
        return list(b.HERD_NO.dropna().sort_values().unique())

    def random_herd(self, event=None, herd=None):
        """Select herd at random"""

        herd = self.random_breakdown_herd()
        self.update(herd=herd)
        self.add_to_recent(herd)
        return

    def herd_info(self):
        """Get herd metrics"""

        herd_no = self.herd
        data = tools.herd_metrics(herd_no, self.meta, self.moves, self.testing,
                               self.feedlots, self.lpis, self.lpis_cent, cell)
        s = pd.Series(data)
        self.herdinfo_pane.object = pd.DataFrame(s,columns=['value'])
        return

    def update(self, event=None, herd=None):
        """Update display"""

        if herd is None:
            herd = self.herd
        else:
            self.herd = herd
        if herd is None:
            return
        #print (herd)
        lpis = self.lpis
        meta = self.meta
        dist = self.dist_input.value
        parcelcol = self.colorby_input.value
        #herd context gets everything about the herd in a dict
        hc = tools.get_herd_context(self.herd, self.meta, self.moves, self.testing,
                               self.feedlots, self.lpis, self.lpis_cent, self.iregrid, dist=dist)
        hdata = hc['data']
        s = pd.Series(hc['metrics'])
        self.herdinfo_pane.object = pd.DataFrame(s,columns=['value'])
        sub = hdata['herd_isolates']

        #get neighbours and others
        #self.parcels = pcl = lpis[lpis.SPH_HERD_N==herd].copy()
        #pcl['HERD_NO'] = pcl.SPH_HERD_N
        self.parcels = pcl = hdata['herd_parcels']
        #self.cont_parcels = tools.shared_borders(pcl, lpis)
        self.neighbours = hdata['neighbour_parcels']
        self.cont_parcels = hdata['contiguous_parcels']

        if self.dist_method.value=='within any parcel':
            nb = self.neighbours
        else:
            nb = self.cont_parcels
        #total pos tests over time
        nb['sr_total'] = nb.apply(lambda x: tools.get_testing_total(x,self.sr),1)
        #print (nb.columns)
        #nearby badgers
        bdg = hdata['near_badger']

        #then get any known strains present in area, including in herd itself
        qry = list(nb.SPH_HERD_N) + list(bdg.HERD_NO) + [herd]
        found = meta[meta.HERD_NO.isin(qry)]
        found['color'],cm = tools.get_color_mapping(found, 'short_name', cmap='Set1')
        self.found = found
        p = self.plot_neighbours(pcl,nb,column=parcelcol)
        bokeh_plot.plot_selection(found, legend=True, ms=20, lw=2, p=p)
        p.title = herd
        p.title.text_font_size = '20pt'
        self.plot_pane.object = p
        self.plot_herd_testing(herd)
        self.plot_movements_summary(herd)

        #nearest sampled herd
        self.nearest = tools.find_nearest_point(pcl.iloc[0].geometry, found)
        #tree
        self.update_tree(sub=found, col='short_name', labelcol='HERD_NO')
        #plot fragments
        mp = pcl.iloc[0]
        G,ns = tools.fragments_to_graph(mp)
        fig = tools.plot_herd_graph(G, gdf=found, title=herd)
        self.fragments_pane.object = fig
        #grid cells plot
        cell = hdata['grid_cell']
        fig = tools.plot_grid_cells(cell, gdf=self.meta, parcels=pcl, neighbours=nb)
        self.grid_pane.object = fig
        #samples
        cols = ['sample','HERD_NO','Animal_ID','Species','short_name','last_move']
        self.samples_pane.value = found[cols]
        #moves
        mov = hdata['isolate_moves']
        if mov is not None:
            self.moves_pane.value = mov.reset_index().drop(columns=['geometry'])

        priority = tools.get_sequencing_priority(hc)
        self.seq_priority_pane.object = priority

        pathway_scores = source_attribution.run_pathways(sub, self.meta,
                                                self.moves, self.sr, self.feedlots,
                                                self.lpis, self.lpis_cent, self.iregrid)
        self.pathways_pane.object = pathway_scores
        return

    def send_to_query_dashboard(self, event=None):
        """Send herd query to samples query dashboard
        - will need to use pn.state.location?
        """

        dash = self.parent.query_dashboard
        herds = ','.join(list(self.found.HERD_NO))
        dash.quick_search(query=herds)
        return

    def quick_search(self, event=None, query=None):
        """Search lpis for a single herd"""

        if query == None:
            query = self.search_input.value
        if len(query)<=1:
            return

        found = self.lpis[self.lpis.SPH_HERD_N==query]
        if len(found)==0:
            return
        self.update(herd=query)
        self.add_to_recent(query)
        return

    def get_test_stats(self):
        """Summary stats from testing file"""

        te = self.testing
        sr = te.filter(regex='^Sr')
        lp = te.filter(regex='^Lp')
        hs = te.filter(regex='^Hs')
        sr_total=sr.iloc[:,-6:].sum(1)
        lp_total=lp.iloc[:,-6:].sum(1)
        hs_med=hs.iloc[:,-6:].median(1)
        teststats = pd.concat([sr_total,lp_total,hs_med],axis=1)
        teststats.columns=['Sr','Lp','Size']
        self.lpis = self.lpis.merge(teststats,left_on='SPH_HERD_N',right_index=True, how='left')
        return

    def plot_herd_testing(self, herd):
        """Get testing plot"""

        te = self.testing
        sr = te.filter(regex='^Sr')
        df = sr.loc[herd]
        x = list(df.index)
        vals = df.values
        p = figure(x_range=x, height=250, title="Std Reactors",
                toolbar_location=None, tools="")
        p.vbar(x=x, top=vals, width=0.9)
        self.testingplot_pane.object = p
        return p

    def plot_movements_summary(self, herd):
        """Get moves summary"""

        m = movement.query_herd_moves_in(herd, date(2012, 1, 1), date(2022, 1, 1))
        if m is None:
            return figure()
        df = m.groupby('year')['tag'].count()
        x=list(df.index)
        vals = df.values
        p = figure(height=250, title="Moves",
                toolbar_location=None, tools="")
        p.vbar(x=x, top=vals, width=0.9)
        self.movesplot_pane.object = p
        return p

    def random_breakdown_herd(self):

        te = self.testing
        sr = te.filter(regex='^Sr').copy()
        sr['total'] = sr.iloc[:,-5:].sum(1)
        #sr=sr.sort_values('sum')
        x=sr[sr.total>5]
        return x.sample(1).index[0]

    def plot_breakdown_history(self, herds):
        """Track multiple herds tests as line graph?"""

        return

class TestingDashboard(FullDashboard):
    """Testing app is a wrapper for various test dashboards"""
    def __init__(self, **kwargs):

        self.query_dashboard = QueryDashboard(parent=self, **kwargs)
        self.herdselect_dashboard = HerdSelectionDashboard(parent=self, **kwargs)
        self.advanced_dashboard = FullDashboard(parent=self, **kwargs)
        self.advanced_dashboard.lpis = self.herdselect_dashboard.lpis
        super(TestingDashboard, self).__init__(**kwargs)
        return

    def setup_widgets(self):
        advanced_pane = self.advanced_dashboard.show()
        query_pane = self.query_dashboard.show()
        herdselect_pane = self.herdselect_dashboard.show()
        app = pn.Row(
            pn.Tabs(('Herd Selection', herdselect_pane),
                    ('Sample Query', query_pane),
                    ('Advanced', advanced_pane)),
                    #('Moves',herds_pane)),
            max_width=2600,min_height=600)

        app.sizing_mode='stretch_both'
        return app

class MovesDashboard(Dashboard):
    """Moves dedicated dashboard"""

    def __init__(self, **kwargs):
        super(MovesDashboard, self).__init__(**kwargs)
        self.feedlots = kwargs['feedlots']
        self.iregrid = kwargs['ireland_grid']
        movement.BASEDIR = self.settings['moves_lake']
        p = bokeh_plot.init_figure()
        self.networkplot_pane.object = p
        return

    def setup_widgets(self):
        """Create widgets"""

        w=140
        self.recents = []
        self.herd_pane = pn.pane.Bokeh()
        #self.networkplot_pane = pn.pane.Matplotlib(width=600,height=600)
        self.networkplot_pane = pn.pane.Bokeh(height=600)
        self.summaryplot_pane = pn.pane.Matplotlib(width=600,height=300)
        self.summaryplot_pane2 = pn.pane.Matplotlib(width=600,height=300)
        search_widgets = self.search_widgets()
        rnd_btn = pnw.Button(name='Random Herd',button_type='primary',width=w)
        refresh_btn = pnw.Button(name='Refresh',button_type='primary',width=w)
        pn.bind(self.random_herd, rnd_btn, watch=True)
        pn.bind(self.update, refresh_btn, watch=True)

        #button toolbar
        self.parcels_btn = pnw.Toggle(icon=get_icon('parcels'), icon_size=icsize)
        self.parcellabel_btn = pnw.Toggle(icon=get_icon('parcel-label'), icon_size=icsize)
        toolbar = pn.GridBox(self.parcels_btn,self.parcellabel_btn,width=w,ncols=3)
        self.parcels_btn.param.watch(self.update, 'value')
        self.parcellabel_btn.param.watch(self.update, 'value')

        self.date_range_slider = pn.widgets.DateRangeSlider(
            name='Date Range',
            start=datetime(2008, 1, 1), end=datetime(2022, 1, 1),
            value=(datetime(2020, 1, 1), datetime(2022, 1, 1)),
            step=30,format='%Y-%m',width=200
        )
        self.info_pane = pn.pane.Markdown('', styles={'color': "red",'font_size':'18px'})
        widgets = pn.WidgetBox(rnd_btn, refresh_btn, toolbar, width=w+20)
        self.main_moves_table = pnw.Tabulator(show_index=False,disabled=True,page_size=200,
                                    frozen_columns=['tag'],stylesheets=[stylesheet],
                                    sizing_mode='stretch_both')
        widgets2 = pn.Row(self.date_range_slider,self.info_pane)
        app = pn.Row(pn.Column(search_widgets,widgets),
                    pn.Row(
                        pn.Tabs(('Network',pn.Row(
                            pn.Column(self.networkplot_pane,widgets2,sizing_mode='stretch_both'),
                            self.herd_pane,sizing_mode='stretch_both')),
                                ('Moves',self.main_moves_table),
                                dynamic=True,
                                sizing_mode='stretch_both'),
                        pn.Tabs(('Data',pn.Column(self.summaryplot_pane,self.summaryplot_pane2)),
                                sizing_mode='stretch_both'),
                         ),
                max_width=2600,min_height=600)
        app.sizing_mode='stretch_both'
        return app

    def update(self, event=None, moves=None, herd=None):
        """draw individual farm and then moves over time."""

        col='snp7'
        provider = "CartoDB Positron"
        cmap = 'Dark2'
        labelsize = 12
        legsize = 15
        self.ms = ms = 6
        lw = .5
        legend = True
        dist = 100

        if moves is None:
            moves = self.moves
        else:
            self.moves = moves
        if moves is None:
            return

        #filter by dates
        start, end = self.date_range_slider.value
        start = pd.to_datetime(start)
        end = pd.to_datetime(end)
        moves = self.moves[(self.moves.event_date>start) & (self.moves.event_date<=end)]

        #get parcels etc for all herds in moves
        herds = moves.herd.unique()
        pcl = self.lpis[self.lpis.SPH_HERD_N.isin(herds)]
        pcl['color'], cm = tools.get_color_mapping(pcl,'SPH_HERD_N', cmap)
        #p = self.plot_neighbours(pcl,nb)#,column=parcelcol)
        labels = self.parcellabel_btn.value
        if len(pcl)>0 and self.parcels_btn.value is True:
            p = bokeh_plot.plot_lpis(pcl,fill_alpha=0.5, line_width=0.2)
        else:
            p = bokeh_plot.init_figure(provider=provider)
        #moves = tools.get_moves_locations(moves, self.lpis_cent)

        G,pos = movement.create_herd_network(moves, herd, self.lpis_cent)
        #self.plot_herd_network(G, pos)
        self.plot_herd_network_bokeh(G, pos, p)
        p.title = herd
        self.networkplot_pane.object = p

        #sub = self.meta[self.meta.HERD_NO.isin(herds)]

        moves_in = moves[(moves.herd==herd) & (moves.data_type!='Birth')]
        births = moves[moves.data_type=='Birth']
        self.info_pane.object = f"moves in: {len(moves_in)}; births: {len(births)}"
        self.main_moves_table.object = moves_in
        self.plot_summaries(moves)
        return

    def quick_search(self, event=None, query=None):
        """Search lpis for a single herd"""

        if query == None:
            query = self.search_input.value
        if len(query)<=1:
            return

        #start, end = self.date_range_slider.value
        start = date(2008, 1, 1)
        end = date(2022, 1, 1)
        moves = movement.query_herd_moves_all(query, start, end)
        if len(moves)==0 or moves is None:
            return
        self.update(moves=moves, herd=query)
        self.add_to_recent(query)
        return

    def random_herd(self, event=None, herd=None):
        """Select herd at random"""

        herd = self.lpis.sample(1).iloc[0].SPH_HERD_N
        #start, end = self.date_range_slider.value
        start = date(2008, 1, 1)
        end = date(2022, 1, 1)
        st=time.time()
        moves = movement.query_herd_moves_all(herd, start, end)
        print (f'got herd moves in {time.time()-st}')
        print (f'{len(moves)} moves found')
        self.update(moves=moves, herd=herd)
        self.add_to_recent(herd)
        return

    def plot_summaries(self, moves):
        """Get moves summary"""

        m = moves
        if m is None:
            return figure()
        '''df = m.groupby('year')['tag'].count()
        p1 = figure(width=500,height=250, title="Moves by year",
                toolbar_location=None, tools="")
        p1.vbar(x=list(df.index), top=df.values, width=0.9)
        df = m.groupby('data_type')['tag'].count()
        p2 = figure(x_range=list(df.index), height=250, title="Moves by type",
                toolbar_location=None, tools="")
        p2.vbar(x=list(df.index), top=df.values, width=0.9)
        #df = m.groupby('event_date')['tag'].count().reset_index()
        #p3 = figure(height=250, title="Moves by date",  x_axis_type="datetime",
        #        toolbar_location=None, tools="")
        #p3.line(df.event_date, df.tag, width=0.9, color='navy')
        p = column(p1, p2)
        self.summaryplot_pane.object = p'''
        sns.set_context('poster')
        cg = sns.catplot(m,x='year',hue='data_type',kind='count',aspect=2)
        cg.ax.set_title('Moves by year')
        self.summaryplot_pane.object = cg.fig
        plt.tight_layout()
        plt.close(cg.fig)
        coords_df = self.lpis_cent
        df_vectors = movement.get_movement_vectors(m, coords_df)
        movement.categorize_moves(df_vectors)
        cg=sns.catplot(df_vectors, x='move_category', kind='count',aspect=2)
        self.summaryplot_pane2.object = cg.fig
        plt.close(cg.fig)
        return

    def plot_herd_network(self, G, pos, counties=None):
        """Plot herd network"""

        import networkx as nx
        fig,ax = plt.subplots(1,1,figsize=(10, 10))
        nx.draw_networkx_nodes(G, pos, node_size=50, node_color='blue', alpha=0.7, ax=ax)
        nx.draw_networkx_edges(G, pos, arrowstyle='->', arrowsize=10,
                            edge_color='gray', alpha=0.6, width=0.6, ax=ax)
        plt.title("Herd Movement Network")
        plt.axis('off')
        if counties is not None:
            counties.plot(ax=ax)
        self.networkplot_pane.object = fig
        return

    def transform_coords_to_mercator(self, pos, from_epsg=29902):
        """Transform graph positions to mercartor before plotting"""

        from pyproj import Transformer
        transformer = Transformer.from_crs(from_epsg, 3857, always_xy=True)
        mercator_pos = {}
        for herd_id, (x, y) in pos.items():
            m_x, m_y = transformer.transform(x, y)
            mercator_pos[herd_id] = (m_x, m_y)
        return mercator_pos

    def plot_herd_network_bokeh(self, G, pos, p=None):
        """Plot herd network with bokeh and map overlay"""

        from bokeh.models import (BoxZoomTool, HoverTool, MultiLine, Plot, Range1d,
                                ResetTool, Circle, NodesAndLinkedEdges, TapTool)
        from bokeh.plotting import from_networkx
        from bokeh.transform import factor_cmap

        pos = self.transform_coords_to_mercator(pos)
        if p is None:
            p = bokeh_plot.init_figure(provider='CartoDB Positron')
        # Create the Bokeh graph renderer from NetworkX
        graph_renderer = from_networkx(G, pos, scale=1, center=(0, 0))
        mapper = factor_cmap('link', palette=['red', 'gray'], factors=['direct', 'indirect'])
        graph_renderer.node_renderer.glyph = Circle(radius=5, radius_units='screen',
                                                    fill_color=mapper, fill_alpha=0.6, line_color='black')
        graph_renderer.node_renderer.hover_glyph = Circle(radius=8, fill_alpha=0.6, radius_units='screen')
        graph_renderer.edge_renderer.glyph = MultiLine(line_color="#746969", line_alpha=0.6, line_width=1)
        graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color='green', line_width=3)
        # Add Interactive Hover Policy
        graph_renderer.inspection_policy = NodesAndLinkedEdges()

        hover_tool = HoverTool(tooltips=[("Herd", "@herd")])
        p.add_tools(hover_tool)
        p.renderers.append(graph_renderer)
        p.xaxis.visible = False
        p.yaxis.visible = False
        p.toolbar.logo = None
        bokeh_plot.add_scalebar(p)
        return p

    def plot_moves_animation(self, event=None):
        """Show moves over a period"""

        p = self.plot_pane.object
        lines = p.select({'name': 'moves'})
        #p.renderers.remove(lines)
        df = self.related_moves
        start, end = self.date_range_slider.value
        #interval
        n_days = 30
        # Generate a list of dates with n-day intervals
        date_range = pd.date_range(start=start, end=end, freq=f"{n_days}D")
        # Iterate over date ranges
        for start_date in date_range:
            end_date = start_date + pd.Timedelta(days=n_days - 1)
            df_slice = df[(df['move_date'] >= start_date.date()) & (df['move_date'] <= end_date.date())]
            print(f"Data from {start_date} to {end_date}:\n", len(df_slice))
            print (df_slice)
            bokeh_plot.plot_moves(p, df_slice, self.lpis_cent, limit=200, name='related')
        return

class ABMDashboard():
    """ABM dashboard"""

    def setup_widgets(self):
        """Create widgets"""

        w=140
        self.recents = []
        self.plot_pane = pn.pane.Bokeh()

        graph_types = ['default','watts_strogatz','erdos_renyi','barabasi_albert',
                       'powerlaw_cluster','random_geometric']
        farm_types = ['mixed','beef','dairy','suckler']
        cmaps = ['Blues','Reds','Greens','RdBu','coolwarm','summer','winter','icefire','hot','viridis']
        #grid_pane = pn.pane.Matplotlib(plt.Figure(),tight=True,width=900,height=620)
        self.grid_pane = pn.pane.Bokeh(sizing_mode='stretch_both')
        self.plot_pane1 = pn.pane.Matplotlib(plt.Figure(),height=300)
        self.plot_pane2 = pn.pane.Matplotlib(plt.Figure(),height=300)
        self.tree_pane = pn.pane.HTML()
        self.str_pane = pnw.TextAreaInput(disabled=True,height=600,width=400)
        self.df_pane = pnw.Tabulator(show_index=False,disabled=True,height=600,stylesheets=[stylesheet])
        self.df2_pane = pnw.Tabulator(show_index=False,disabled=True,height=600,stylesheets=[stylesheet])

        w=140
        colorby = ['num_infected','perc_infected','herd_size','loc_type','herd_class','strain']
        go_btn = pnw.Button(name='run',width=w,button_type='success')
        stop_btn = pnw.Button(name='stop',width=w,button_type='danger')
        go_btn.param.watch(self.run_model, 'clicks')
        self.inputs = {}
        self.inputs['F'] = pnw.IntSlider(name='farms',value=20,start=5,end=1000,step=1,width=w)
        self.inputs['C'] = pnw.IntSlider(name='cows',value=400,start=10,end=5000,step=10,width=w)
        self.inputs['S'] = pnw.IntSlider(name='setts',value=5,start=1,end=100,step=1,width=w)
        #self.inputs['farmtypes'] = pnw.Select(name='farm types',options=farm_types,width=w)
        self.inputs['cctrans'] = pnw.FloatSlider(name='CC trans',value=1,step=.1,start=0,end=5,width=w)
        self.inputs['bctrans'] = pnw.FloatSlider(name='BC trans',value=1,step=.1,start=0,end=5,width=w)
        self.inputs['mean_stay_time'] = pnw.FloatSlider(name='mean stay time',value=100,step=1,start=5,end=1000,width=w)
        self.inputs['mean_inf_time'] = pnw.FloatSlider(name='mean inf. time',value=60,step=1,start=5,end=600,width=w)
        self.inputs['mean_latency_time'] = pnw.FloatSlider(name='mean latency time',value=100,step=1,start=10,end=600,width=w)
        self.inputs['infected_start'] = pnw.FloatSlider(name='start infected',value=5,step=1,start=1,end=500,width=w)
        self.steps_input = pnw.IntSlider(name='steps',value=10,start=1,end=2000,width=w)
        self.refresh_input = pnw.IntSlider(name='refresh rate',value=1,start=1,end=100,width=w)
        self.delay_input = pnw.FloatSlider(name='step delay',value=0,start=0,end=3,step=.2,width=w)
        graph_input = pnw.Select(name='graph type',options=graph_types,width=w)
        graph_seed_input = pnw.IntInput(name='graph seed',value=10,width=w)
        #seed_input = pnw.Select(name='graph seed',options=['random'],width=w)
        self.colorby_input = pnw.Select(name='color by',value='strain',options=colorby,width=w)
        self.cmap_input = pnw.Select(name='colormap',options=cmaps,width=w)
        self.nodesize_input = pnw.Select(name='node size',value='herd_size',options=colorby[:3],width=w)
        self.labels_input = pnw.Checkbox(name='node labels',value=False,width=w)
        self.progress = pn.indicators.Progress(name='Progress', value=0, width=600, bar_color='primary')

        widgets = pn.Column(pn.Tabs(('model',pn.WidgetBox(*self.inputs.values(),self.steps_input)),
                                    ('options',pn.WidgetBox(graph_input,graph_seed_input,self.colorby_input,
                                                            self.cmap_input,self.nodesize_input,self.labels_input))), width=w+30)

        app = pn.Column(pn.Row(go_btn,self.progress),
                pn.Row(widgets,self.grid_pane,
                 pn.Tabs(('plots',pn.Column(self.plot_pane1,self.plot_pane2)), ('tree',self.tree_pane),
                         ('inf_data',self.df_pane), ('herd_data',self.df2_pane)),
                 sizing_mode='stretch_both',max_width=2600,min_height=600))

        return app

    def run_model(self, event=None):

        from btbabm import models, utils
        def callback(x):
            self.str_pane.value += str(x)+'\n'
        ns = self.nodesize_input.value
        col = self.colorby_input.value
        params = {k:self.inputs[k].value for k in self.inputs}
        print (params)
        steps = self.steps_input.value
        refresh = self.refresh_input.value
        delay = self.delay_input.value

        model = models.FarmPathogenModel(**params, callback=callback)
        self.str_pane.value = ''
        callback(model)
        fig1,ax1 = plt.subplots(1,1,figsize=(15,10))
        #self.grid_pane.object = fig1
        fig2,ax2 = plt.subplots(1,1,figsize=(8,6))
        self.plot_pane1.object = fig2
        fig3,ax3 = plt.subplots(1,1,figsize=(8,6))
        self.plot_pane2.object = fig3
        self.progress.max=steps
        self.progress.value = 0

        showsteps = list(range(1,steps+1,refresh))
        #step through the model and plot at each step
        for i in range(1,steps+1):
            model.step()
            plt.close()
            if i in showsteps:
                #ax1.clear()

                y=model.year
                mov=len(model.get_moves())
                deaths=model.deaths
                total = len(model.get_animals())
                text='day=%s year=%s moves=%s deaths=%s animals=%s' %(i,y,mov,deaths,total)
                utils.plot_grid(model,ax=ax1,pos=model.pos,
                          title=text, colorby=col, cmap=self.cmap_input.value,
                          ns=ns)#, with_labels=labels_input.value)
                p = utils.plot_grid_bokeh(model, title=text)
                self.grid_pane.object = p
                #self.grid_pane.param.trigger('object')

                ax2.clear()
                #s = model.circulating_strains()
                #d=model.get_infected_data()
                #df_pane.value = d
                hd = model.get_herds_data()
                self.df2_pane.value = hd
                #fig2 = utils.plot_inf_data(model)
                #fig2 = utils.plot_by_species(model)
                #self.plot_pane1.object = fig2
                #self.plot_pane1.param.trigger('object')
                df = model.get_column_data()
                ax3.clear()
                df.plot(ax=ax3)
                ax3.set_xlim(0,steps)
                self.plot_pane2.param.trigger('object')
                #html=html_tree(model)
                #self.tree_pane.object = html
                out = model.G.nodes

            self.progress.value += 1
            time.sleep(delay)
        plt.clf()
        return

    def html_tree(self, model):

        result = model.make_phylogeny(removed=True,redundant=False)
        if result==None:
            return '<p>no tree</p>'
        cl = model.get_clades('tree.newick')
        idf = model.get_animal_data(removed=True)
        x=idf.merge(cl,left_on='id',right_on='SequenceName')
        x=x.set_index('SequenceName')
        x.index = x.index.map(str)
        #tre=toytree.tree('tree.newick')
        col='strain'
        #canvas = utils.draw_tree('tree.newick',x,col,tip_labels=False,width=500)

        return

    def set_stop(self, event):
        global stop
        stop = True
        print ('STOP')

