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
from bokeh.models import Range1d
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

speciescolors = {'Bovine':'blue','Badger':'red','Ovine':'green'}
speciesmarkers = {'Bovine':'circle','Badger':'square','Ovine':'diamond',None:'x'}
card_style = {
    'background': '#f9f9f9',
    'border': '1px solid #bcbcbc',
    'padding': '5px',
    'margin': '5px',
    'box-shadow': '4px 4px 4px #bcbcbc'
}
colormaps = ['Paired', 'Dark2', 'Set1', 'Set2', 'Set3',
            'tab10', 'tab20', 'tab20b', 'tab20c']

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

def draw_toytree(treefile, df, **kwargs):
    """Toytree plot"""

    import toyplot
    from tracebtb import trees
    with open(treefile, 'r') as f:
        s = f.read()
    if len(df)<800:
        canvas = trees.draw_tree(s, df, **kwargs)
        toyplot.html.render(canvas, "temp.html")
    else:
        return '<h2>too many tips</b>'
    with open('temp.html', 'r') as f:
        html = f.read()
    return html

def phylocanvas_tree(treefile, df, col):
    """Draw newick tree"""

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
    sub['marker'] = sub.Species.map(speciesmarkers)
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
    #splitview = pn.Row(bokeh_plot.split_view(sub, 'HERD_NO', sp),height=800, styles=CARD_STYLE)

    #main = pn.Column(pn.Row(plot, treepane, overview, height=600),pn.Row(table,timeline))
    main = pn.FlexBox(pn.Row(overview, plot, treepane, height=600), table, timeline)
    report = pn.Column(header, main)
    return report

def dashboard(meta, parcels, moves=None, lpis_cent=None,
              snpdist=None, lpis_master_file=None,
              selections={}):
    """
    Dashboard app with panel for tracebtb.
    Args:
        meta: geodataframe of meta with locations of samples
        parcels: parcels for relevant samples from lpis
        moves: dataframe of moves
        lpis_cent: centroids from LPIS
    """

    lpis = None
    coords = None
    view_history = []
    def update(event=None, sub=None):
        """Update selection"""

        #optimise by checking if sub has changed from before?
        provider = provider_input.value
        cmap = cmap_input.value
        col = colorby_input.value
        ms = markersize_input.value
        legend = legend_btn.value

        global selected
        if sub is None:
            sub = selected

        if len(sub[col].unique())>20:
            cmap=None
        sub['color'],c = tools.get_color_mapping(sub, col, cmap)
        sub['marker'] = sub.Species.map(speciesmarkers)
        #set global selected
        selected = sub
        view_history.append(sub)
        #filters
        sub = apply_filters(sub)
        info_pane.object = f'**{len(sub)} samples**'
        update_overview(sub=sub)

        #get moves first so we know which parcels to show
        herds = list(sub.HERD_NO)
        mov = tools.get_moves_bytag(sub, moves, lpis_cent)
        if mov is not None and moves_btn.value == True:
            herds.extend(mov.move_to)
        sp = parcels[parcels.SPH_HERD_N.isin(herds)].copy()
        pcmap='tab20'
        if len(sp.SPH_HERD_N.unique())>20:
            pcmap=None
        sp['color'],c = tools.get_color_mapping(sp, 'SPH_HERD_N',pcmap)

        if parcels_btn.value == True:
            pcls=sp
            if mov is not None:
                herds.extend(mov.move_to)
        else:
            pcls=None
        labels = parcellabel_btn.value
        p = bokeh_plot.plot_selection(sub, pcls, provider=provider, ms=ms,
                                      col=col, legend=legend, labels=labels)
        if showcounties_btn.value == True:
             bokeh_plot.plot_counties(p)
        if kde_btn.value == True:
            bokeh_plot.kde_plot_groups(sub, p, col, 6)
        global lpis
        if neighbours_btn.value == True and lpis is not None:
            #get neighbours
            #nbr = find_neighbours(sub, 800, lpis_cent, lpis)
            shb = shared_borders(sp, lpis)
            bokeh_plot.plot_lpis(shb, p)

        mov = tools.get_moves_bytag(sub, moves, lpis_cent)
        if moves_btn.value == True:
            bokeh_plot.plot_moves(p, mov, lpis_cent)
        if mov is not None:
            moves_pane.value = mov.reset_index().drop(columns=['geometry'])
        else:
            moves_pane.value = pd.DataFrame()
        plot_pane.object = p

        #change selection table
        selected_pane.value = sub.drop(columns=['geometry'])
        def highlight(x):
            color = speciescolors[x]
            return 'background-color: %s' % color
        selected_pane.style.apply(highlight, subset=['Species'], axis=1)

        #timeline
        herdcolors = dict(zip(sp.SPH_HERD_N,sp.color))
        p = bokeh_plot.plot_moves_timeline(mov, herdcolors)
        timeline_pane.object = p
        timeline_pane.param.trigger('object')

        if tree_btn.value ==True:
            update_tree(sub=sub)
        #summary plots
        update_catplot(sub=sub)
        #herds
        update_herd_summary()
        #cluster summary
        update_cluster_summary()
        #dist matrix
        #dm = snpdist.loc[sub.index,sub.index]
        #p = bokeh_plot.heatmap(dm)
        #snpdist_pane.object = p
        #snpdist_pane.param.trigger('object')

        #locked coords
        '''global coords
        if lockbtn.value == True:
            if coords != None:
                print (coords)
                xmin, xmax, ymin, ymax = coords
                p.x_range.start = xmin
                p.x_range.end = xmax
                #p.y_range.start = ymin
                #p.y_range.end = ymax
            else:
                coords = get_figure_coords(p)
                print (coords)
        else:
            coords = None'''

        return

    def update_herd_summary(event=None):
        """herd summary"""

        global selected
        h = tools.herd_summary(selected, moves, snpdist)
        herds_pane.value = h
        return

    def update_cluster_summary(event=None):

        global selected
        col = groupby_input.value
        cl = tools.cluster_summary(selected, col, 5, snpdist)
        clusters_pane.value = cl
        return

    def update_overview(event=None, sub=None):

        fig = plot_overview(sub)
        overview_pane.object = fig
        plt.close(fig)
        return

    def load_lpis(event=None):
        global lpis
        if lpis_master_file != None:
            lpis = gpd.read_file(lpis_master_file).set_crs('EPSG:29902')
        return

    def apply_filters(df):
        """Filter from widgets"""

        minsize = clustersizeslider.value
        key = groupby_input.value
        groups = df[key].value_counts()
        groups = groups[groups>=minsize]
        df = df[df[key].isin(groups.index)].copy()
        #print (df)
        start, end = timeslider.value
        df = df[(df.Year>=start) & (df.Year<=end)].copy()
        if homebredbox.value == True:
            df = df[df.Homebred=='yes'].copy()
        return df

    def set_provider(event=None):
        """Change map provider"""

        p = plot_pane.object
        #remove old tile
        p.renderers = [x for x in p.renderers if not str(x).startswith('TileRenderer')]
        provider = provider_input.value
        p.add_tile(provider, retina=True)

    def update_groups(event=None):
        """Change group choices"""

        groupby = groupby_input.value
        vals = pd.DataFrame(meta[groupby].value_counts())
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
        from pyproj import Transformer
        transformer_to_latlon = Transformer.from_crs("EPSG:3857", "EPSG:29902", always_xy=True)
        xmin, ymin = transformer_to_latlon.transform(xmin, ymin)
        xmax, ymax = transformer_to_latlon.transform(xmax, ymax)
        sub = meta.cx[xmin:xmax, ymin:ymax]
        update(sub=sub)
        return

    def select_radius(event=None):
        """Select within radius of a center"""

        global selected
        point = selected.union_all().centroid
        distances = selected.geometry.apply(lambda x: x.distance(point))
        radius = distances.median()
        sub = meta[meta.geometry.distance(point) <= radius].copy()
        update(sub=sub)
        return

    def select_from_table(event=None):
        """Select samples from table"""

        df = meta_pane.selected_dataframe
        sub = meta.loc[df.index]
        update(sub=sub)
        return

    def select_related(event=None, df=None):
        """Find related samples to selected indexes i.e. within n snps"""

        df = meta_pane.selected_dataframe
        idx = list(df.index)
        dist = threshold_input.value
        names=[]
        for i in idx:
            found = tools.get_within_distance(snpdist, i, dist)
            names.extend(found)
        names = list(set(names))
        sub = meta.loc[names]
        update(sub=sub)
        return

    def find_related(event=None):
        """Find related samples"""

        global selected
        col = groupby_input.value
        cl = list(selected[col].unique())
        sub = meta[(meta[col].isin(cl))].copy()
        update(sub=sub)
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

    def update_tree(event=None, sub=None):
        """Update with toytree"""

        #replace this with just selecting tips on a main tree using another package..
        if len(sub)<=1 or len(sub)>500:
            html = '<h1><2 or too many samples</h1>'
            tree_pane.object = html
            return
        treefile = get_tree(snpdist, sub.index)
        col = colorby_input.value
        html = draw_toytree(treefile, sub,
                            tiplabelcol=tiplabel_input.value, markercol='Species', height=400)
        #html = phylocanvas_tree(treefile, sub, col)
        tree_pane.object = html
        return

    def do_search(event=None):

        query = search_input.value
        col = searchcol_select.value
        found = meta[meta[col]==query]
        meta_pane.value = found.drop(columns=['geometry'])
        return

    def reset_table(event=None):
        meta_pane.value = meta
        return

    def create_report(event=None):

        global selected
        sub = selected
        sp = parcels[parcels.SPH_HERD_N.isin(sub.HERD_NO)].copy()
        sp['color'],c = tools.get_color_mapping(sp, 'SPH_HERD_N')
        mov = tools.get_moves_bytag(sub, moves, lpis_cent)
        col = colorby_input.value
        result = report(selected, sp, mov, col, lpis_cent, snpdist)
        #report_file = f'report_{datetime.today().date()}.html'
        result.save(report_file)
        return

    def load_selection(event=None):

        global selected
        name = selections_input.value
        idx = selections[name]['indexes']
        selected = meta.loc[idx]
        update(sub=selected)
        return

    def delete_selection(event=None):

        global selected
        name = selections_input.value
        del (selections[name])
        selections_input.options = list(selections.keys())
        with open(selections_file,'w') as f:
            f.write(json.dumps(selections))
        return

    def save_selection(event=None):
        """Save a selection"""

        global selected
        name = selectionname_input.value
        if name == '':
            print ('empty name')
            return
        selections[name] = {}
        selections[name]['indexes'] = list(selected.index)
        #add to menu
        selections_input.options = list(selections.keys())
        #save to file
        with open(selections_file,'w') as f:
            f.write(json.dumps(selections))
        return

    '''def update_scatter(event=None):
        """Update categorical strip plot"""

        global selected
        p = bokeh_plot.cat_plot(selected, row_input.value, col_input.value,
                                ms=mscat_input.value, marker='marker')
        scatter_pane.object = p
        return'''

    def update_catplot(event=None, sub=None):
        """Update categorical plot"""

        if sub is None:
            global selected
            sub = selected
        x = catx_input.value
        y = caty_input.value
        hue = cathue_input.value
        kind = catkind_input.value
        cmap = cmap_input.value
        row = None
        import seaborn as sns
        cg = sns.catplot(data=sub, x=x, y=y, hue=hue, kind=kind, aspect=2,
                         palette=cmap, dodge=True)
        catplot_pane.object = cg.fig
        plt.close(cg.fig)
        #print (fig)
        return

    w=140
    #main panes
    plot_pane = pn.pane.Bokeh()
    #overview_pane = pn.pane.Bokeh(height=300)
    overview_pane = pn.pane.Matplotlib(height=300)
    split_pane = pn.pane.Bokeh(sizing_mode='stretch_both')
    timeline_pane = pn.pane.Bokeh(sizing_mode='stretch_height')
    info_pane = pn.pane.Markdown('', styles={'color': "red"})

    #main table
    meta_pane = pnw.Tabulator(disabled=True,page_size=100,
                              frozen_columns=['sample'],
                              sizing_mode='stretch_both')
    showselected_btn = pnw.Button(name='Select Samples', button_type='primary', align="end")
    pn.bind(select_from_table, showselected_btn, watch=True)
    selectrelated_btn = pnw.Button(name='Find Related', button_type='primary', align="end")
    pn.bind(select_related, selectrelated_btn, watch=True)
    threshold_input = pnw.IntInput(name='Threshold', value=7, step=1, start=2, end=20,width=60)
    #search
    scols = ['sample','Year','HERD_NO','Animal_ID','Species','County','IE_clade','Region']
    search_input = pnw.TextInput(name="Search", value='',width=150)
    searchcol_select = pnw.Select(name='Column',value='HERD_NO',options=scols,width=100)
    search_btn = pnw.Button(icon=get_icon('search'), icon_size='1.8em', align="end")
    pn.bind(do_search, search_btn, watch=True)
    reset_btn = pnw.Button(icon=get_icon('refresh'), icon_size='1.8em', align="end")
    pn.bind(reset_table, reset_btn, watch=True)

    table_widgets = pn.Row(showselected_btn, selectrelated_btn, threshold_input,
                           search_input, searchcol_select, search_btn, reset_btn,
                           sizing_mode='stretch_width')
    table_pane = pn.Column(meta_pane,table_widgets,sizing_mode='stretch_both')

    #selected table
    selected_pane = pnw.Tabulator(show_index=False,disabled=True,page_size=100,
                                  frozen_columns=['sample'],sizing_mode='stretch_both')

    moves_pane = pnw.Tabulator(show_index=False,disabled=True,page_size=100,
                                  frozen_columns=['tag'],sizing_mode='stretch_both')
    herds_pane = pnw.Tabulator(show_index=False,disabled=True,page_size=100,
                               frozen_columns=['HERD_NO'],sizing_mode='stretch_both')
    clusters_pane = pnw.Tabulator(show_index=False,disabled=True,page_size=100,
                                  frozen_columns=['cluster'],sizing_mode='stretch_both')

    tree_pane = pn.pane.HTML(height=400)
    network_pane = pn.pane.Bokeh()
    snpdist_pane = pn.pane.Bokeh(height=400)

    cols = [None]+tools.get_ordinal_columns(meta)
    #selections pane
    selections_input = pnw.Select(name='Selections',options=list(selections.keys()),value='',width=w)
    loadselection_btn = pnw.Button(name='Load Selection', button_type='primary',width=w)
    pn.bind(load_selection, loadselection_btn, watch=True)
    deleteselection_btn = pnw.Button(name='Delete Selection', button_type='primary',width=w)
    pn.bind(delete_selection, deleteselection_btn, watch=True)
    saveselection_btn = pnw.Button(name='Save Selection', button_type='primary',width=w)
    pn.bind(save_selection, saveselection_btn, watch=True)
    selectionname_input = pnw.TextInput(value='',width=w)

    card2 = pn.Column('## Selections', pn.Row(pn.Column(selections_input, selectionname_input),
                      pn.Column(loadselection_btn, deleteselection_btn, saveselection_btn)),
                       width=340, styles=card_style)
    #reporting
    doreport_btn = pnw.Button(name='Generate', button_type='primary')
    savereport_btn = pnw.FileDownload(file=report_file, button_type='success', auto=False,
                                    embed=False, name="Download Report")
    card3 = pn.Column('## Reporting', doreport_btn, savereport_btn, width=340, styles=card_style)
    pn.bind(create_report, doreport_btn, watch=True)

    #lpis
    loadlpis_btn = pnw.Button(name='Load LPIS', button_type='primary')
    pn.bind(load_lpis, loadlpis_btn, watch=True)
    card4 = pn.Column('## LPIS', loadlpis_btn, width=340, styles=card_style)

    utils_pane = pn.FlexBox(*[card2, card3, card4], flex_direction='column', min_height=400,
                             styles={'margin': '10px'}, sizing_mode='stretch_both')

    #widgets
    groupby_input = pnw.Select(name='group by',options=cols,value='snp7',width=w)
    groups_table = pnw.Tabulator(disabled=True, widths={'index': 70}, layout='fit_columns',pagination=None, height=250, width=w)
    colorby_input = pnw.Select(name='color by',options=cols,value='snp7',width=w)
    cmap_input = pnw.Select(name='colormap',options=colormaps,value='Set1',width=w)
    provider_input = pnw.Select(name='provider',options=['']+bokeh_plot.providers,value='CartoDB Positron',width=w)
    markersize_input = pnw.FloatInput(name='marker size', value=10, step=1, start=2, end=100,width=w)
    tiplabel_input = pnw.Select(name='tip label',options=list(meta.columns),value='sample',width=w)
    widgets = pn.Column(pn.WidgetBox(groupby_input,groups_table,colorby_input,cmap_input,tiplabel_input,
                                     provider_input,markersize_input),info_pane,width=w+30)
    #toolbar
    #split_btn = pnw.Button(icon=get_icon('plot-grid'), description='split view', icon_size='1.8em')
    selectregion_btn = pnw.Button(icon=get_icon('plot-region'), description='select in region', icon_size='1.8em')
    selectradius_btn = pnw.Button(icon=get_icon('plot-centroid'), description='select within radius', icon_size='1.8em')
    tree_btn = pnw.Toggle(icon=get_icon('tree'), icon_size='1.8em')
    parcels_btn = pnw.Toggle(icon=get_icon('parcels'), icon_size='1.8em')
    moves_btn = pnw.Toggle(icon=get_icon('plot-moves'), icon_size='1.8em')
    legend_btn = pnw.Toggle(icon=get_icon('legend'), icon_size='1.8em')
    neighbours_btn = pnw.Toggle(icon=get_icon('neighbours'), icon_size='1.8em')
    parcellabel_btn = pnw.Toggle(icon=get_icon('parcel-label'), icon_size='1.8em')
    showcounties_btn = pnw.Toggle(icon=get_icon('counties'), icon_size='1.8em')
    kde_btn = pnw.Toggle(icon=get_icon('contour'), icon_size='1.8em')
    #lockbtn = pnw.Toggle(icon=get_icon('lock'), icon_size='1.8em')
    toolbar = pn.Column(pn.WidgetBox(selectregion_btn,selectradius_btn,tree_btn,
                                     parcels_btn,parcellabel_btn,showcounties_btn,moves_btn,legend_btn,
                                     neighbours_btn,kde_btn),width=70)

    #option below plot
    timeslider = pnw.IntRangeSlider(name='Time',width=150,
                    start=2000, end=2024, value=(2000, 2024), step=1)
    clustersizeslider = pnw.IntSlider(name='Min. Cluster Size',width=150,
                    start=1, end=20, value=1, step=1)
    homebredbox = pnw.Checkbox(name='Homebred',value=False)
    findrelated_btn = pnw.Button(name='Find Related', button_type='primary', align="end")
    pn.bind(find_related, findrelated_btn, watch=True)

    filters = pn.Row(findrelated_btn,timeslider,clustersizeslider,homebredbox)

    groupby_input.param.watch(update_groups, 'value')
    groups_table.param.watch(select_group, 'selection')
    provider_input.param.watch(set_provider, 'value')
    colorby_input.param.watch(update, 'value')
    cmap_input.param.watch(update, 'value')
    tiplabel_input.param.watch(update, 'value')
    tree_btn.param.watch(update, 'value')
    parcels_btn.param.watch(update, 'value')
    parcellabel_btn.param.watch(update, 'value')
    showcounties_btn.param.watch(update, 'value')
    moves_btn.param.watch(update, 'value')
    legend_btn.param.watch(update, 'value')
    neighbours_btn.param.watch(update, 'value')
    kde_btn.param.watch(update, 'value')
    timeslider.param.watch(update, 'value')
    clustersizeslider.param.watch(update, 'value')
    homebredbox.param.watch(update, 'value')

    #pn.bind(update_tree, treebtn, watch=True)
    #pn.bind(split_view, split_btn, watch=True)
    pn.bind(select_region, selectregion_btn, watch=True)
    pn.bind(select_radius, selectradius_btn, watch=True)

    #categorical plot widget
    catplot_pane = pn.pane.Matplotlib(height=350,sizing_mode='stretch_both')
    catx_input = pnw.Select(name='x',options=cols,value='Species',width=w)
    caty_input = pnw.Select(name='y',options=cols,value=None,width=w)
    cathue_input = pnw.Select(name='hue',options=cols,value=None,width=w)
    kinds = ['count','bar','strip','swarm']
    catkind_input = pnw.Select(name='kind',options=kinds,
                               value='count',width=w)
    cathue_input = pnw.Select(name='hue',options=cols,value='County',width=w)
    catx_input.param.watch(update_catplot, 'value')
    caty_input.param.watch(update_catplot, 'value')
    cathue_input.param.watch(update_catplot, 'value')
    catkind_input.param.watch(update_catplot, 'value')
    analysis_pane1 = pn.Column(pn.Row(catx_input,caty_input,cathue_input,catkind_input),
                               catplot_pane,height=600)

    app = pn.Column(
                pn.Row(widgets,
                       toolbar,
                pn.Column(
                    pn.Tabs(('Map',pn.Column(plot_pane,filters)),
                            #('Split View',pn.Column(split_pane, height=300)),
                            ('Summary I',analysis_pane1),
                            ('Main Table', table_pane),
                            ('Selected', selected_pane),
                            ('Moves', moves_pane),
                            ('Herds',herds_pane),('Groups',clusters_pane),
                            ('Tools',utils_pane),
                             dynamic=True,
                             sizing_mode='stretch_both'),
                          ),
                pn.Column(pn.Tabs(('Overview',overview_pane),
                                  ('Tree',tree_pane),
                                  #('Network',network_pane),
                                  dynamic=True,width=500),
                                  timeline_pane, width=500),
             max_width=2000,min_height=600
    ))
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

    #load config file
    if not os.path.exists(configfile):
        d = {'dashboard':{'lpis_master_file':''}}
        with open(configfile, "w") as outfile:
            json.dump(d, outfile)
        lpis_master_file = None
    else:
        with open(configfile) as f:
            jsondata = json.load(f)
        print('found settings file')
        lpis_master_file = jsondata['dashboard']['lpis_master_file']

    data = pickle.load(open(args.project,'rb'))
    meta = data['meta']#.to_crs('EPSG:3857')
    moves = data['moves']
    lpis_cent = data['lpis_cent']
    parcels = data['parcels']
    snpdist = data['snpdist']
    #selections = data['selections']
    print (selections_file)
    if os.path.exists(selections_file):
        selections = json.load(open(selections_file,'r'))
    else:
        selections = {}
    #create template
    bootstrap = pn.template.BootstrapTemplate(title='TracebTB',
                        favicon=logoimg,logo=logoimg,header_color='green')
    #pn.config.sizing_mode = 'stretch_width'
    app = dashboard(meta, parcels, moves, lpis_cent, snpdist, lpis_master_file, selections)
    app.project_file = args.project
    bootstrap.main.append(app)
    bootstrap.servable()
    pn.serve(bootstrap, port=5010)
             #basic_auth={'guest':"mbovis"}, cookie_secret='cookie_secret')
             #websocket_origin=['bola.ucd.ie:80','localhost:5010'])

if __name__ == '__main__':
    main()
