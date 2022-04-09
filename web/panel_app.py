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
import json
import numpy as np
import pandas as pd
import string
import sqlite3
from collections import OrderedDict
import toytree, toyplot
from snipgenie import snp_typing, trees, tools, plotting
import geopandas as gpd

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
#DATABASE = '../notebooks/'

parcels = gpd.read_file('/storage/btbgenie/monaghan/metadata/lpis_monaghan_10km_buff.shp').set_crs('EPSG:29902')
parcels = parcels.to_crs(epsg=3395)

meta = pd.read_csv('ireland_test_data.csv')
meta = meta.fillna('')
snpcols=['snp3','snp5','snp12','snp20','snp50']
for i in snpcols:
    meta[i]=meta[i].astype(str)

species_colors = {'Cow':'green', 'Badger':'blue', 'Deer':'red','Dog':'orange'}
c,sb_colors = plotting.get_color_mapping(meta, 'SB', cmap='Set1')
c,county_colors = plotting.get_color_mapping(meta, 'county', cmap=None)
c,snp_colors = plotting.get_color_mapping(meta, 'snp12', cmap=None)

cmaps = {'species': species_colors,'SB':sb_colors,'county':county_colors,'snp12':snp_colors}
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

def draw_tree(tre, df, colorby=None, layout='r', font_size=10, node_size=8, root=None, labels=True, **kwargs):
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
        node_sizes=[0 if i else node_size for i in tre.get_node_values(None, 1, 0)]
        node_colors = [colormap[mapping[n]] if n in mapping else 'gray' for n in tre.get_node_values('name', True, True)]
        if len(tipnames)>40 or labels == False:
            tipnames = ['' for i in tipnames]
    else:
        tip_colors = None
    tip_labels_style={
            "font-size": "%spx" %font_size,
            "-toyplot-anchor-shift": "13px",
        }
    #render to html
    canvas, axes, mark = tre.draw(tip_labels=tipnames, tip_labels_colors=tip_colors,tip_labels_style=tip_labels_style,
                                  layout=layout,node_colors=node_colors,node_sizes=node_sizes,node_hover=True,
                                  scalebar=True, width=400, height=500, **kwargs);
    toyplot.html.render(canvas, "tree.html");
    return canvas

def dist_matrix(df):
    """Distance matrix from nuc SNP matrix"""

    snpmat = df.snps.apply(snp_typing.decode_snps)
    names = snpmat.index = df.name
    M=[]
    for i,r in snpmat.iterrows():
        x=[]
        s1 = ''.join(r)
        for i,r in snpmat.iterrows():
            s2=''.join(r)
            m = [1 for x,y in zip(s1,s2) if x!=y]
            #print (m)
            x.append(sum(m))
        M.append(x)
    M = pd.DataFrame(M,index=names,columns=names)
    return M

def dist_plot(dist):
    """distance matrix plot"""

    import seaborn as sns
    annot=False
    if len(dist)<12:
        annot=True
    g=sns.clustermap(dist,annot=annot,fmt='.0f',xticklabels=True,yticklabels=True,cmap='Blues',figsize=(10,10))
    return g.fig

def get_cluster_samples(cl, col='snp12'):
    df=meta[meta[col].isin(cl)]
    return list(df.name)

def get_geodatasource(gdf):
    """Get getjsondatasource from geopandas object"""

    json_data = json.dumps(json.loads(gdf.to_json()))
    return GeoJSONDataSource(geojson = json_data)

def bokeh_geodataframe(gdf, p):

    geosource = get_geodatasource(gdf)
    p.patches('xs','ys', source=geosource, fill_alpha=.4, line_width=0.5, fill_color='red', line_color='black')
    return

def wgs84_to_web_mercator(df, lon="LON", lat="LAT"):
     """convert mat long to web mercartor"""

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

def get_moves(df):
    """get moves for a sample"""

    cols=['ANIMAL_ID','HERD_NO','move_from','move_date','time_from_last_bd']
    t = df.merge(allmov,left_on='ANIMAL_ID',right_on='tag',how='inner')[cols]

    t = t.merge(lpis_cent,left_on='move_from',right_on='SPH_HERD_N')
    t = t.sort_values('move_date')
    if len(t)==0:
        return
    x = lpis_cent[lpis_cent.SPH_HERD_N.isin(df.HERD_NO)]
    t = pd.concat([t,x])
    t = t.drop(columns='geometry')
    t = wgs84_to_web_mercator(t, lon="LONG", lat="LAT")
    return t

def bokeh_map(df=None, long=None, lat=None, height=600,
              tile_provider='CARTODBPOSITRON', colorby='County',
              labels=None, arrows=None):
    """Bokeh map"""

    tile_provider = get_provider(tile_provider)
    tools = "pan,wheel_zoom,hover,tap,lasso_select,reset,save"
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
    df['size'] = 10
    source = ColumnDataSource(df)

    if len(df)==1:
        x_range=(x-200000, x+200000)
        y_range=(y-200000, y+200000)
    else:
        x_range=None; y_range=None

    #draw figure
    p = figure(tools=tools, #x_range=x_range,y_range=y_range,
               #x_axis_type="mercator", y_axis_type="mercator",
               plot_width=height, plot_height=height, sizing_mode=sizing_mode)#, active_scroll='wheel_zoom')
    p.add_tile(tile_provider)
    p.circle(x='x', y='y', size='size', alpha=0.7, color='color', source=source,
                nonselection_fill_alpha=0.7, selection_fill_color="red", selection_line_color="black")
                #, legend_group=colorby)

    labels = LabelSet(x='x', y='y', text='label',text_font_size='10pt',
                     x_offset=5, y_offset=5, source=source, render_mode='canvas')
    p.add_layout(labels)

    #arrows for moves
    a_source = ColumnDataSource(data=dict(
            x=[], y=[], x2=[], y2=[] ))

    arrows = Arrow(end=VeeHead(fill_color="black", size=10, line_alpha=.9),
                       x_start='x', y_start='y',
                       x_end='x2', y_end='y2',
                       source=source)
    p.renderers.append(arrows)
    p.add_layout(arrows)

    #p.legend.location = "top_left"
    #p.legend.click_policy="mute"
    p.toolbar.logo = None
    p.match_aspect = True
    hover = p.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([
        ("name", "@sample"),
        ("species", "@Species"),
        ("SB", "@SB"),
        ("snp12", "@snp12"),
        ("snp50", "@snp50"),
        ("county", "@County"),
    ])
    return p


def map_dash(df):
    """Map dashboard"""

    names = sorted(list(df['sample'].unique()))
    snp12 = sorted(list(df.snp12.unique()))
    tre = None
    sel = None
    cols = df.columns[:6]
    cats=['County','Species','snp12','snp50','SB']
    labels=['','name','snp12','snp50','SB']
    counties=['All','Wicklow','Monaghan','NI']
    map_pane = pn.pane.Bokeh(width=400)
    tree_pane = pn.pane.HTML(width=300)
    plot_pane = pn.pane.Matplotlib(height=500)

    tile_select = pnw.Select(name='tile layer',options=providers,width=200)
    colorby_select = pnw.Select(name='color by',options=cats,width=200)
    label_select = pnw.Select(name='label',options=labels,width=200)
    name_select = pnw.MultiSelect(name='name',options=names,size=6,width=200)
    cluster_select = pnw.MultiSelect(name='cluster',options=snp12,size=6,width=200)
    county_select = pnw.Select(name='county',options=counties,width=200)
    show_parcels_box = pnw.Checkbox(name="show land parcels",value=False,width=200)

    #buttons
    reset_btn = pnw.Button(name='Reset', button_type='primary',width=200)
    related_btn = pnw.Button(name='Find Related Isolates', button_type='primary',width=200)
    related_clusters_btn = pnw.Button(name='Related Clusters', button_type='primary',width=200)
    moves_btn = pnw.Button(name='Show Moves', button_type='primary',width=200)
    outliers_btn = pnw.Button(name='Find outliers', button_type='primary',width=200)
    help_btn = pnw.Button(name='Help', button_type='primary',width=200)

    style1 = {'background':'lightgray','padding':'5px','font-family':'monospace'}
    info_pane = pn.pane.HTML(style=style1, width=200, height=200,sizing_mode='stretch_both')
    df_pane = pn.pane.DataFrame(df[cols],width=500,height=200,sizing_mode='scale_both',max_rows=20,index=False)
    snps_pane = pnw.Tabulator()
    meta_pane = pnw.Tabulator(df[cols],pagination='remote', page_size=20, width=500,height=600)

    #settings pane
    tree_layout_select = pnw.Select(name='tree layout',options=['r','c','d'],width=200)
    root_select = pnw.Select(name='root on',options=[''],width=200)
    tip_align_box = pnw.Checkbox(name="tip labels align",value=False)
    edge_type_select = pnw.Select(name= "edge type",options=['p','b','c'],width=100)
    node_size_entry = pnw.IntInput(name="tree node size", step=1, start=1, end=30, value=8)
    point_size_entry = pnw.IntInput(name="map point size", step=1, start=2, end=30, value=10)
    settings_pane = pn.Column(tree_layout_select,root_select,tip_align_box,edge_type_select,node_size_entry,point_size_entry,width=100)
    left_tab_pane = pn.Tabs(('tree', tree_pane), ('dists',plot_pane), ('dist matrix',snps_pane), ('meta', meta_pane), ('settings',settings_pane), ('debug',info_pane))

    empty_pane = pn.pane.HTML(width=300,style=style1,sizing_mode='scale_height')
    empty_pane.object = 'test'

    loading = pn.indicators.LoadingSpinner(value=False, width=100, height=100)

    def update_info(attr,new,old):
        #print(new,old)
        info_pane.object = '<p>%s,%s</p>' %(int(new),int(old))

    def zoom_to_points(sel, p, pad=200000):
        #zoom with aspect conserved
        if len(sel)==1:
            x=sel.iloc[0].x
            y=sel.iloc[0].y
        else:
            x=-779236
            y=7076025
        p.x_range.update(start=x-pad,end=x+pad)
        p.y_range.update(start=y-pad,end=y+pad)
        return

    def names_selected(event):
        names = name_select.value
        items_selected(names)

    def clusters_selected(event):
        clusters = cluster_select.value
        names = get_cluster_samples(df, clusters)
        items_selected(names)

    def items_selected(names):

        global tre, sel
        root_select.options = ['']+names
        p = map_pane.object
        source = p.renderers[1].data_source
        colorby = colorby_select.value
        colormap = cmaps[colorby]
        df['color'] = [colormap[i] if i in colormap else 'gray' for i in df[colorby]]
        df['size'] = point_size_entry.value
        #selected data
        sel = df[df.name.isin(names)]
        info_pane.object = sel[cols]
        #show these points only on map
        source.data = dict(sel)

        #draw parcels?
        if show_parcels_box.value==1:
            lp = parcels.sample(200) #[lpis_p.SPH_HERD_N.isin(sel.HERD_NO)]
            bokeh_geodataframe(lp, p)

        #zoom to points selected
        zoom_to_points(sel, p)
        update_tile()

        #get a tree
        if len(sel)>=3:
            loading.value = True
            tre = get_tree(sel)
            canvas = draw_tree(tre, sel, colorby, layout=tree_layout_select.value,
                               tip_labels_align=tip_align_box.value, edge_type=edge_type_select.value,
                               root=root_select.value)
            tree_pane.object = canvas
            #dist matrix
            m = dist_matrix(sel)
            plot_pane.object = dist_plot(m)
            snps_pane.object = m
            loading.value = False
        else:
            tree_pane.object = ''
            snps_pane.object = ''
        return

    def county_selected(event):

        global sel
        county = county_select.value
        if county == 'All':
            sel = df
        else:
            sel = df[df.County==county]
        update_map(event)
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
            canvas = draw_tree(tre, sel, colorby_select.value, layout=tree_layout_select.value,
                               tip_labels_align=tip_align_box.value, edge_type=edge_type_select.value,
                               root=root_select.value)
            tree_pane.object = canvas
            #dist matrix
            m = dist_matrix(sel)
            plot_pane.object = dist_plot(m)
            snps_pane.object = m
            loading.value = False
            #info_pane.object = sel
        return

    def tap_callback(event):
        """tap tool callback"""

        p = map_pane.object
        source = p.renderers[1].data_source
        ind = source.selected.indices
        info_pane.object = str(ind)
        df = pd.DataFrame(source.data)
        sel = df.iloc[ind]
        return

    def find_related_callback(event):

        global sel
        p = map_pane.object
        source = p.renderers[1].data_source
        ind = source.selected.indices
        info_pane.object = str(ind)
        df = pd.DataFrame(source.data)
        sel = df.iloc[ind]
        if len(sel)>0:
            names = sel['sample']
            sel = find_related(names)
            if sel is None:
                return
        update_map(event)
        return

    def find_related(names):
        """Find related isolates"""

        cl = df[df.name.isin(names)].snp12.unique()
        if len(cl)==0:
            return
        sub = df[df.snp12.isin(cl)]
        return sub

    def draw_map(event):
        """Redraw from scratch"""

        global sel
        sel = None
        cluster_select.value = []
        p = map_pane.object = bokeh_map(df)
        #p.x_range.on_change('start', update_info)
        source = p.renderers[1].data_source
        source.selected.on_change('indices', points_selected)
        p.on_event('tap',tap_callback)
        tree_pane.object = ''
        plot_pane.object = None

        return

    def update_tile(event=None):

        p = map_pane.object
        p.renderers = [x for x in p.renderers if not str(x).startswith('TileRenderer')]
        rend = renderers.TileRenderer(tile_source= get_provider(tile_select.value))
        p.renderers.insert(0, rend)

    def update_map(event):
        """Updated colors or labels without replotting"""

        global sel,df
        p = map_pane.object
        source = p.renderers[1].data_source
        colorby = colorby_select.value
        colormap = cmaps[colorby]
        if sel is not None:
            d = sel
        else:
            d = df
        info_pane.object = d[['sample','snp12']]
        d['color'] = [colormap[i] if i in colormap else 'gray' for i in d[colorby]]
        d['size'] = point_size_entry.value
        if label_select.value != '':
            d['label'] = d[label_select.value]
        else:
            d['label'] = ''

        source.data = dict(d)
        update_tree(event)
        return

    def update_tree(event):

        global tre
        #use subset of samples if selected
        names = name_select.value
        if tre != None:
            sel = df[df.name.isin(names)]
            canvas = draw_tree(tre, sel, colorby=colorby_select.value, layout=tree_layout_select.value,
                               tip_labels_align=tip_align_box.value, edge_type=edge_type_select.value,
                               node_size=node_size_entry.value,
                               root=root_select.value)
            tree_pane.object = canvas
        return

    def show_moves(event):

        global sel
        p = map_pane.object
        #arrows source
        if sel is not None:
            d = sel
        else:
            d = df
        a_source = p.renderers[2].source
        t = get_moves(sel)
        if t is None:
            return
        info_pane.object = t
        coords = get_coords_data(t)
        info_pane.object = str(coords)
        a_source.data = coords
        return

    draw_map(None)
    reset_btn.on_click(draw_map)
    moves_btn.on_click(show_moves)
    related_btn.on_click(find_related_callback)
    help_btn.js_on_click(args={'target':'https://github.com/dmnfarrell/btbgenietools' }, code='window.open(target)')

    #label_box = pnw.Checkbox(name='Show labels')
    tile_select.param.watch(update_tile,'value')
    colorby_select.param.watch(update_map,'value')
    label_select.param.watch(update_map,'value')
    name_select.param.watch(names_selected,'value')
    cluster_select.param.watch(clusters_selected,'value')
    county_select.param.watch(county_selected,'value')
    #show_parcels_box.param.watch(update_map,'value')
    tree_layout_select.param.watch(update_tree,'value')
    root_select.param.watch(update_tree,'value')

    #layout dashboard
    app = pn.Column(pn.Row(reset_btn,related_btn,moves_btn,outliers_btn,help_btn),
                    pn.Row(pn.Column(tile_select,colorby_select,label_select,name_select,cluster_select,county_select,show_parcels_box,
                                     background='whitesmoke',sizing_mode='stretch_height'),
                           pn.Column(map_pane,width=600,sizing_mode='stretch_both'),pn.Column(left_tab_pane,width=350,sizing_mode='stretch_height'),loading),
                                     sizing_mode='stretch_both')
    return app

bootstrap = pn.template.BootstrapTemplate(title='BTBGenIE WGS Mapper',
            favicon='static/logo.png',logo='static/logo.png',header_color='blue')
pn.config.sizing_mode = 'stretch_width'
app = map_dash(meta)
bootstrap.main.append(app)
bootstrap.servable()

if __name__ == '__main__':
	pn.serve(bootstrap, port=5000)
