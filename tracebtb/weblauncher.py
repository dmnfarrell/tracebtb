#!/usr/bin/env python3

"""
    TraceBTB panel dashboard launcher module
    Created Aug 2024
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

module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
data_path = os.path.join(module_path,'data')
logo_path = os.path.join(module_path,'logos')
logoimg = os.path.join(logo_path, 'logo.png')

home = os.path.expanduser("~")
configpath = os.path.join(home, '.config','tracebtb')
report_file = 'report.html'
selections_file = os.path.join(configpath,'selections.json')
layers_file = os.path.join(configpath,'layers.gpkg')

from .core import config
from tracebtb import dashboards, bokeh_plot

def layers_from_file(layers_file):
    """load layers from file"""

    import fiona
    layers = {}
    names = fiona.listlayers(layers_file)
    for name in names:
        gdf = gpd.read_file(layers_file, layer=name)
        layers[name] = gdf
    return layers

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

def create_bootstrap_layout_cls(dashboard_cls, title, bkgr, logo=None, **kwargs):
    """
    Returns a callable that creates a fresh BootstrapTemplate
    with the given dashboard class and parameters.

    dashboard_cls: class of your dashboard (e.g., FullDashboard)
    title: page title
    bkgr: header background color
    kwargs: passed to the dashboard constructor
    """
    if logo == None:
        logo = logoimg

    logo_css = """
    .bk-header-logo {
        height: 60px !important;
        width: auto !important;
    }
    """

    def _layout():
        # create a fresh dashboard instance for this session
        app = dashboard_cls(**kwargs)
        layout = app.show()
        bootstrap = pn.template.BootstrapTemplate(
            site='TracebTB',
            title=title,
            favicon=logo,
            logo=logo,
            header_color='white',
            header_background=bkgr
        )
        menu_html = pn.pane.HTML("""
        <div style="display:flex; align-items:center; gap:1.2rem; font-size:.9rem;">
        <a href="/" style="color:white; text-decoration:none;">Main Dashboard</a>
        <a href="/herd" style="color:white; text-decoration:none;">Herd Query</a>
        <a href="/moves" style="color:white; text-decoration:none;">Movements</a>
        </div>
        """ , width=150)
        bootstrap.header.append(menu_html)
        bootstrap.main.append(layout)
        return bootstrap
    return _layout

def main():
    "Run the application"
    from argparse import ArgumentParser
    parser = ArgumentParser(description='TracebTB')
    parser.add_argument("-f", "--proj", dest="project",default=None,
                            help="load project file", metavar="FILE")
    parser.add_argument("-s", "--settings", dest="settings",default='',
                            help="load a json settings file", metavar="FILE")
    parser.add_argument("-p", "--port", dest="port",default=5010,
                            help="port to run server on")
    parser.add_argument("-nl", "--nolpis", dest="nolpis", action="store_true", default=False,
                            help="don't load LPIS (testing mode)")
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
        data = pickle.load(open(args.project,'rb'))
        meta = data['meta']
        moves = data['moves']
        lpis_cent = data['lpis_cent']
        parcels = data['parcels']
        snpdist = data['snpdist']
        testing = data['testing']
        #load config file
        lpis = None
        if os.path.exists(args.settings):
           #we can overwrite settings here
           pass

        settings = config.settings
        lpis_master_file = settings['lpis_master_file']
        treefile = settings['tree_file']

        if args.nolpis == True:
            #for testing only so server launches faster
            lpis_master_file = None
        else:
            lpis = gpd.read_file(lpis_master_file).set_crs('EPSG:29902')
            print ('loaded LPIS')
        if os.path.exists(selections_file):
            selections = json.load(open(selections_file,'r'))
        else:
            selections = {}
        if os.path.exists(layers_file):
            layers = layers_from_file(layers_file)
        else:
            layers = {}

        #Generate multiple dashboard instances and put them on different URLs
        bkgr='#4B7CC1'
        nav_links = pn.Row(
            pn.pane.Markdown("[🏠 Main](/)", sizing_mode="fixed", width=100),
            pn.pane.Markdown("[🐄 Herd](/herd)", sizing_mode="fixed", width=100),
            pn.pane.Markdown("[🚚 Moves](/moves)", sizing_mode="fixed", width=100),
            align="center",
            margin=(10, 20),
        )
        layout1 = create_bootstrap_layout_cls(
            dashboards.FullDashboard,
            'Main',
            '#4B7CC1',
            layers=layers,
            treefile=treefile,
            lpis=lpis,
            selections=selections,
            settings=settings,
            **data
        )
        layout2 = create_bootstrap_layout_cls(
            dashboards.HerdQueryDashboard,
            'Herds',
            "#438328",
            logo=os.path.join(logo_path, 'cow.png'),
            layers=layers,
            treefile=treefile,
            lpis=lpis,
            settings=settings,
            **data
        )
        layout3 = create_bootstrap_layout_cls(
            dashboards.MovesDashboard,
            'Moves',
            "#903E3E",
            logo=os.path.join(logo_path, 'moves.png'),
            layers=layers,
            treefile=treefile,
            lpis=lpis,
            settings=settings,
            **data
        )
        layouts = {'/':layout1, '/herd':layout2, '/moves':layout3}
        port = int(args.port)
        if 'reverse_proxy' in settings:
            #using nginx proxy with basic auth
            s = settings['reverse_proxy']
            pn.serve(layouts, port=port, prefix=s['prefix'],
                websocket_origin=s['origin'], #basic_auth='credentials.json',
                cookie_secret='cookie_secret', show=False)
        else:
            #default - no security
            if port in settings:
                port = settings['port']
            origin = f"localhost:{port}"
            pn.serve(layouts, port=port, websocket_origin=[origin],
                     cookie_secret='cookie_secret', show=False)

if __name__ == '__main__':
    main()
