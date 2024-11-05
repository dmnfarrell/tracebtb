#!/usr/bin/env python3

"""
    TraceBTB panel dashboard launcher module
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

from tracebtb import dashboards, bokeh_plot

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
defaults = {'dashboard':{'lpis_master_file':'','tree_file':None}}

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

def check_settings_file(filename):
    """Check settings file"""

    print (filename)
    with open(filename) as f:
        settings = json.load(f)['dashboard']
    for key in defaults['dashboard']:
        if not key in settings:
            settings[key] = defaults['dashboard'][key]
    print (settings)
    return settings

def main():
    "Run the application"
    from argparse import ArgumentParser
    parser = ArgumentParser(description='TracebTB')
    parser.add_argument("-p", "--proj", dest="project",default=None,
                            help="load project file", metavar="FILE")
    parser.add_argument("-s", "--settings", dest="settings",default=configfile,
                            help="load a json settings file", metavar="FILE")
    parser.add_argument("-i", "--interface", dest="interface",default='full',
                            help="type of dashboard to launch", metavar="FILE")
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
            settings = check_settings_file(args.settings)
            lpis_master_file = settings['lpis_master_file']
            treefile = settings['tree_file']

        data = pickle.load(open(args.project,'rb'))
        meta = data['meta']
        moves = data['moves']
        lpis_cent = data['lpis_cent']
        parcels = data['parcels']
        snpdist = data['snpdist']
        testing = data['testing']

        if os.path.exists(selections_file):
            selections = json.load(open(selections_file,'r'))
        else:
            selections = {}
        if os.path.exists(layers_file):
            layers = layers_from_file(layers_file)
        else:
            layers = {}

        def create_app():
            #Generate a new dashboard instance per session
            if args.interface == 'full':
                title='TracebTB'
                bkgr='#4B7CC1 '
                app = dashboards.FullDashboard(meta, parcels, moves, lpis_cent,
                                               snpdist, lpis_master_file,
                                                treefile, testing, selections, layers)
            elif args.interface == 'simple':
                title = 'TracebTB Query'
                bkgr='#30833C'
                app = dashboards.SimpleQueryDashboard(meta, parcels, moves, lpis_cent,
                                               snpdist, lpis_master_file,
                                                treefile, testing, selections, layers)
            app.project_file = args.project
            app.settings = settings
            app.treefile = treefile
            layout = app.show()
            # Create a session-specific app function
            bootstrap = pn.template.BootstrapTemplate(
                title=title,
                favicon=logoimg,
                logo=logoimg,
                header_color='white',
                header_background=bkgr
            )
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
            #default - no security
            pn.serve(create_app, port=5010, websocket_origin=["localhost:5010"],
                     cookie_secret='cookie_secret')

if __name__ == '__main__':
    main()


