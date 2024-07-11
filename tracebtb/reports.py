#!/usr/bin/env python

"""
    tracebtb report module.
    Created June 2023
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

import sys,os,shutil,subprocess,time
import glob,random
from datetime import datetime
import tempfile
import pandas as pd
import numpy as np
import pylab as plt
import matplotlib as mpl
from Bio.Align import MultipleSeqAlignment
from . import widgets, plotting, trees, gui

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
data_path = os.path.join(module_path,'data')
logoimg = os.path.join(module_path, 'logo.svg')

css_file = os.path.join(data_path, 'style.css')

def cluster_report(sub, parcels, lpis_cent, moves=None, treefile=None,
                    cmap='tab20', colorcol='snp7', labelcol='Animal_ID',
                    basemap=False, title='TEST', outfile='report.pdf'):
    """Generate html and pdf report for cluster"""

    from weasyprint import HTML

    s = '<html>'
    s += '<head><link rel="stylesheet" href="%s"><head>' %css_file
    #s += '<div class="title"><h3>cluster report (%s)</h3></div>'%datetime.now().strftime("%Y-%m-%d %H:%M")

    cols = ['Animal_ID','HERD_NO','County','Year','Species','']
    #plot map
    fig,ax=plt.subplots(1,1,figsize=(8,8))
    parcels.plot(column='SPH_HERD_N',alpha=0.6,lw=1,cmap=cmap,ax=ax)

    #idx = sub[sub.duplicated('geometry')].index
    gui.plot_single_cluster(sub,col=colorcol,cmap=cmap,ax=ax)
    gui.show_labels(sub, labelcol, ax)
    if basemap == True:
        import contextily as cx
        cx.add_basemap(ax, crs=cent.crs,
                attribution=False, source=gui.providers['OSM'])

    plt.tight_layout()
    tempname = tempfile.mktemp()
    fig.savefig(tempname, format='svg', dpi=150)

    #tree
    treehtml = ''
    import toyplot
    #canvas = trees.draw_tree(treefile, sub, col=colorcol, cmap=cmap, width=200, height=300)
    #toyplot.html.render(canvas, "temp.html")
    #with open('temp.html', 'r') as f:
    #    treehtml = f.read()

    #plot timeline of moves
    if moves is not None:
        fig,ax=plt.subplots(1,1,figsize=(8,3))
        herdcolors = dict(zip(parcels.SPH_HERD_N,parcels.color))
        gui.plot_moves_timeline(moves,herdcolors,ax)
        tempname2 = tempfile.mktemp()
        fig.savefig(tempname2, format='svg', dpi=150)

    s += '<div class="title"><h3>%s</h3></div>' %title
    s += '<div class="container">'
    s += '<div class="column"> <img  class="scaled" src="%s">  </div>' % tempname
    s += '<div class="column">%s</div>' %treehtml
    s += '</div>'
    s += '<div class="title"><h4>moves timeline</h4></div>'
    s += '<div><img src="%s" class="scaled"></div>' % tempname2
    #tables
    table = '<div> %s</div>' %sub[cols].reset_index().to_html(classes="mytable",index=False)
    s+=table

    with open('report.html', 'w') as f:
        f.write(s)

    #html = HTML(string=s)
    html = HTML(string=open('report.html', 'rb').read(), base_url='./')
    html.write_pdf(outfile, stylesheets=[css_file],
                   optimize_images=True, jpeg_quality=80, dpi=150)
    return
