#!/usr/bin/env python

"""
    Phylogeny viewer.
    Created Feb 2021
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

import sys, os, io, random
import numpy as np
import pandas as pd
import pylab as plt
import string
from .qt import *
from . import tools, widgets
import toyplot, toytree

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
iconpath = os.path.join(module_path, 'plugins', 'icons')
logoimg = os.path.join(iconpath, 'tree.svg')
version = '0.2'

widgetstyle = '''
    QLabel {
        font-size: 14px;
    }
    QPlainTextEdit {
        max-height: 80px;
    }
    QComboBox {
        combobox-popup: 0;
        max-height: 30px;
        max-width: 300px;
    }
    '''

def colormap_from_labels(colormap_name, labels):
    """Get dict of colors mapping to labels using toyplot colormap"""

    n = len(labels)
    #colormap = toyplot.color.brewer.map(colormap_name,domain_min=0,domain_max=n,count=n)
    colormap = toyplot.color.CategoricalMap(toyplot.color.brewer.palette(colormap_name))
    cm = {labels[i]:colormap.css(i) for i in range(n)}
    return cm

def colormap_from_values(colormap_name, vals):
    """Toyplot color mapping from continuous values"""

    labels = vals.unique()
    colormap = toyplot.color.brewer.map(colormap_name, domain_min=vals.min(), domain_max=vals.max(), reverse=True)
    cm = {i:colormap.css(i) for i in labels}
    return cm

def add_legend(canvas, cmap, label='', shape='o'):
    """Category or continuous value legend.
    canvas: canvas
    cmap: mapping of labels to colors
    label: label for legend, optional
    """

    legendmarkers = [(i,toyplot.marker.create(shape=shape, mstyle={"fill":cmap[i]}, size=14)) for i in cmap]
    l = len(legendmarkers)
    s = 50
    min = 15
    e = min+s+l*15
    legend = canvas.legend(
        legendmarkers,
        bounds=("70%", "100%", s, e),
        label=label,
        margin=100,
        )
    return

def add_color_scale(canvas, colormap_name, min=0, max=1, label=''):

    colormap = toyplot.color.brewer.map(name=colormap_name, domain_min=min, domain_max=max, reverse=True)
    legend = canvas.color_scale(colormap,
                                x1="90%", y1="-20%", x2="90%", y2="20%",
                                label=label)
    return

def get_node_colors(tre, df, col, cmap):
    """Node colors from tree and dataframe labels"""

    node_colors=[]
    for n in tre.get_node_values('name', True, True):
        if n in df.index:
            val = df.loc[n][col]
            if val in cmap:
                node_colors.append(cmap[val])
            else:
                node_colors.append('black')
        else:
            node_colors.append('black')
    return node_colors

class PhyloApp(QMainWindow):
    """Tree viewer with Toytree, wraps TreeViewer widget"""
    def __init__(self, filename=None, metafile=None):

        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("Tree Viewer")
        self.setWindowIcon(QIcon(logoimg))
        self.setGeometry(QtCore.QRect(200, 200, 1000, 800))
        self.setMinimumHeight(150)
        if metafile != None:
            meta = pd.read_csv(metafile).set_index('sample')
        else:
            meta = None
        self.main = TreeViewer(self, meta=meta)
        self.main.setFocus()
        self.setCentralWidget(self.main)
        if filename != None:
            self.main.load_tree(filename)
        return

class TreeViewer(QWidget):
    """Phylogeny viewer widget with toytree"""
    def __init__(self, parent=None, filename=None, meta=None):

        super(TreeViewer, self).__init__(parent)
        self.setGeometry(QtCore.QRect(200, 200, 1000, 300))
        self.setMinimumHeight(150)
        self.parent = parent
        self.tree = None
        self.meta = meta
        self.width = 500
        self.height = 800
        self.ts = ''
        import toytree
        self.colors = {}
        self.node_colors = {}
        self.node_sizes = 10
        self.showtiplabels = True
        self.default_style = {
            "layout":'r',
            "edge_type": 'p',
            "edge_widths": 1,
            "edge_style": {
                "stroke": 'black',
                "stroke-width": 2,
            },
            "tip_labels": True,
            "tip_labels_align": False,
            "tip_labels_colors": 'black',
            "tip_labels_style": {
                "font-size": "14px"
            },
            "node_labels": False,
            "node_colors": toytree.colors[2],
            "node_sizes": self.node_sizes,
            "node_markers":"o",
            "node_style": { "stroke": 'black','stroke-width':0},
            "use_edge_lengths":True
        }

        self.style = self.default_style.copy()
        self.add_widgets()
  
        return

    def add_widgets(self):
        """Add widgets"""

        layout = self.layout = QVBoxLayout(self)
        self.main = QSplitter()
        layout.addWidget(self.main)

        self.browser = QWebEngineView()
        self.main.addWidget(self.browser)

        toolswidget = QWidget()
        toolswidget.setMaximumWidth(300)
        self.main.addWidget(toolswidget)
        l = QVBoxLayout(toolswidget)

        self.main.setSizes([400,150])
        self.main.setStretchFactor(1,0)

        w = QLabel('zoom')
        l.addWidget(w)
        self.zoomslider = w = QSlider(QtCore.Qt.Horizontal)
        w.setSingleStep(5)
        w.setMinimum(5)
        w.setMaximum(50)
        w.setValue(10)
        l.addWidget(w)
        w.valueChanged.connect(self.zoom)

        w = QLabel('hscale')
        l.addWidget(w)
        self.widthslider = w = QSlider(QtCore.Qt.Horizontal)
        w.setSingleStep(10)
        w.setMinimum(20)
        w.setMaximum(2000)
        w.setValue(600)
        l.addWidget(w)
        w.valueChanged.connect(self.hscale)

        w = QLabel('vscale')
        l.addWidget(w)
        self.heightslider = w = QSlider(QtCore.Qt.Horizontal)
        w.setSingleStep(10)
        w.setMinimum(20)
        w.setMaximum(1200)
        w.setValue(500)
        l.addWidget(w)
        w.valueChanged.connect(self.vscale)

        return

    def show_tree_menu(self, pos):
        """Show right click tree menu"""

        item = self.tipitems.itemAt( pos )
        menu = QMenu(self.tipitems)
        colorAction = menu.addAction("Set Label Color")
        #nodecolorAction = menu.addAction("Set Node Color")
        rootAction = menu.addAction("Root On")
        dropAction = menu.addAction("Drop Tips")
        action = menu.exec_(self.tipitems.mapToGlobal(pos))
        if action == rootAction:
            self.root_tree()
        elif action == colorAction:
            self.set_color()
        #elif action == nodecolorAction:
        #    self.set_color('node')
        elif action == dropAction:
            self.drop_tips()

    def load_tree(self, filename=None):
        """Load the tree"""

        import toytree
        if filename == None:
            options = QFileDialog.Options()
            filter = "newick files (*.newick);;All files (*.*)"
            filename, _ = QFileDialog.getOpenFileName(self,"Open tree file",
                                        "",filter=filter,selectedFilter =filter, options=options)
            if not filename:
                return
        self.set_tree(toytree.tree(filename))
        self.height = 200+self.tree.ntips*10
        #self.update()
        self.fit_to_window()
        return

    def clear(self):

        self.tree = None
        self.update()
        return

    def update(self):
        """Update the plot"""

        if self.tree == None:
            self.browser.setHtml('')
            self.tipitems.clear()
            return

        canvas = trees.draw_tree(treefile, self.sub, colorcol)
        toyplot.html.render(canvas, "temp.html")
        with open('temp.html', 'r') as f:
            html = f.read()
            self.browser.setHtml(html)
        self.canvas = canvas
        return

    def root_tree(self):

        item = self.tipitems.selectedItems()[0]
        row = self.tipitems.selectedIndexes()[0].row()
        name = item.text(0)
        self.tree = self.tree.root(name).ladderize()
        self.update()
        return

    def unroot_tree(self):
        self.tree = self.tree.unroot().ladderize()
        self.update()
        return

    def export_image(self):
        """Save tree as image"""

        options = QFileDialog.Options()
        filter = "png files (*.png);;pdf files (*.pdf);;All files (*.*)"
        filename, _ = QFileDialog.getSaveFileName(self,"Save Project",
                                    "",filter=filter,selectedFilter =filter, options=options)
        if not filename:
            return

        ext = os.path.splitext(filename)
        print (ext)
        from toyplot import png
        png.render(self.canvas, filename, width=(4, "inches"))
        return

    def fit_to_window(self):
        """Fit tree to frame"""

        self.width = self.browser.geometry().width()-50
        self.height = self.browser.geometry().height()-50
        self.browser.setZoomFactor(1)
        self.update()
        return

    def zoom(self):
        zoom = self.zoomslider.value()/10
        self.browser.setZoomFactor(zoom)

    def hscale(self):
        self.width = self.widthslider.value()
        self.update()
        return

    def vscale(self):
        self.height = self.heightslider.value()
        self.update()
        return
 
    def tree_style_options(self):
        """Tree style dialog"""

        fonts = ['%spx' %i for i in range (6,34)]
        markers = ['o','s','d','^','>']
        nlabels = ['','idx','support']
        tip_labels_style = self.style['tip_labels_style']
        tree_styles = [None,'n','s','c','o','p','d']
        opts = {#'ts': {'type':'combobox','default':self.ts,'items':tree_styles},
                'layout': {'type':'combobox','default':self.style['layout'],'items':['r','d']},
                'edge_type': {'type':'combobox','default':self.style['edge_type'],'items':['p','b','c']},
                'edge_widths': {'type':'spinbox','default':self.style['edge_widths'],'range':(1,5),'interval':1},
                'use_edge_lengths':{'type':'checkbox','default':self.style['use_edge_lengths'] },
                'tip_labels':{'type':'checkbox','default': self.showtiplabels},
                'tip_labels_align':{'type':'checkbox','default':self.style['tip_labels_align'] },
                #'node_labels':{'type':'combobox','default':self.style['node_labels'],'items': nlabels},
                'node_sizes':{'type':'spinbox','default':self.node_sizes,'range':(2,40),'interval':1},
                'node_markers': {'type':'combobox','default':self.style['node_markers'],'items':markers},
                'font_size':{'type':'combobox','default':tip_labels_style['font-size'],'items':fonts}
                }

        dlg = widgets.MultipleInputDialog(self, opts, title='Tree Style', width=300)
        dlg.exec_()
        if not dlg.accepted:
            return False
        kwds = dlg.values
        self.set_style(kwds)
        self.update()
        return

    def set_style(self, kwds):
        """Set style"""

        omit = ['width','height','font_size','ts']
        for k in kwds:
            if k not in omit:
                self.style[k] = kwds[k]
        #if kwds['node_labels'] == '':
        #    self.style['node_labels'] = False
        if 'font_size' in kwds:
            self.style['tip_labels_style']['font-size'] = kwds['font_size']
        #self.ts = kwds['ts']
        self.node_sizes = kwds['node_sizes']
        return

    def reset_style(self):

        self.style = self.default_style
        self.colors = {}
        print (self.style)
        self.update()


    def drop_tips(self):

        items = self.tipitems.selectedItems()
        names = [i.text(0) for i in items]
        self.tree = self.tree.drop_tips(names=names).ladderize()
        self.update()
        return


def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='snipgenie gui tool')
    parser.add_argument("-f", "--file", dest="filename",default=[],
                        help="input tree file", metavar="FILE")
    parser.add_argument("-m", "--meta", dest="metafile",default=None,
                        help="meta data csv file", metavar="FILE")
    args = vars(parser.parse_args())
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_ShareOpenGLContexts)
    app = QApplication(sys.argv)
    aw = PhyloApp(**args)
    aw.show()
    app.exec_()

if __name__ == '__main__':
    main()
