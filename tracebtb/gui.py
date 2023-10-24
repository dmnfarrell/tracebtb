#!/usr/bin/env python

"""
    tracebtb GUI.
    Created Mar 2022
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

import sys,os,traceback,subprocess
import glob,platform,shutil
import pickle
import threading,time
import math
from .qt import *
import pandas as pd
import numpy as np
import pylab as plt
from Bio import SeqIO, AlignIO
import matplotlib as mpl
from . import core, widgets, webwidgets, tables, tools, plotting, treeview, trees
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from matplotlib_scalebar.scalebar import ScaleBar
import contextily as cx

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
data_path = os.path.join(module_path,'data')
logoimg = os.path.join(module_path, 'logo.svg')
iconpath = os.path.join(module_path, 'icons')
#settingspath = os.path.join(homepath, '.config','tracebtb')

counties_gdf = gpd.read_file(os.path.join(data_path,'counties.shp')).to_crs("EPSG:29902")
counties = ['Clare','Cork','Cavan','Monaghan','Louth','Kerry','Meath','Wicklow']
cladelevels = ['snp3','snp12','snp50','snp200','snp500']
providers = {'None':None,
            'OSM':cx.providers.OpenStreetMap.Mapnik,
            'CartoDB':cx.providers.CartoDB.Positron,
            'CartoDB.Voyager': cx.providers.CartoDB.Voyager}
#colormaps = sorted(m for m in plt.cm.datad if not m.endswith("_r"))
colormaps = ['Paired', 'Dark2', 'Set1', 'Set2', 'Set3',
            'tab10', 'tab20', 'tab20b', 'tab20c',
            'twilight', 'twilight_shifted', 'hsv',
            'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg',
            'gist_rainbow', 'rainbow', 'jet', 'turbo', 'nipy_spectral']
style = '''
    QWidget {
        max-width: 130px;
        font-size: 12px;
    }
    QPlainTextEdit {
        max-height: 80px;
    }
    QScrollBar:vertical {
         width: 15px;
         margin: 1px 0 1px 0;
     }
    QScrollBar::handle:vertical {
         min-height: 20px;
     }
    QComboBox {
        combobox-popup: 0;
        max-height: 30px;
        max-width: 120px;
    }
    '''

dockstyle = '''
    QDockWidget {
        max-width:1000px;
    }
    QDockWidget::title {
        background-color: #80bfff;
    }
    QScrollBar:vertical {
         width: 15px;
         margin: 1px 0 1px 0;
     }
    QScrollBar::handle:vertical {
         min-height: 20px;
     }
'''

#module level functions

def show_labels(df, col, ax):
    """Add labels to plot"""

    if col == '': return
    for x, y, label in zip(df.geometry.x, df.geometry.y, df[col]):
        ax.annotate(label, xy=(x, y), xytext=(5, 0), textcoords="offset points",
                    fontsize=8)
    return

def get_largest_poly(x):
    if type(x) is MultiPolygon:
        return max(x.geoms, key=lambda a: a.area)
    else:
        return x

def calculate_parcel_centroids(parcels):
    """Get centroids of lpis parcels"""

    largest = parcels.geometry.apply(get_largest_poly)
    cent = largest.geometry.centroid
    cent = gpd.GeoDataFrame(geometry=cent,crs='EPSG:29902')
    cent['SPH_HERD_N'] = parcels.SPH_HERD_N
    return cent

def plot_single_cluster(df, col=None, cmap=None, margin=None, ms=40,
                        legend=False, title='', ax=None):
    """Plot a single map view of a set of points/farms"""

    df = df[~df.is_empty]
    df = df[df.geometry.notnull()]

    if len(df) == 0:
        ax.clear()
        ax.set_title('no locations available')
        return
    minx, miny, maxx, maxy = df.total_bounds
    maxx+=10
    maxy+=10
    df = df[df.Species.isin(['Bovine','Badger'])]
    df['color'] = df.Species.map({'Bovine':'blue','Badger':'orange'})

    if col == None or col == '':
        df.plot(color=df.color,ax=ax,alpha=0.6,markersize=ms,linewidth=.5,label='farm/badger',legend=legend)
    else:
        cow = df[df.Species=='Bovine']
        badger = df[df.Species=='Badger']
        cow.plot(column=col,ax=ax,alpha=0.6,markersize=ms,linewidth=1,ec='black',cmap=cmap,legend=legend)
        badger.plot(color='orange',ax=ax,alpha=0.6,marker='s',markersize=ms,linewidth=1,ec='black')
    ax.set_title(title)
    ax.axis('off')
    if margin == None:
        margin = (maxx-minx)*0.3
    ax.set_xlim(minx-margin,maxx+margin)
    ax.set_ylim(miny-margin,maxy+margin)
    ax.add_artist(ScaleBar(dx=1,location=3))
    #ax.legend(fontsize=12)
    return

def calculate_grid_dimensions(n):
    """Calculate the number of rows and columns that best approximate a square layout"""

    sqrt_n = math.ceil(n**0.5)
    rows = math.ceil(n / sqrt_n)
    columns = math.ceil(n / rows)
    return rows, columns

def plot_grid(gdf, parcels=None, moves=None, fig=None, source=None):
    """
    Plot herds on multiple axes. Useful for overview of farms and reports.
    """

    groups = gdf.groupby('HERD_NO')
    n = len(groups)+1
    rows, cols = calculate_grid_dimensions(n)
    if fig == None:
        fig,ax = plt.subplots(rows,rows,figsize=(8,8))
        axs = ax.flat
    else:
        axs = widgets.add_subplots_to_figure(fig, rows, cols)
        print (axs)
    i=0
    for idx,df in groups:
        ax=axs[i]
        r=df.iloc[0]
        if parcels is not None:
            subp = parcels[parcels.SPH_HERD_N==idx]
            subp.plot(color='red',alpha=0.6,lw=1,ax=ax)
        df.plot(color='black',ax=ax)            
        ax.axis('off')
        ax.add_artist(ScaleBar(dx=1,location=3))
        ax.set_title(r.HERD_NO)       
        if source != None:     
            cx.add_basemap(ax, crs=df.crs, attribution=False, source=source)
        i+=1
    ax=axs[i]
    gdf.plot(ax=ax)
    counties_gdf.plot(color='none',lw=.2,ax=ax)
    if moves is not None:
        plot_moves(moves,lpis_cent,ax)
    ax.axis('off')
    plt.tight_layout()
    return ax

def plot_moves(moves, lpis_cent, ax):
    """Show moves as lines on plot"""

    colors = plotting.random_colors(250, seed=12)
    i=0
    if moves is None:
        return
    moves = moves[moves.geometry.notnull()]
    for tag,t in moves.groupby('Animal_ID'):
        if t is not None:
            #print (t[cols])
            moved = lpis_cent[lpis_cent.SPH_HERD_N.isin(t.move_to)]
            coords = tools.get_coords_data(t)
            if len(coords)>0:
                mlines = gpd.GeoDataFrame(geometry=coords)
                mlines.plot(color='black',linewidth=.5,ax=ax)
                moved.plot(color='none',ec='black',marker='s',
                            markersize=80,linewidth=.8,alpha=0.5,ax=ax)
                i+=1
    return

def plot_moves_timeline(df, ax=None):
    """Timeline from moves"""

    from datetime import datetime, timedelta
    from matplotlib.patches import Rectangle
    import matplotlib.dates as mdates

    df['move_date'] = pd.to_datetime(df.move_date)
    groups = df.groupby('Animal_ID')
    if ax == None:
        fig,ax=plt.subplots(1,1,figsize=(10,4))
    i=.1
    tags = groups.groups.keys()
    clrs,cmap = plotting.get_color_mapping(df, 'move_to')
    leg = {}
    for tag,t in groups:
        if t is None:
            continue
        d = get_move_dates(t)
        for r,row in d.iterrows():
            herd = row.move_to
            start = mdates.date2num(row.move_date)
            end = mdates.date2num(row.end_date)
            width = end - start
            rect = Rectangle((start, i), width, .8, color=cmap[herd], ec='black')
            ax.add_patch(rect)
            leg[herd] = rect
        i+=1

    ax.set_xlim((14000, 19000))
    ax.set_ylim((0, i))
    locator = mdates.AutoDateLocator(minticks=3)
    formatter = mdates.AutoDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)

    ax.set_yticks(np.arange(len(tags))+0.5)
    ax.set_yticklabels(tags)
    ax.grid(axis='y')
    plt.subplots_adjust(left=0.3)
    #ax.legend(leg.values(), leg.keys(), loc='upper right', bbox_to_anchor=(1.2, 1))

    ax.tick_params(axis='both', labelsize=8)
    return

def jitter_points(r, scale=1):
    """Jitter GeoDataFrame points, vector based function"""

    a=np.random.normal(0,scale)
    b=np.random.normal(0,scale)
    if (r.geometry.is_empty): return r.geometry
    x,y = r.geometry.x+a,r.geometry.y+b
    return Point(x,y)

def get_moves_bytag(df, move_df, lpis_cent):
    """Get moves and coords for one or more samples.
    """

    cols=['Animal_ID','HERD_NO','move_to','move_date','data_type','breed','dob']
    t = df.merge(move_df,left_on='Animal_ID',right_on='tag',how='inner')[cols]
    m = t.merge(lpis_cent,left_on='move_to',right_on='SPH_HERD_N', how='left')

    if len(m)==0:
        return
    x = lpis_cent[lpis_cent.SPH_HERD_N.isin(df.HERD_NO)]
    m = pd.concat([m,x]).dropna(subset='Animal_ID')
    m = m.sort_values(by=['Animal_ID','move_date'])
    m = gpd.GeoDataFrame(m)
    return m

def get_move_dates(df):

    df['end_date'] = df.move_date.shift(-1)
    return df[['move_date','end_date','move_to']][:-1]

class CustomTreeWidgetItem( QTreeWidgetItem ):
    def __init__(self, parent=None):
        QTreeWidgetItem.__init__(self, parent)

    def __lt__(self, otherItem):
        column = self.treeWidget().sortColumn()
        try:
            return float( self.text(column) ) < float( otherItem.text(column) )
        except ValueError:
            return self.text(column) > otherItem.text(column)

class App(QMainWindow):
    """GUI Application using PySide2 widgets"""
    def __init__(self, project=None):

        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("TracebTB")

        self.setWindowIcon(QIcon(logoimg))
        self.create_menu()
        self.main = QSplitter(self)
        screen_resolution = QGuiApplication.primaryScreen().availableGeometry()
        width, height = screen_resolution.width()*0.9, screen_resolution.height()*.8
        if screen_resolution.width()>1920:
            self.setGeometry(QtCore.QRect(100, 100, int(width), int(height)))
        else:
            self.showMaximized()
        self.setMinimumSize(400,300)

        self.recent_files = ['']
        self.scratch_items = {}
        self.opentables = {}
        self.lpis_master = None
        self.lpis_cent = None
        self.parcels = None
        self.neighbours = None
        self.sub = None
        self.main.setFocus()
        self.setCentralWidget(self.main)
        self.create_tool_bar()
        self.load_settings()
        self.setup_gui()
        self.show_recent_files()

        self.load_base_data()
        self.new_project()
        self.running = False
        self.title = None
        #self.load_test()

        if project != None:
            self.load_project(project)
        self.threadpool = QtCore.QThreadPool()
        self.redirect_stdout()
        return

    def redirect_stdout(self):
        """redirect stdout"""
        self._stdout = StdoutRedirect()
        self._stdout.start()
        self._stdout.printOccur.connect(lambda x : self.info.insert(x))
        return

    def load_settings(self):
        """Load GUI settings"""

        s = self.settings = QtCore.QSettings('tracebtb','default')
        try:
            winsize = s.value('window_size')
            if winsize != None:
                self.resize(s.value('window_size'))
                self.move(s.value('window_position'))
                core.FONT = s.value("font")
                core.FONTSIZE = int(s.value("fontsize"))
                core.DPI = int(s.value("dpi"))
                import matplotlib as mpl
                mpl.rcParams['savefig.dpi'] = int(core.DPI)
                core.ICONSIZE = int(s.value("iconsize"))
                self.setIconSize(QtCore.QSize(core.ICONSIZE, core.ICONSIZE))
                r = s.value("recent_files")
                if r != '':
                    rct = r.split(',')
                    self.recent_files = [f for f in rct if os.path.exists(f)]

        except Exception as e:
            print (e)
        return

    def apply_settings(self):
        """Apply settings to GUI when changed"""

        self.setIconSize(QtCore.QSize(core.ICONSIZE, core.ICONSIZE))
        for i in self.opentables:
            table = self.opentables[i]
            table.fontname = core.FONT
            table.fontsize = core.FONTSIZE
            table.updateFont()
            table.refresh()
        import matplotlib as mpl
        mpl.rcParams['savefig.dpi'] = core.DPI
        return

    def save_settings(self):
        """Save GUI settings"""

        self.settings.setValue('window_size', self.size())
        self.settings.setValue('window_position', self.pos())
        self.settings.setValue('iconsize', core.ICONSIZE)
        self.settings.setValue('font', core.FONT)
        self.settings.setValue('fontsize', core.FONTSIZE)
        self.settings.setValue('dpi', core.DPI)
        self.settings.setValue('recent_files',','.join(self.recent_files))
        #print (self.settings)
        self.settings.sync()
        return

    def load_base_data(self):
        #reference map of counties
        self.counties = counties_gdf#gpd.read_file(os.path.join(data_path,'counties.shp')).to_crs("EPSG:29902")
        return

    def create_tool_bar(self):
        """Create main toolbar"""

        items = {'New project': {'action': lambda: self.new_project(ask=True),'file':'document-new'},
                 'Open': {'action':self.load_project,'file':'document-open'},
                 'Save': {'action': lambda: self.save_project(),'file':'save'},
                 'Zoom out': {'action':self.zoom_out,'file':'zoom-out'},
                 'Zoom in': {'action':self.zoom_in,'file':'zoom-in'},
                 'Scratchpad': {'action':self.show_scratchpad,'file':'scratchpad'},
                 'Settings': {'action':self.preferences,'file':'settings'},
                 'Herd Summary':{'action':self.herd_summary,'file':'cow'},
                 'Cluster Report':{'action':self.cluster_report ,'file':'cluster_report'},
                 #'Make Simulated Data':{'action':self.make_test_data ,'file':'simulate'},
                 'Quit': {'action':self.quit,'file':'application-exit'}
                }

        toolbar = QToolBar("Main Toolbar")
        self.addToolBar(toolbar)
        #toolbar.setOrientation(QtCore.Qt.Vertical)
        widgets.addToolBarItems(toolbar, self, items)
        return

    def add_dock(self, widget, name, side='left', scrollarea=True):
        """Add a dock widget"""

        dock = QDockWidget(name)
        dock.setStyleSheet(dockstyle)
        if scrollarea == True:
            area = QScrollArea()
            area.setWidgetResizable(True)
            dock.setWidget(area)
            area.setWidget(widget)
        else:
            dock.setWidget(widget)
        if side == 'left':
            self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)
        else:
            self.addDockWidget(QtCore.Qt.RightDockWidgetArea, dock)
        if side == 'floating':
            dock.setFloating(True)
        self.docks[name] = dock
        dock.setSizePolicy( QSizePolicy.Expanding, QSizePolicy.Minimum )
        return dock

    def create_option_widgets(self):
        """Set up map view options"""

        self.widgets = {}
        m = QWidget()
        m.setStyleSheet(style)
        m.setMaximumWidth(140)
        m.setMinimumWidth(140)
        l = QVBoxLayout()
        l.setAlignment(QtCore.Qt.AlignTop)
        m.setLayout(l)

        #select cluster level
        self.cladelevelw = w = QComboBox(m)
        w.addItems(cladelevels)
        l.addWidget(QLabel('Clade level:'))
        l.addWidget(w)
        w.currentIndexChanged.connect(self.update_clades)
        self.widgets['cladelevel'] = w
        #select clade/cluster
        l.addWidget(QLabel('Clade:'))
        t = self.cladew = QTreeWidget()
        t.setHeaderItem(CustomTreeWidgetItem(["name","size"]))
        t.setColumnWidth(0, 50)
        t.setColumnWidth(1, 30)
        t.setSortingEnabled(True)
        t.setMinimumHeight(30)
        #t.setSelectionMode(QAbstractItemView.ExtendedSelection)
        #t.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        #t.customContextMenuRequested.connect(self.update_clades)
        t.itemSelectionChanged.connect(self.plot_selected_clade)
        l.addWidget(t)
        self.widgets['clade'] = t

        #zoom to county
        self.countyw =w = QComboBox(m)
        w.addItems(counties)
        l.addWidget(QLabel('County:'))
        l.addWidget(w)
        w.currentIndexChanged.connect(self.plot_county)
        #labels
        self.labelsw = w = QComboBox(m)
        l.addWidget(QLabel('Labels:'))
        l.addWidget(w)
        self.widgets['labels'] = w
        #color points by
        self.colorbyw = w = QComboBox(m)
        l.addWidget(QLabel('Color by:'))
        l.addWidget(w)
        w.setMaxVisibleItems(12)
        self.widgets['colorby'] = w
        #colormaps
        self.cmapw = w = QComboBox(m)
        l.addWidget(QLabel('Colormap:'))
        l.addWidget(w)
        w.addItems(colormaps)
        w.setCurrentText('Paired')
        self.widgets['colormap'] = w
        #toggle cx
        self.contextw = w = QComboBox(m)
        l.addWidget(QLabel('Context map:'))
        l.addWidget(w)
        w.addItems(providers.keys())
        self.widgets['context'] = w
        self.markersizew = w = QSpinBox(m)
        w.setRange(1,300)
        w.setValue(50)
        l.addWidget(QLabel('Marker size:'))
        l.addWidget(w)
        l.addStretch()
        self.widgets['markersize'] = w
        return m

    def map_buttons(self):
        """Create map buttons"""

        m = QWidget()
        m.setMaximumWidth(60)
        m.setMinimumWidth(50)
        l = QVBoxLayout()
        l.setAlignment(QtCore.Qt.AlignTop)
        m.setLayout(l)
        b = widgets.createButton(m, None, self.update, 'refresh', core.ICONSIZE)
        l.addWidget(b)
        self.showfoliumb = b = widgets.createButton(m, None, self.show_folium, 'folium', core.ICONSIZE, 'interactive view')
        #b.setCheckable(True)
        l.addWidget(b)
        self.multiaxesb= b = widgets.createButton(m, None, self.plot_farms, 'plot-grid', core.ICONSIZE, 'farm split view')
        l.addWidget(b)
        b = widgets.createButton(m, None, self.plot_in_region, 'plot-region', core.ICONSIZE, 'show all in region')
        l.addWidget(b)
        self.parcelsb = b = widgets.createButton(m, None, self.update, 'parcels', core.ICONSIZE, 'show parcels')
        b.setCheckable(True)
        l.addWidget(b)
        self.widgets['showparcels'] = b
        self.movesb = b = widgets.createButton(m, None, self.update, 'plot-moves', core.ICONSIZE, 'show moves')
        b.setCheckable(True)
        l.addWidget(b)
        self.widgets['showmoves'] = b
        self.timelineb = b = widgets.createButton(m, None, self.update, 'timeline', core.ICONSIZE, 'moves timeline')
        b.setCheckable(True)
        l.addWidget(b)
        self.widgets['showtimeline'] = b
        self.legendb = b = widgets.createButton(m, None, self.update, 'legend', core.ICONSIZE, 'show legend')
        b.setCheckable(True)
        l.addWidget(b)
        self.widgets['showlegend'] = b
        self.colorcountiesb = b = widgets.createButton(m, None, self.update, 'counties', core.ICONSIZE, 'color counties')
        b.setCheckable(True)
        l.addWidget(b)
        self.jitterb = b = widgets.createButton(m, None, self.update, 'jitter', core.ICONSIZE, 'jitter points')
        b.setCheckable(True)
        l.addWidget(b)
        self.showneighboursb = b = widgets.createButton(m, None, self.update, 'neighbours', core.ICONSIZE, 'show neighbours')
        b.setCheckable(True)
        l.addWidget(b)
        self.widgets['showneighbours'] = b
        self.showtreeb = b = widgets.createButton(m, None, self.update, 'tree', core.ICONSIZE, 'show tree')
        b.setCheckable(True)
        l.addWidget(b)
        self.widgets['showtree'] = b
        self.showmstb = b = widgets.createButton(m, None, self.update, 'mst', core.ICONSIZE, 'show MST')
        b.setCheckable(True)
        l.addWidget(b)
        b = widgets.createButton(m, None, self.save_to_scratchpad, 'snapshot', core.ICONSIZE, 'take snapshot')
        l.addWidget(b)
        return m

    def update_clades(self):
        """Update the clade tree widget"""

        level = self.cladelevelw.currentText()
        clades = self.cent[level].value_counts()
        clades = clades[clades>1]
        t = self.cladew
        t.clear()

        for cl,size in clades.items():
            item = QTreeWidgetItem(t)
            item.setText(0, cl)
            item.setText(1, str(size))
        return

    def setup_gui(self):
        """Add all GUI elements"""

        self.docks = {}
        self.main = main = QWidget()
        self.m = QSplitter()

        self.main.setFocus()
        self.setCentralWidget(self.main)
        l = QVBoxLayout(main)
        l.addWidget(self.m)

        self.meta_table = tables.SampleTable(self, dataframe=pd.DataFrame(),
                                           font=core.FONT, fontsize=core.FONTSIZE, app=self)
        t = self.table_widget = tables.DataFrameWidget(parent=self, table=self.meta_table,
                                        toolbar=False)
        self.add_dock(self.table_widget, 'meta data', scrollarea=False)
        self.add_dock_item('meta data')
        #self.m.addWidget(self.table_widget)
        self.opentables['main'] = self.meta_table

        w = self.create_option_widgets()
        self.m.addWidget(w)
        w = self.map_buttons()
        self.m.addWidget(w)

        self.tabs = QTabWidget(main)
        self.m.addWidget(self.tabs)
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self.close_tab)

        self.plotview = widgets.CustomPlotViewer(self, controls=False, app=self)
        #self.m.addWidget(self.plotview)
        idx = self.tabs.addTab(self.plotview, 'Map')
        self.tabs.setCurrentIndex(idx)

        self.treeview = treeview.TreeViewer()
        #self.m.addWidget(self.treeview)

        self.info = widgets.Editor(main, readOnly=True, fontsize=10)
        self.add_dock(self.info, 'log', 'right')
        self.foliumview = webwidgets.FoliumViewer(main)
        #self.foliumview.show()
        idx = self.tabs.addTab(self.foliumview, 'Interactive')
        #test map
        #tmp = webwidgets.FoliumViewer(main)
        #idx = self.tabs.addTab(tmp, 'test')
        #tmp.test_map()

        self.info.append("Welcome\n")
        self.statusBar = QStatusBar()
        self.projectlabel = QLabel('')
        self.projectlabel.setStyleSheet('color: blue')
        self.projectlabel.setAlignment(Qt.AlignLeft)
        self.statusBar.addWidget(self.projectlabel, 1)

        self.progressbar = QProgressBar()
        self.progressbar.setRange(0,1)
        self.statusBar.addWidget(self.progressbar, 3)
        self.progressbar.setAlignment(Qt.AlignRight)
        self.setStatusBar(self.statusBar)

        #redirect stdout
        #self._stdout = StdoutRedirect()
        #self._stdout.start()
        #self._stdout.printOccur.connect(lambda x : self.info.insert(x))

        #add dock menu items
        for name in ['log']:
            action = self.docks[name].toggleViewAction()
            self.dock_menu.addAction(action)
            action.setCheckable(True)
        return

    def add_dock_item(self, name):

        action = self.docks[name].toggleViewAction()
        self.dock_menu.addAction(action)
        action.setCheckable(True)
        return

    @QtCore.Slot(int)
    def close_tab(self, index):
        """Close current tab"""

        index = self.tabs.currentIndex()
        name = self.tabs.tabText(index)
        if name != 'Map':
            self.tabs.removeTab(index)
        return

    def create_menu(self):
        """Create the menu bar for the application. """

        self.file_menu = QMenu('File', self)
        #self.file_menu.addAction('Load Folder', lambda: self.load_folder())
        self.file_menu.addAction('Load Samples', lambda: self.load_samples())
        #icon = QIcon(os.path.join(iconpath,'shapefile.svg'))
        #self.tools_menu.addAction(icon,'Load Shapefile', self.load_shapefile)
        self.file_menu.addAction('Load Moves', lambda: self.load_moves())
        self.file_menu.addAction('Load Alignment', lambda: self.load_alignment())
        icon = QIcon(os.path.join(iconpath,'document-new.svg'))
        self.file_menu.addAction(icon, 'New Project', lambda: self.new_project(ask=True))
        icon = QIcon(os.path.join(iconpath,'document-open.svg'))
        self.file_menu.addAction(icon, 'Open Project', self.load_project_dialog)
        self.recent_files_menu = QMenu("Recent Projects", self.file_menu)
        self.file_menu.addAction(self.recent_files_menu.menuAction())
        icon = QIcon(os.path.join(iconpath,'save.svg'))
        self.file_menu.addAction(icon, '&Save Project', self.save_project,
                QtCore.Qt.CTRL + QtCore.Qt.Key_S)
        self.file_menu.addAction('Save Project As', self.save_project_dialog)
        icon = QIcon(os.path.join(iconpath,'application-exit.svg'))
        self.file_menu.addAction(icon, 'Quit', self.quit)
        self.menuBar().addMenu(self.file_menu)

        self.edit_menu = QMenu('Edit', self)
        self.menuBar().addMenu(self.edit_menu)
        icon = QIcon(os.path.join(iconpath,'settings.svg'))
        self.edit_menu.addAction(icon, 'Preferences', self.preferences)

        self.view_menu = QMenu('View', self)
        self.menuBar().addMenu(self.view_menu)
        icon = QIcon(os.path.join(iconpath,'zoom-in.svg'))
        self.view_menu.addAction(icon, 'Zoom In', self.zoom_in,
                QtCore.Qt.CTRL + QtCore.Qt.Key_Equal)
        icon = QIcon(os.path.join(iconpath,'zoom-out.svg'))
        self.view_menu.addAction(icon, 'Zoom Out', self.zoom_out,
                QtCore.Qt.CTRL + QtCore.Qt.Key_Minus)

        self.tools_menu = QMenu('Tools', self)
        self.menuBar().addMenu(self.tools_menu)
        icon = QIcon(os.path.join(iconpath,'shapefile.svg'))
        self.tools_menu.addAction(icon, 'Load LPIS file', self.set_lpis_file)
        icon = QIcon(os.path.join(iconpath,'parcels.svg'))
        self.tools_menu.addAction(icon, 'Extract LPIS data', self.get_lpis_centroids)
        icon = QIcon(os.path.join(iconpath,'neighbours.svg'))
        self.tools_menu.addAction(icon,'Extract neighbouring parcels', self.get_neighbouring_parcels)

        self.tools_menu.addSeparator()
        icon = QIcon(os.path.join(iconpath,'mbovis.svg'))
        self.tools_menu.addAction(icon, 'Strain Typing', self.strain_typing)
        icon = QIcon(os.path.join(iconpath,'cow.svg'))
        self.tools_menu.addAction(icon, 'Show Herd Summary', self.herd_summary)
        icon = QIcon(os.path.join(iconpath,'simulate.svg'))
        self.tools_menu.addAction(icon, 'Make Simulated Data', self.make_test_data)
        icon = QIcon(os.path.join(iconpath,'pdf.svg'))
        self.tools_menu.addAction(icon,'Cluster Report', self.cluster_report)

        self.scratch_menu = QMenu('Scratchpad', self)
        self.menuBar().addMenu(self.scratch_menu)
        icon = QIcon(os.path.join(iconpath,'scratchpad.svg'))
        self.scratch_menu.addAction(icon,'Show Scratchpad', lambda: self.show_scratchpad())
        icon = QIcon(os.path.join(iconpath,'snapshot.svg'))
        self.scratch_menu.addAction(icon,'Plot to Scratchpad', lambda: self.save_to_scratchpad())

        self.dock_menu = QMenu('Docks', self)
        self.menuBar().addMenu(self.dock_menu)

        self.help_menu = QMenu('Help', self)
        self.menuBar().addMenu(self.help_menu)
        self.help_menu.addAction('&Help', self.online_documentation)
        self.help_menu.addAction('About', self.about)

    def show_recent_files(self):
        """Populate recent files menu"""

        from functools import partial
        if self.recent_files == None:
            return
        for fname in self.recent_files:
            self.recent_files_menu.addAction(fname, partial(self.load_project, fname))
        self.recent_files_menu.setEnabled(len(self.recent_files))
        return

    def add_recent_file(self, fname):
        """Add file to recent if not present"""

        fname = os.path.abspath(fname)
        if fname and fname not in self.recent_files:
            self.recent_files.insert(0, fname)
            if len(self.recent_files) > 5:
                self.recent_files.pop()
        self.recent_files_menu.setEnabled(len(self.recent_files))
        return

    def update_dock_items(self, data):
        """Update dock status from saved menu items"""

        for action in self.dock_menu.actions():
            name = action.text()
            if name in data:
                val = data[name]
                if val is False:
                    self.docks[name].hide()
        return

    def save_project(self):
        """Save project"""

        if self.proj_file == None:
            self.save_project_dialog()

        filename = self.proj_file
        data={}
        keys = ['cent','sub','moves','parcels','lpis_cent','neighbours','aln']
        for k in keys:
            if hasattr(self, k):
                data[k] = self.__dict__[k]
        data['scratch_items'] = self.scratch_items
        data['fig'] = self.plotview.fig
        data['widget_values'] = widgets.getWidgetValues(self.widgets)
        dock_items = {}
        for action in self.dock_menu.actions():
            dock_items[action.text()] = action.isChecked()
        data['dock_items'] = dock_items
        pickle.dump(data, open(filename,'wb'))
        self.add_recent_file(filename)
        return

    def save_project_dialog(self):
        """Save as project"""

        options = QFileDialog.Options()
        filename, _ = QFileDialog.getSaveFileName(self,"Save Project",
                                                  "","Project files (*.tracebtb);;All files (*.*)",
                                                  options=options)
        if filename:
            if not os.path.splitext(filename)[1] == '.tracebtb':
                filename += '.tracebtb'
            self.proj_file = filename
            self.save_project()
            self.projectlabel.setText(filename)
        return

    def new_project(self, ask=False):
        """Clear all loaded inputs and results"""

        reply=None
        if ask == True:
            reply = QMessageBox.question(self, 'Confirm',
                                "This will clear the current project.\nAre you sure?",
                                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.No:
            return False

        self.outputdir = None
        self.proj_file = None
        self.meta_table.setDataFrame(pd.DataFrame({'sample':[]}))
        self.plotview.clear()
        self.cladew.clear()
        for i in self.opentables:
            w = self.opentables[i]
            w.setDataFrame()
        #if hasattr(self, 'treeview'):
        #    self.treeview.clear()
        return

    def load_project(self, filename=None):
        """Load project"""

        self.new_project()
        data = pickle.load(open(filename,'rb'))
        keys = ['cent','sub','moves','parcels','lpis_cent','neighbours','aln']
        for k in keys:
            if k in data:
                self.__dict__[k] = data[k]

        t = self.meta_table
        t.setDataFrame(self.cent)
        self.table_widget.updateStatusBar()
        self.update_clades()
        cols = ['']+list(self.cent.columns)
        self.labelsw.addItems(cols)
        self.colorbyw.addItems(cols)
        self.colorbyw.setCurrentText('')
        if 'widget_values' in data:
            #print (data['widget_values'])
            widgets.setWidgetValues(self.widgets, data['widget_values'])
        self.proj_file = filename
        self.projectlabel.setText(self.proj_file)
        if 'scratch_items' in data:
            self.scratch_items = data['scratch_items']
        self.add_recent_file(filename)
        self.update()
        if 'dock_items' in data:
            self.update_dock_items(data['dock_items'])
        #t = self.cladew
        #t.setCurrentIndex(t.model().index(0, 4))
        #print (t.selectedIndexes())
        #self.cent = self.cent.set_index('sample')
        return

    def load_project_dialog(self):
        """Load project"""

        filename, _ = QFileDialog.getOpenFileName(self, 'Open Project', './',
                                        filter="Project Files(*.tracebtb);;All Files(*.*)")
        if not filename:
            return
        if not os.path.exists(filename):
            print ('no such file')
        self.load_project(filename)
        return

    def load_shapefile(self):

        return

    def load_test(self):
        """Load test dataset"""

        #metadata
        df = pd.read_csv('testing/metadata.csv')
        #index_col = 'sample'
        #df.set_index(index_col,inplace=True)
        self.cent = df
        t = self.meta_table
        t.setDataFrame(df)
        t.resizeColumns()

        for col in cladelevels:
            self.cent[col] = self.cent[col].astype(str)
        self.update_clades()
        #print (self.cent[:5])

        cols = ['']+list(self.cent.columns)
        self.labelsw.addItems(cols)
        self.colorbyw.addItems(cols)
        self.colorbyw.setCurrentText('')
        self.plot_selected_clade()

        #snps
        #self.coresnps = pd.read_csv('testing/core_snps_mbovis.txt', sep=' ')
        #moves
        self.moves= pd.read_csv('testing/moves.csv')
        return

    def set_lpis_file(self):
        """Set path to LPIS master file"""

        filename, _ = QFileDialog.getOpenFileName(self, 'Open File', './',
                                        filter="Shapefiles(*.shp);;All Files(*.*)")
        if not filename:
            return
        print ('reading LPIS file..')
        def func(progress_callback):
            self.lpis_master = gpd.read_file(filename).set_crs('EPSG:29902')
        self.run_threaded_process(func, self.processing_completed)
        return

    def get_lpis_centroids(self):
        """Calculate LPIS centroids and parcels for samples, run once when metdata
        is updated."""

        if self.lpis_master is None:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("No LPIS master file loaded")
            msg.exec_()
            return

        def completed():
            #add locations to meta data
            print ('updating table..')
            self.set_locations()
            x = self.lpis_master
            self.parcels = x[x.SPH_HERD_N.isin(self.cent.HERD_NO)]
            self.processing_completed()
            return

        def func(progress_callback):
            self.lpis_cent = calculate_parcel_centroids(self.lpis_master)
        print ('calculating centroids..')
        self.run_threaded_process(func, completed)
        return

    def get_neighbouring_parcels(self):
        """Find all neighbouring parcels from LPIS.
        Requires the LPIS master file."""

        if self.lpis_master is None:
            print ('LPIS master file not loaded')
            return

        def completed():
            self.processing_completed()
            self.update()
            print ('found %s parcels' %len(self.neighbours))
            return

        def func(progress_callback):
            d=800
            found = []
            lpis = self.lpis_master
            for x in self.cent.geometry:
                dists = self.lpis_cent.distance(x)
                points = self.lpis_cent[(dists<=d) & (dists>10)]
                found.append(points)

            found = pd.concat(found).drop_duplicates()
            p = lpis[lpis.SPH_HERD_N.isin(found.SPH_HERD_N)]
            self.neighbours = p

        print ('getting parcels..')
        self.run_threaded_process(func, completed)
        return

    def set_locations(self):
        """Add locations to meta data from lpis_cent."""

        if self.lpis_cent is None:
            print ('no centroids')
            return
        df = self.cent
        #remove previous geometry if any
        if type(df) is gpd.GeoDataFrame:
            df = pd.DataFrame(df)
            df.drop(columns=['geometry'], inplace=True)

        df = df.merge(self.lpis_cent,left_on='HERD_NO',right_on='SPH_HERD_N',how='left')#.set_index(idx)
        df = gpd.GeoDataFrame(df,geometry=df.geometry).set_crs('EPSG:29902')
        df['geometry'] = [Point() if x is None else x for x in df.geometry]
        #if 'sample' in df.columns:
        #    df = df.set_index('sample')
        self.cent = df
        self.meta_table.setDataFrame(self.cent)
        return

    def load_samples(self):
        """Load meta data"""

        reply = QMessageBox.question(self, 'Continue?', "This will overwrite the current meta data. Are you sure?",
                                        QMessageBox.Cancel | QMessageBox.Yes)
        if reply == QMessageBox.Cancel:
            return

        filename, _ = QFileDialog.getOpenFileName(self, 'Open Meta Data', './',
                                        filter="csv file(*.csv *.txt);;All Files(*.*)")
        if not filename:
            return
        self.cent = pd.read_csv(filename)
        for col in cladelevels:
            self.cent[col] = self.cent[col].astype(str)
        t = self.meta_table
        t.setDataFrame(self.cent)
        self.update_clades()
        return

    def load_moves(self):
        """Load movement data"""

        filename, _ = QFileDialog.getOpenFileName(self, 'Open CSV file', './',
                                        filter="csv file(*.csv *.txt);;All Files(*.*)")
        if not filename:
            return
        self.moves = pd.read_csv(filename)
        print ('loaded %s rows' %len(self.moves))
        return

    def load_alignment(self):
        """Load sequence alignment"""

        filename, _ = QFileDialog.getOpenFileName(self, 'Open Fasta', './',
                                        filter="fasta file(*.fa *.fasta);;All Files(*.*)")
        if not filename:
            return
        self.aln = AlignIO.read(filename, format='fasta')
        print ('loaded alignment with %s rows' %len(self.aln))
        return

    def load_folder(self, path=None):
        """Load files from a folder"""

        if path == None:
            options = QFileDialog.Options()
            path = QFileDialog.getExistingDirectory(self,"Select folder",
                                                os.getcwd(),
                                                QFileDialog.ShowDirsOnly)
        if not path:
            return
        meta_file = os.path.join(path, 'metadata.csv')
        snp_file = os.path.join(path, 'snpdist.csv')
        gdf_file = os.path.join(path, 'gdf.shp')
        self.load_data(meta_file, snp_file, gdf_file)
        return

    def load_data(self, meta_file, snp_file, gdf_file=None, moves_file=None):
        """Load datasets"""

        df = pd.read_csv(meta_file)
        index_col = 'id'
        df.set_index(index_col,inplace=True)
        t = self.meta_table
        t.setDataFrame(df)
        t.resizeColumns()

        #try to make GeoDataFrame from input metadata
        if gdf_file == None:
            self.cent = self.gdf_from_table(df)
        else:
            gdf = gpd.read_file(gdf_file)
            diffcols = df.columns.difference(gdf.columns)
            self.cent = gdf.merge(df[diffcols], left_on='id', right_index=True)
            self.cent = self.cent.set_index('id')
            #print(self.cent)

        for col in cladelevels:
            self.cent[col] = self.cent[col].astype(str)
        self.update_clades()
        #snps
        #self.coresnps = pd.read_csv(snp_file, sep=' ')
        #self.snpdist = pd.read_csv(snp_file,index_col=0)
        #movement
        #self.moves= pd.read_csv(moves_file)
        self.update()
        return

    def make_test_data(self):
        """Artificial datasets using btbabm"""

        import btbabm
        options = QFileDialog.Options()
        path = QFileDialog.getExistingDirectory(self,"Select folder",
                                            #os.path.expanduser("~"),
                                            os.getcwd(),
                                            QFileDialog.ShowDirsOnly)
        if not path:
            return
        def func(progress_callback):
            btbabm.simulate_test_data(path, steps=500)
        def completed():
            #process meta data
            m = pd.read_csv(os.path.join(path,'meta.csv'))
            gdf = gpd.read_file(os.path.join(path,'gdf.shp'))
            gdf = gdfrename(columns={'id':'Animal_ID','herd':'HERD_NO'})
            #m.to_csv(os.path.join,'meta.csv')
            #process moves
            #moves =
            #self.load_folder(path)
            self.cent = gdf
            self.meta_table.setDataFrame(self.cent)
            self.processing_completed()

        self.run_threaded_process(func, completed)
        return

    def gdf_from_table(self, df, x='X_COORD',y='Y_COORD'):

        cent = gpd.GeoDataFrame(df,geometry=gpd.points_from_xy(df[x], df[y])).set_crs('EPSG:29902')
        return cent

    def get_tabs(self):

        n=[]
        for i in range(self.tabs.count()):
            n.append(self.tabs.tabText(i))
        return n

    def processing_completed(self):
        """Generic process completed"""

        self.progressbar.setRange(0,1)
        self.running = False
        self.meta_table.refresh()
        print ('finished')
        return

    def run_threaded_process(self, process, on_complete):
        """Execute a function in the background with a worker"""

        if self.running == True:
            return
        worker = Worker(fn=process)
        self.threadpool.start(worker)
        worker.signals.finished.connect(on_complete)
        worker.signals.progress.connect(self.progress_fn)
        self.progressbar.setRange(0,0)
        return

    def progress_fn(self, msg):

        self.info.append(msg)
        self.info.verticalScrollBar().setValue(1)
        return

    def get_counties(self):
        self.counties.groupby('sample').bounds()
        return

    def make_phylogeny(self, infile, treefile='tree.newick'):
        """Phylogeny from sequence alignment file"""

        try:
            utils.run_fasttree(infile, treefile, bootstraps=50)
        except Exception as e:
            print ('fasttree error')
            print(e)
            return
        aln = AlignIO.read(infile,'fasta')
        ls = len(aln[0])
        utils.convert_branch_lengths(treefile, treefile, ls)
        return

    def update(self):
        """Update plot"""

        self.plotview.clear()
        ax = self.plotview.ax
        fig = self.plotview.fig

        ms = self.markersizew.value()
        colorcol = self.colorbyw.currentText()
        cmap = self.cmapw.currentText()
        legend = self.legendb.isChecked()
        jitter = self.jitterb.isChecked()

        if self.sub is None or len(self.sub) == 0:
            return
        
        self.plot_counties()
        #land parcels
        if self.parcelsb.isChecked():
            herds = list(self.sub.HERD_NO)
            p = self.parcels[self.parcels.SPH_HERD_N.isin(herds)]
            self.plot_parcels(p, col=colorcol,cmap=cmap)

        if jitter == True:
            idx = self.sub[self.sub.duplicated('geometry')].index
            self.sub['geometry'] = self.sub.apply(lambda x: jitter_points(x,50) if x.name in idx else x.geometry,1)

        if self.showneighboursb.isChecked() and self.neighbours is not None:
            self.neighbours.plot(color='gray',alpha=0.4,ax=ax)

        plot_single_cluster(self.sub,col=colorcol,ms=ms,cmap=cmap,legend=legend,ax=ax)

        labelcol = self.labelsw.currentText()
        if labelcol != '':
            show_labels(self.sub, labelcol, ax)

        #leg = ax.get_legend()
        #leg.set_bbox_to_anchor((0., 0., 1.2, 0.9))

        #moves
        if self.movesb.isChecked():
            mov = get_moves_bytag(self.sub, self.moves, self.lpis_cent)
            plot_moves(mov, self.lpis_cent, ax=ax)
            self.show_moves_table(mov)
            if self.timelineb.isChecked():
                self.show_moves_timeline(mov)

        #tree
        if self.showtreeb.isChecked():
            self.show_tree()
        #mst
        if self.showmstb.isChecked():
            self.plot_mst()

        if self.title != None:
            fig.suptitle(self.title)
        try:
            fig.tight_layout()
        except:
            pass

        #if refreshing using previous plot limits
        lims = self.plotview.lims
        if self.plotview.lims != None:
            ax.set_xlim(lims[0],lims[1])
            ax.set_ylim(lims[2],lims[3])

        #context map
        cxsource = self.contextw.currentText()
        self.add_context_map(providers[cxsource])

        self.plotview.redraw()

        #update subset table
        self.show_selected_table()

        #update folium map
        #if self.showfoliumb.isChecked():
        #    self.foliumview.plot(self.sub, self.parcels, colorcol=colorcol)
        return

    def show_folium(self):
        """Update folium map"""
        
        colorcol = self.colorbyw.currentText()
        cmap = self.cmapw.currentText()
        if self.movesb.isChecked():
            mov = get_moves_bytag(self.sub, self.moves, self.lpis_cent)
        else:
            mov = None
        self.foliumview.plot(self.sub, self.parcels, moves=mov, 
                             lpis_cent=self.lpis_cent, colorcol=colorcol)
        self.tabs.setCurrentIndex(1)
        return

    def plot_farms(self):
        """Seperate farms view in grid"""

        if self.sub is None or len(self.sub) == 0:
            return
        source = providers[self.contextw.currentText()]
        
        herds = list(self.sub.HERD_NO)
        p = self.parcels[self.parcels.SPH_HERD_N.isin(herds)]

        if not hasattr(self, 'gridview'):
            self.gridview = widgets.CustomPlotViewer(self, controls=False, app=self)
            self.m.addWidget(self.gridview)
            idx = self.tabs.addTab(self.gridview, 'Farms')
            self.tabs.setCurrentIndex(idx)

        def plot_completed():
            self.gridview.redraw()
            self.processing_completed()
        def func(progress_callback):
            axs = plot_grid(self.sub, p, fig=self.gridview.fig, source=source)   
        self.run_threaded_process(func, plot_completed) 
        return
    
    def cluster_report(self):
        """Make pdf report"""


        options = QFileDialog.Options()
        filename, _ = QFileDialog.getSaveFileName(self,"Save Report",
                                                  "","pdf files (*.pdf);;All files (*.*)",
                                                  options=options)
        if not filename:
            return

        clades = [item.text(0) for item in self.cladew.selectedItems()]
        colorcol = self.colorbyw.currentText()
        cmap = self.cmapw.currentText()
        legend = self.legendb.isChecked()

        def func(progress_callback):
            from . import reports
            reports.cluster_report(self.cent, self.parcels, self.lpis_cent, moves=self.moves,
                                    alignment=self.aln,
                                    level='snp3', clades=clades, cmap=cmap,
                                    labelcol='Animal_ID', outfile=filename)

        def completed():
            #self.show_pdf(filename)
            self.processing_completed()

        self.run_threaded_process(func, completed)
        return

    def refresh(self):
        """Update with current zoom"""

        self.plotview.redraw()
        return

    def plot_selected_clade(self):
        """Plot points from cluster menu selection"""

        cent = self.cent
        level = self.cladelevelw.currentText()
        clades = [item.text(0) for item in self.cladew.selectedItems()]

        if len(clades) == 0:
            return
        self.sub = cent[cent[level].isin(clades)].copy()
        cl= ','.join(clades)
        self.title = '%s=%s n=%s' %(level,cl,len(self.sub))
        self.plotview.lims = None
        self.update()
        return

    def plot_counties(self):
        """plot county borders"""

        ax = self.plotview.ax
        if self.colorcountiesb.isChecked():
            cty = 'NAME_TAG'
            clr=None
        else:
            cty = None
            clr='none'
        self.counties.plot(edgecolor='gray',column=cty,color=clr,cmap='tab20',alpha=0.4,ax=ax)
        #labels
        if self.colorcountiesb.isChecked():
            c = self.counties
            c['cent'] = c.geometry.centroid
            for x, y, label in zip(c.cent.x, c.cent.y, c["NAME_TAG"]):
                ax.text(x, y, label, fontsize = 12)
        return

    def plot_county(self):
        """Plot all points in a county"""

        cent = self.cent
        county = self.countyw.currentText()
        self.sub = cent[cent.County==county]
        self.plotview.lims = None
        self.title = county
        self.update()
        return

    def plot_table_selection(self):
        """Plot points from table selection"""

        df = self.meta_table.model.df
        rows = self.meta_table.getSelectedRows()
        idx = df.index[rows]
        mask = self.cent.index.isin(idx)
        self.sub = self.cent[mask]
        #self.sub = self.cent.loc[idx]
        self.title = '(table selection) n=%s' %len(self.sub)
        self.plotview.lims = None
        self.update()
        return

    def plot_in_region(self):
        """Show all points in the visible region of plot"""

        xmin,xmax,ymin,ymax = self.plotview.get_plot_lims()
        df = self.cent
        self.sub = df.cx[xmin:xmax, ymin:ymax]
        self.title = ('selected region n=%s' %len(self.sub))
        self.update()
        #ax = self.plotview.ax
        #ax.set_xlim(xmin,xmax)
        #ax.set_ylim(ymin,ymax)
        return

    def plot_parcels(self, parcels, col, cmap='Set1'):
        """Show selected land parcels"""

        if len(parcels) == 0 or parcels is None:
            return
        ax = self.plotview.ax
        parcels.plot(column='SPH_HERD_N',alpha=0.6,lw=1,cmap=cmap,ax=ax)
        return

    def plot_herd_selection(self, herd_no):
        """Plot farm(s)"""

        self.sub = self.cent[self.cent.HERD_NO.isin(herd_no)]
        self.title = '(herd selection) %s' %' '.join(list(herd_no))
        self.parcelsb.setChecked(True)
        self.update()
        #set bounds to parcels
        self.set_bounds(self.parcels)
        return

    def set_bounds(self, gdf):

        margin=10
        ax=self.plotview.ax
        minx, miny, maxx, maxy = gdf.total_bounds
        ax.set_xlim(minx-margin,maxx+margin)
        ax.set_ylim(miny-margin,maxy+margin)
        self.plotview.redraw()
        return

    def add_context_map(self, source=None):
        """Contextily background map"""

        if source == None:
            return
        ax = self.plotview.ax
        cx.add_basemap(ax, crs=self.cent.crs, #zoom=18,
                attribution=False, source=source)
        return

    def show_moves_timeline(self, df):
        """Show moves timeline plot"""

        fig,ax = plt.subplots(1,1)
        #w = widgets.CustomPlotViewer(self, controls=False, app=self)
        w = widgets.PlotWidget(self)
        plot_moves_timeline(df,w.ax)
        #w.redraw()
        self.show_dock_object(w, 'timeline')
        return

    def show_tree(self):
        """Show phylogeny for selected subset"""

        colorcol = self.colorbyw.currentText()
        cmap = self.cmapw.currentText()
        labelcol = self.labelsw.currentText()

        from Bio.Align import MultipleSeqAlignment
        import toyplot

        idx = list(self.sub.index)
        seqs = [rec for rec in self.aln if rec.id in idx]
        aln = MultipleSeqAlignment(seqs)
        #print (aln)
        treefile = trees.tree_from_aln(aln)
        w = QWebEngineView()
        if labelcol == '':
            labelcol=None
        canvas = trees.draw_tree(treefile, self.sub, colorcol, cmap, tiplabelcol=labelcol)
        toyplot.html.render(canvas, "temp.html")
        with open('temp.html', 'r') as f:
            html = f.read()
            w.setHtml(html)
        self.show_dock_object(w, 'tree')
        return

    def plot_mst(self):
        """Show min spanning tree for selected subset"""

        cmap = self.cmapw.currentText()
        colorcol = self.colorbyw.currentText()

        from Bio.Align import MultipleSeqAlignment
        idx = list(self.sub.index)
        seqs = [rec for rec in self.aln if rec.id in idx]
        aln = MultipleSeqAlignment(seqs)
        D = tools.snp_dist_matrix(aln)
        w = widgets.PlotWidget(self)
        tools.dist_matrix_to_mst(D,self.sub,colorcol, cmap=cmap,ax=w.ax)
        self.show_dock_object(w, 'mst')
        return

    def show_scratchpad(self):
        """Show the scratchpad"""

        if not hasattr(self, 'scratchpad'):
            self.scratchpad = widgets.ScratchPad()
            try:
                self.scratchpad.resize(self.settings.value('scratchpad_size'))
            except:
                pass
        self.scratchpad.update(self.scratch_items)
        self.scratchpad.show()
        self.scratchpad.activateWindow()
        return

    def save_to_scratchpad(self, label=None):
        """Save plot to scratchpad"""

        name = self.title
        if label == None or label is False:
            t = time.strftime("%H:%M:%S")
            label = name+'-'+t
        #get the current figure and make a copy of it by using pickle

        fig = self.plotview.fig
        p = pickle.dumps(fig)
        fig = pickle.loads(p)

        self.scratch_items[label] = fig
        if hasattr(self, 'scratchpad'):
            self.scratchpad.update(self.scratch_items)
        return

    def get_selected(self):
        """Get selected rows of fastq table"""

        df = self.meta_table.model.df
        rows = self.meta_table.getSelectedRows()
        if len(rows) == 0:
            print ('no samples selected')
            return
        data = df.iloc[rows]
        return data

    def herd_summary(self):
        """Summary by herd"""

        res=[]
        for herd, df in self.cent.groupby('HERD_NO'):
            #print (herd)
            clades = len(df.snp3.unique())
            m = self.moves[self.moves.move_to == herd]
            res.append([herd,len(df),clades,len(m)])
        res = pd.DataFrame(res, columns=['HERD_NO','isolates','strains','moves'])
        res = res.sort_values('strains',ascending=False)
        w = tables.HerdTable(self, dataframe=res,
                    font=core.FONT, fontsize=core.FONTSIZE, app=self)
        self.show_dock_object(w, 'herd summary')
        return

    def show_dock_object(self, widget, name, side='right'):
        """Show a table in the dock"""

        if not name in self.docks:
            dock = self.add_dock(widget,name,side)
            self.docks[name] = dock
            self.add_dock_item(name)
        else:
            self.docks[name].setWidget(widget)
        if type(widget) is tables.SampleTable:
            self.opentables[name] = widget
        return

    def show_selected_table(self):
        """Show selected samples in separate table"""

        w = tables.SampleTable(self, dataframe=self.sub,
                font=core.FONT, fontsize=core.FONTSIZE, app=self)
        self.show_dock_object(w, 'selected', 'left')
        return

    def show_moves_table(self, df):
        """Show moves for samples in separate table"""

        if df is None:
            df = pd.DataFrame()
        w = tables.MovesTable(self, dataframe=df,
                            font=core.FONT, fontsize=core.FONTSIZE, app=self)
        self.show_dock_object(w, 'moves')
        return

    def sample_details(self, data):
        """Show sample details"""

        w = tables.DataFrameTable(self, pd.DataFrame(data),
                                font=core.FONT, fontsize=core.FONTSIZE)
        w.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.show_dock_object(w, 'sample details')
        return

    def strain_typing(self):

        from . import strain_typing
        self.st = strain_typing.StrainTypingTool(self)
        self.add_dock(self.st, 'strain typing', 'right')
        self.add_dock_item('strain typing')
        return

    def zoom_in(self):

        core.FONTSIZE+=1
        for i in self.opentables:
            w=self.opentables[i]
            w.zoomIn(core.FONTSIZE)
        self.info.zoomIn()
        return

    def zoom_out(self):

        core.FONTSIZE-=1
        for i in self.opentables:
            w=self.opentables[i]
            w.zoomOut(core.FONTSIZE)
        self.info.zoomOut()
        return

    def show_browser_tab(self, link, name):
        """Show browser"""

        browser = QWebEngineView()
        browser.setUrl(QUrl(link))
        idx = self.tabs.addTab(browser, name)
        self.tabs.setCurrentIndex(idx)
        return

    def show_info(self, msg, color=None):

        if color != None:
            #self.info.appendHtml("<p style=\"color:red\">" + msg + "</p>")
            self.info.append("<font color=%s>%s</font>" %(color,msg))
        else:
            self.info.append(msg)
        self.info.verticalScrollBar().setValue(
            self.info.verticalScrollBar().maximum())
        return

    def get_tab_names(self):
        return {self.tabs.tabText(index):index for index in range(self.tabs.count())}

    def preferences(self):
        """Preferences dialog"""

        opts = {}
        for k in core.defaults.keys():
            opts[k] = getattr(core,k)
        dlg = widgets.PreferencesDialog(self, opts)
        dlg.exec_()
        return

    def quit(self):
        self.close()
        return

    def closeEvent(self, event=None):

        if self.proj_file != None and event != None:
            reply = QMessageBox.question(self, 'Confirm', "Save the current project?",
                                            QMessageBox.Cancel | QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Cancel:
                event.ignore()
                return
            elif reply == QMessageBox.Yes:
                self.save_project()
        self.save_settings()
        event.accept()
        return

    def show_pdf(self, filename):
        """Show a pdf"""

        self.show_browser_tab(filename, filename)
        return

    def online_documentation(self,event=None):
        """Open the online documentation"""

        link='https://github.com/dmnfarrell/tracebtb/wiki'
        self.show_browser_tab(link, 'Help')
        return

    def about(self):

        from . import __version__
        import matplotlib
        try:
            import PySide2
            pyqtver = 'PySide2 v'+PySide2.QtCore.__version__
        except:
            from PyQt5.QtCore import PYQT_VERSION_STR
            pyqtver = 'PyQt5 v'+PYQT_VERSION_STR
        pandasver = pd.__version__
        pythonver = platform.python_version()
        mplver = matplotlib.__version__

        text='TracebTB\n'\
            +'version '+__version__+'\n'\
            +'Copyright (C) Damien Farrell 2022-\n'\
            +'This program is free software; you can redistribute it and/or '\
            +'modify it under the terms of the GNU GPL '\
            +'as published by the Free Software Foundation; either '\
            +'version 3 of the License, or (at your option) any '\
            +'later version.\n'\
            +'Using Python v%s, %s\n' %(pythonver, pyqtver)\
            +'pandas v%s, matplotlib v%s' %(pandasver,mplver)

        msg = QMessageBox.about(self, "About", text)
        return

#https://www.learnpyqt.com/courses/concurrent-execution/multithreading-pyqt-applications-qthreadpool/
class Worker(QtCore.QRunnable):
    """Worker thread for running background tasks."""

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()
        self.kwargs['progress_callback'] = self.signals.progress

    @QtCore.Slot()
    def run(self):
        try:
            result = self.fn(
                *self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)
        finally:
            self.signals.finished.emit()

class WorkerSignals(QtCore.QObject):
    """
    Defines the signals available from a running worker thread.
    Supported signals are:
    finished
        No data
    error
        `tuple` (exctype, value, traceback.format_exc() )
    result
        `object` data returned from processing, anything
    """
    finished = QtCore.Signal()
    error = QtCore.Signal(tuple)
    result = QtCore.Signal(object)
    progress = QtCore.Signal(str)

class StdoutRedirect(QObject):
    printOccur = Signal(str, str, name="print")

    def __init__(self, *param):
        QObject.__init__(self, None)
        self.daemon = True
        self.sysstdout = sys.stdout.write
        self.sysstderr = sys.stderr.write

    def stop(self):
        sys.stdout.write = self.sysstdout
        sys.stderr.write = self.sysstderr

    def start(self):
        sys.stdout.write = self.write
        sys.stderr.write = lambda msg : self.write(msg, color="red")

    def write(self, s, color="black"):
        sys.stdout.flush()
        self.printOccur.emit(s, color)

class AppOptions(widgets.BaseOptions):
    """Class to provide a dialog for global plot options"""

    def __init__(self, parent=None):
        """Setup variables"""
        self.parent = parent
        self.kwds = {}
        cpus = os.cpu_count()
        self.groups = {'general':['threads','labelsep','overwrite'],
                       }
        self.opts = {'threads':{'type':'spinbox','default':4,'range':(1,cpus)},
                    }
        return

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='TracebTB')
    parser.add_argument("-p", "--proj", dest="project",default=None,
                        help="load project file", metavar="FILE")
    args = vars(parser.parse_args())
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_ShareOpenGLContexts)
    app = QApplication(sys.argv)
    aw = App(**args)
    aw.show()
    app.exec_()

if __name__ == '__main__':
    main()
