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

#fix for browser display
os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--no-sandbox"
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
    if len(df) == 1:
        minx-=1000
        miny-=1000
        maxx+=1000
        maxy+=1000
    else:
        maxx+=10
        maxy+=10

    #df['color'] = df.Species.map({'Bovine':'blue','Badger':'orange'})
    #map color to col here properly?
    df['color'] = 'blue'
    ncols = 1
    if col in df.columns and len(df.groupby(col).groups)>15:
        ncols=2
    legfmt = {'title':col,'fontsize':'small','frameon':False,'draggable':True,'ncol':ncols}
    if col == None or col == '':
        df.plot(color=df.color,ax=ax,alpha=0.6,markersize=ms,linewidth=.5,
                legend=legend,legend_kwds=legfmt)
    else:
        df.plot(column=col,ax=ax,alpha=0.6,markersize=ms,linewidth=.5,cmap=cmap,
                legend=legend,legend_kwds=legfmt)
    cow = df[df.Species=='Bovine']
    badger = df[df.Species=='Badger']
    #cow.plot(column=col,ax=ax,alpha=0.6,markersize=ms,linewidth=1,ec='black',cmap=cmap,legend=legend)
    #outline badgers with shape
    if len(badger)>0:
        badger.plot(color='none',ax=ax,alpha=0.8,marker='^',markersize=ms+40,linewidth=.5,ec='black')
    ax.set_title(title)
    ax.axis('off')
    if margin == None:
        margin = (maxx-minx)*0.3
    ax.set_xlim(minx-margin,maxx+margin)
    ax.set_ylim(miny-margin,maxy+margin)
    ax.add_artist(ScaleBar(dx=1,location=3))

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
        #print (axs)
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

def plot_moves(moves, lpis_cent, ax, ms=80):
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
                            markersize=ms,linewidth=.8,alpha=0.5,ax=ax)
                i+=1
    return

def plot_parcels(parcels, ax, col=None, cmap='Set1'):
    """Show selected land parcels"""

    if len(parcels) == 0 or parcels is None:
        return
    if col == '' or col == None:
        parcels.plot(color='none',lw=.5,ec='black',cmap=cmap,ax=ax)
    else:
        parcels.plot(column=col,alpha=0.6,lw=.3,ec='black',cmap=cmap,ax=ax)
    return

def plot_moves_timeline(df, ax=None):
    """Timeline from moves"""

    from datetime import datetime, timedelta
    from matplotlib.patches import Rectangle
    import matplotlib.dates as mdates

    df['move_date'] = pd.to_datetime(df.move_date)
    groups = df.groupby('Animal_ID')
    if (len(groups))>25:
        print ('too many moves to plot, reduce selection')
        return
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

def apply_jitter(gdf, radius=100):
    """Apply jitter to points on same farm"""

    def circular_jitter(centroid, points, radius):
        angles = np.linspace(0, 2 * np.pi, len(points), endpoint=False)
        jittered_points = []
        for angle in angles:
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            jittered_points.append(Point(centroid.x + x, centroid.y + y))
        return jittered_points

    #Group by 'HERD_NO' and apply jitter to overlapping points
    for herd, group in gdf.groupby('HERD_NO'):
        #print (group)
        if len(group) <= 1: continue
        centroid = group['geometry'].iloc[0]
        if centroid == None or centroid.is_empty: continue
        jittered = circular_jitter(centroid, group['geometry'], radius)
        #print (jittered)
        gdf.loc[group.index, 'geometry'] = jittered
    return gdf

def get_moves_bytag(df, move_df, lpis_cent):
    """
    Get moves and coords for one or more samples.
    """

    cols=['Animal_ID']+list(move_df.columns)
    t = df.merge(move_df,left_on='Animal_ID',right_on='tag',how='inner')[cols]
    #print (t)
    #merge result with parcels to get coords of moved_to farms
    m = t.merge(lpis_cent,left_on='move_to',right_on='SPH_HERD_N', how='left')

    if len(m)==0:
        return
    #add in source farms - is this needed?
    #x = lpis_cent[lpis_cent.SPH_HERD_N.isin(df.HERD_NO)]
    #m = pd.concat([m,x])#.dropna(subset='Animal_ID')
    #print (m.iloc[0])
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

        self.load_base_data()
        self.recent_files = ['']
        self.scratch_items = {}
        self.selections = {}
        self.opentables = {}

        self.lpis_master = None
        self.lpis_master_file = None
        self.lpis_cent = None
        self.parcels = None
        self.neighbours = None
        self.sub = None
        #self.moves = None

        self.main.setFocus()
        self.setCentralWidget(self.main)
        self.create_tool_bar()
        self.load_settings()
        self.setup_gui()
        self.show_recent_files()
        self.update_selections_menu()

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
                core.THREADS = int(s.value('threads'))
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
        self.settings.setValue('threads', core.THREADS)
        self.settings.setValue('recent_files',','.join(self.recent_files))
        #print (self.settings)
        self.settings.sync()
        return

    def load_base_data(self):
        """Set reference map of counties"""

        self.counties = counties_gdf
        return

    def create_tool_bar(self):
        """Create main toolbar"""

        items = {'New project': {'action': lambda: self.new_project(ask=True),'file':'document-new'},
                 'Open': {'action':self.load_project,'file':'document-open'},
                 'Save': {'action': lambda: self.save_project(),'file':'save'},
                 'Zoom out': {'action':self.zoom_out,'file':'zoom-out'},
                 'Zoom in': {'action':self.zoom_in,'file':'zoom-in'},
                 'Scratchpad': {'action':self.show_scratchpad,'file':'scratchpad'},
                 'Filter': {'action':self.show_filter,'file':'filter'},
                 'Settings': {'action':self.preferences,'file':'settings'},
                 'Build Tree': {'action':self.show_tree,'file':'tree'},
                 'Show MST': {'action':self.plot_mst,'file':'mst'},
                 'Herd Summary':{'action':self.herd_summary,'file':'cow'},
                 'Cluster Report':{'action':self.cluster_report ,'file':'cluster_report'},
                 'Make Simulated Data':{'action':self.simulate_data ,'file':'simulate'},
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

        df = self.meta_table.model.df
        self.widgets = {}
        m = QWidget()
        m.setStyleSheet(style)
        m.setMaximumWidth(140)
        m.setMinimumWidth(140)
        l = QVBoxLayout()
        l.setAlignment(QtCore.Qt.AlignTop)
        m.setLayout(l)

        #select grouping column
        self.groupbyw = w = QComboBox(m)
        #w.addItems(cols)
        l.addWidget(QLabel('Group by:'))
        l.addWidget(w)
        w.currentIndexChanged.connect(self.update_groups)
        self.widgets['cladelevel'] = w
        #select clade/cluster
        l.addWidget(QLabel('Group:'))
        t = self.groupw = QTreeWidget()
        t.setHeaderItem(CustomTreeWidgetItem(["name","size"]))
        t.setColumnWidth(0, 60)
        t.setColumnWidth(1, 20)
        t.setSortingEnabled(True)
        t.setMinimumHeight(30)
        #t.setSelectionMode(QAbstractItemView.ExtendedSelection)
        #t.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        #t.customContextMenuRequested.connect(self.update_groups)
        t.itemSelectionChanged.connect(self.plot_selected_group)
        l.addWidget(t)
        self.widgets['clade'] = t

        #zoom to county
        #self.countyw =w = QComboBox(m)
        #w.addItems(['']+counties)
        #l.addWidget(QLabel('County:'))
        #l.addWidget(w)
        #w.currentIndexChanged.connect(self.plot_county)
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
        #color parcels by
        self.colorparcelsbyw = w = QComboBox(m)
        l.addWidget(QLabel('Color parcels by:'))
        l.addWidget(w)
        w.setMaxVisibleItems(12)
        self.widgets['colorparcelsby'] = w
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
        m.setMaximumWidth(65)
        m.setMinimumWidth(65)
        l = QVBoxLayout()
        l.setAlignment(QtCore.Qt.AlignTop)
        m.setLayout(l)
        b = widgets.createButton(m, None, self.update, 'refresh', core.ICONSIZE)
        l.addWidget(b)
        self.showfoliumb = b = widgets.createButton(m, None, self.update, 'folium', core.ICONSIZE, 'interactive view')
        b.setCheckable(True)
        l.addWidget(b)
        #self.multiaxesb= b = widgets.createButton(m, None, self.plot_farms, 'plot-grid', core.ICONSIZE, 'farm split view')
        #l.addWidget(b)
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
        self.widgets['colorcounties'] = b
        #self.jitterb = b = widgets.createButton(m, None, self.update, 'jitter', core.ICONSIZE, 'jitter points')
        #b.setCheckable(True)
        #l.addWidget(b)
        self.showneighboursb = b = widgets.createButton(m, None, self.get_neighbouring_parcels,
                                                        'neighbours', core.ICONSIZE, 'show neighbours')
        #b.setCheckable(True)
        l.addWidget(b)
        #self.widgets['showneighbours'] = b
        b = widgets.createButton(m, None, self.get_neighbours_in_region, 'parcels-region',
                                 core.ICONSIZE, 'show all parcels in region')
        l.addWidget(b)
        #self.showtreeb = b = widgets.createButton(m, None, self.update, 'tree', core.ICONSIZE, 'show tree')
        #b.setCheckable(True)
        #l.addWidget(b)
        #self.widgets['showtree'] = b
        #self.showmstb = b = widgets.createButton(m, None, self.update, 'mst', core.ICONSIZE, 'show MST')
        #b.setCheckable(True)
        #l.addWidget(b)
        b = widgets.createButton(m, None, self.save_to_scratchpad, 'snapshot', core.ICONSIZE, 'take snapshot')
        l.addWidget(b)
        return m

    def update_groups(self):
        """Update the 'groups' widget. Called when grouping column is changed."""

        df = self.meta_table.model.df
        groupby = self.groupbyw.currentText()
        if groupby not in df.columns:
            return
        vals = df[groupby].value_counts()
        #print (vals)
        t = self.groupw
        t.clear()
        for g,size in vals.items():
            item = CustomTreeWidgetItem(t)
            item.setTextAlignment(0, QtCore.Qt.AlignLeft)
            item.setText(0, str(g))
            item.setText(1, str(size))
        t.sortItems(0, QtCore.Qt.SortOrder(0))
        return

    def get_ordinal_columns(self):
        """Try to get ordinal columns from table"""

        ignore = ['sample','Aliquot','geometry','Animal_ID','X_COORD','Y_COORD']
        df = self.meta_table.model.df
        ocols = []
        for col in df.columns:
            count = df[col].nunique()
            #print (col, count, len(df[col]))
            if count==len(df[col]) or col in ignore:
                continue
            ocols.append(col)
        return ocols

    def update_widgets(self):
        """Update widgets when new table loaded"""

        df = self.meta_table.model.df
        cols = ['']+list(df.columns)
        ocols = self.get_ordinal_columns()

        self.groupbyw.clear()
        self.groupbyw.addItems(ocols)
        self.labelsw.clear()
        self.labelsw.addItems(cols)
        self.colorbyw.clear()
        self.colorbyw.addItems(cols)
        self.colorbyw.setCurrentText('')
        self.colorparcelsbyw.clear()
        if self.parcels is not None:
            cols = ['']+list(self.parcels.columns)
            self.colorparcelsbyw.addItems(cols)
        self.colorparcelsbyw.setCurrentText('')
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
        self.tabs.setTabsClosable(False)
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

        w = self.data_statusbar()
        self.statusBar.addWidget(w)

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

    def data_statusbar(self):
        """Add a widget to show if specific data is present"""

        w = QWidget()
        layout = QHBoxLayout(w)
        values = {'lpis_master':'parcels-master','snpdist':'snp-dist','aln':'alignment'}
        S = self.status_icons = {}

        for name, filename in values.items():
            icon = QLabel()
            #set the icon based on the status
            iconfile = os.path.join(iconpath,filename)
            pixmap = QPixmap(iconfile).scaled(QtCore.QSize(20,20),
                                              Qt.KeepAspectRatio, Qt.SmoothTransformation)
            icon.setPixmap(pixmap)
            icon.setAlignment(Qt.AlignCenter)
            icon.setEnabled(False)
            layout.addWidget(icon)
            S[name] = icon
        return w

    def update_data_status(self):
        """Update icons that show staus of loaded data"""

        S = self.status_icons
        for name in S:
            icon = S[name]
            if not hasattr(self, name) or self.__dict__[name] is None:
                icon.setEnabled(False)
            else:
                icon.setEnabled(True)
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
        icon = QIcon(os.path.join(iconpath,'document-new.svg'))
        self.file_menu.addAction('New Project', lambda: self.new_project(ask=True))
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

        self.data_menu = QMenu('Data', self)
        self.menuBar().addMenu(self.data_menu)
        self.data_menu.addAction('Load Samples', lambda: self.load_samples())
        icon = QIcon(os.path.join(iconpath,'shapefile.svg'))
        self.data_menu.addAction(icon,'Load Parcels', self.load_parcels)
        icon = QIcon(os.path.join(iconpath,'plot-moves.svg'))
        self.data_menu.addAction(icon, 'Load Moves', lambda: self.load_moves())
        self.data_menu.addAction('Load Alignment', lambda: self.load_alignment())
        icon = QIcon(os.path.join(iconpath,'snp-dist.svg'))
        self.data_menu.addAction(icon, 'Load SNP Distance Matrix', lambda: self.load_snp_dist())
        self.data_menu.addAction('Load Simulated Data', lambda: self.load_folder())
        self.data_menu.addSeparator()
        icon = QIcon(os.path.join(iconpath,'parcels-master.svg'))
        self.data_menu.addAction(icon, 'Set Master Parcels file', self.set_lpis_file)
        icon = QIcon(os.path.join(iconpath,'shapefile.svg'))
        self.data_menu.addAction(icon, 'Load Master Parcels', self.load_lpis_master)
        icon = QIcon(os.path.join(iconpath,'parcels.svg'))
        self.data_menu.addAction(icon, 'Extract Parcels/Centroids', self.get_lpis_centroids)
        #icon = QIcon(os.path.join(iconpath,'neighbours.svg'))
        #self.data_menu.addAction(icon,'Extract Neighbouring Parcels', self.get_neighbouring_parcels)
        icon = QIcon(os.path.join(iconpath,'clusters.svg'))
        self.data_menu.addAction(icon, 'Get Clusters from SNPs', self.get_clusters)
        icon = QIcon(os.path.join(iconpath,'cow.svg'))
        self.data_menu.addAction(icon, 'Count Animal Moves', self.count_animal_moves)
        self.data_menu.addSeparator()

        self.tools_menu = QMenu('Tools', self)
        self.menuBar().addMenu(self.tools_menu)
        icon = QIcon(os.path.join(iconpath,'tree.svg'))
        self.tools_menu.addAction(icon, 'Build Tree', self.show_tree)
        icon = QIcon(os.path.join(iconpath,'mst.svg'))
        self.tools_menu.addAction(icon, 'Show MST', self.plot_mst)
        icon = QIcon(os.path.join(iconpath,'mbovis.svg'))
        self.tools_menu.addAction(icon, 'Strain Typing', self.strain_typing)
        icon = QIcon(os.path.join(iconpath,'cow.svg'))
        self.tools_menu.addAction(icon, 'Show Herd Summary', self.herd_summary)
        icon = QIcon(os.path.join(iconpath,'simulate.svg'))
        self.tools_menu.addAction(icon, 'Make Simulated Data', self.simulate_data)
        icon = QIcon(os.path.join(iconpath,'pdf.svg'))
        self.tools_menu.addAction(icon,'Cluster Report', self.cluster_report)

        self.scratch_menu = QMenu('Scratchpad', self)
        self.menuBar().addMenu(self.scratch_menu)
        icon = QIcon(os.path.join(iconpath,'scratchpad.svg'))
        self.scratch_menu.addAction(icon,'Show Scratchpad', lambda: self.show_scratchpad())
        icon = QIcon(os.path.join(iconpath,'snapshot.svg'))
        self.scratch_menu.addAction(icon,'Plot to Scratchpad', lambda: self.save_to_scratchpad())

        self.selections_menu = QMenu('Selections', self)
        self.menuBar().addMenu(self.selections_menu)
        self.selections_menu.addAction('Save Selection', lambda: self.save_selection())
        self.selections_menu.addAction('Clear Selections', lambda: self.clear_selections(ask=True))
        self.selections_menu.addSeparator()

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
        data['meta'] = self.meta_table.model.df
        keys = ['sub','moves','parcels','lpis_cent','aln','snpdist','selections','lpis_master_file']
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
        #self.meta_table.setDataFrame(pd.DataFrame({'sample':[]}))

        self.plotview.clear()
        self.foliumview.clear()
        self.groupw.clear()
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
        keys = ['sub','moves','parcels','lpis_cent','aln','snpdist','selections','lpis_master_file']
        for k in keys:
            if k in data:
                self.__dict__[k] = data[k]

        #load main table
        t = self.meta_table
        t.setDataFrame(data['meta'])
        self.table_widget.updateStatusBar()
        self.update_groups()
        self.update_widgets()

        if 'widget_values' in data:
            #print (data['widget_values'])
            widgets.setWidgetValues(self.widgets, data['widget_values'])
        self.proj_file = filename
        self.projectlabel.setText(self.proj_file)
        if 'scratch_items' in data:
            self.scratch_items = data['scratch_items']
        self.add_recent_file(filename)
        self.update_selections_menu()
        self.update()
        if 'dock_items' in data:
            self.update_dock_items(data['dock_items'])
        self.update_data_status()
        #t = self.groupw
        #t.setCurrentIndex(t.model().index(0, 4))
        #print (t.selectedIndexes())
        #self.parcels['HERD_NO'] = self.parcels.SPH_HERD_N
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

    def load_test(self):
        """Load test dataset"""

        #metadata
        df = pd.read_csv('testing/metadata.csv')
        #index_col = 'sample'
        #df.set_index(index_col,inplace=True)

        for col in cladelevels:
            df[col] = df[col].astype(str)
        t = self.meta_table
        t.setDataFrame(df)
        t.resizeColumns()
        self.update_groups()
        self.update_widgets()
        self.plot_selected_group()

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
        print ('setting LPIS file to %s' %filename)
        self.lpis_master_file = filename
        return

    def load_lpis_master(self):
        """Load the master parcels file"""

        if self.lpis_master_file == None:
            print ('no master file set')
            return
        print ('reading master shapefile..')
        def completed():
            self.processing_completed()
            self.update_data_status()

        def func(progress_callback):
            self.lpis_master = gpd.read_file(self.lpis_master_file).set_crs('EPSG:29902')
        self.run_threaded_process(func, completed)
        return

    def get_lpis_centroids(self):
        """Calculate LPIS centroids and parcels for samples, run once when the metadata
           is updated using land parcel data."""

        if self.lpis_master is None:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("No LPIS master file loaded")
            msg.exec_()
            return

        def completed():
            #add locations to meta data
            print ('updating table..')
            df = self.meta_table.model.df
            #herds to find
            herds = list(df.HERD_NO)
            #also get intermediate herds from moves if present
            if self.moves is not None:
                herds.extend(self.moves.move_to)
            self.set_locations()
            x = self.lpis_master
            self.parcels = x[x.SPH_HERD_N.isin(herds)]
            print ('found %s parcels' %len(self.parcels))
            self.processing_completed()
            return

        def func(progress_callback):
            self.lpis_cent = calculate_parcel_centroids(self.lpis_master)
        print ('calculating centroids..')
        self.run_threaded_process(func, completed)
        return

    def get_neighbouring_parcels(self):
        """
        Find all neighbouring parcels from LPIS.
        Requires the LPIS master file.
        """

        df = self.sub
        if self.lpis_master is None:
            print ('LPIS master file is not loaded.')
            return

        lpis = self.lpis_master
        #add farms that are in current moves data aswell
        mov = get_moves_bytag(self.sub, self.moves, self.lpis_cent)
        m = lpis[lpis.SPH_HERD_N.isin(mov.SPH_HERD_N)]
        #combine both
        df = pd.concat([df,m])

        def completed():
            self.processing_completed()
            self.update()
            print ('found %s parcels' %len(self.neighbours))
            return

        def func(progress_callback):
            d=1200
            found = []
            for x in df.geometry:
                dists = self.lpis_cent.distance(x)
                points = self.lpis_cent[(dists<=d) & (dists>10)]
                found.append(points)

            found = pd.concat(found).drop_duplicates()
            p = lpis[lpis.SPH_HERD_N.isin(found.SPH_HERD_N)]
            p['color'] = p.apply(tools.random_grayscale_color, 1)
            self.neighbours = p

        print ('getting parcels..')
        self.run_threaded_process(func, completed)
        return

    def get_neighbours_in_region(self):
        """Get all parcels in region from lpis"""

        if self.lpis_master is None:
            print ('LPIS master file is not loaded.')
            return

        xmin,xmax,ymin,ymax = self.plotview.get_plot_lims()
        lpis = self.lpis_master
        p = lpis.cx[xmin:xmax, ymin:ymax]
        p = p[~p.SPH_HERD_N.isin(self.sub.HERD_NO)]
        p['color'] = p.apply(tools.random_grayscale_color, 1)
        #p['color'] = tools.random_colormap_colors('GnBu',len(p))
        self.neighbours = p
        self.update()
        return

    def set_locations(self):
        """Add locations to meta data from lpis_cent. Overwrites any geometry present."""

        if self.lpis_cent is None:
            print ('no centroids')
            return
        df = self.meta_table.model.df
        #remove previous geometry? - ask user
        if type(df) is gpd.GeoDataFrame:
            reply = QMessageBox.question(self, 'Replace?', "Locations already present. Replace?",
                                            QMessageBox.No | QMessageBox.Yes, QMessageBox.No)
            if reply == QMessageBox.No:
                return
            df = pd.DataFrame(df)
            df.drop(columns=['geometry'], inplace=True)

        df = df.merge(self.lpis_cent,left_on='HERD_NO',right_on='SPH_HERD_N',how='left')#.set_index(idx)
        df = gpd.GeoDataFrame(df,geometry=df.geometry).set_crs('EPSG:29902')
        df['geometry'] = [Point() if x is None else x for x in df.geometry]

        #jitter any points in same farm
        print ('jittering points')
        df = apply_jitter(df, radius=100)
        self.meta_table.setDataFrame(df)
        return

    def get_clusters(self):
        """Get clades/clusters for samples"""

        from . import clustering
        def func(progress_callback):
            if not hasattr(self, 'snpdist'):
                if self.aln is not None:
                    print ('calculating snp matrix. may take some time..')
                    self.snpdist = tools.snp_dist_matrix(self.aln)
                else:
                    print ('no alignment or dist matrix found')
                    return
            clusts,members = clustering.get_cluster_levels(self.snpdist)

        def completed():
            print ()
            self.processing_completed()
        self.run_threaded_process(func, completed)
        return

    def count_animal_moves(self):
        """Count animal moves from movement data"""

        if self.moves is None:
            return

        gdf = self.meta_table.model.df
        cols=['Animal_ID']+list(self.moves.columns)
        df = gdf.merge(self.moves,left_on='Animal_ID',right_on='tag',how='inner')[cols]
        df = df[df.data_type=='F_to_F']
        g = df.Animal_ID.value_counts()
        g = df.groupby('Animal_ID').count()['id'].reset_index()
        g = g.rename(columns={'id':'moves'})
        print (g)
        #put into main table

        gdf = gdf.merge(g, on='Animal_ID', how='left')
        gdf['moves'] = gdf.moves.fillna(0).astype(str)
        self.meta_table.setDataFrame(gdf)
        self.update_widgets()
        return

    def load_samples(self, filename=None, index=None):
        """Load meta data"""

        if filename == None:
            reply = QMessageBox.question(self, 'Continue?', "This will overwrite any current meta data. Are you sure?",
                                            QMessageBox.Cancel | QMessageBox.Yes)
            if reply == QMessageBox.Cancel:
                return
            filename, _ = QFileDialog.getOpenFileName(self, 'Open Meta Data', './',
                                            filter="csv file(*.csv *.txt);;All Files(*.*)")
            if not filename:
                return
        df = pd.read_csv(filename)

        for col in cladelevels:
            if col in df.columns:
                df[col] = df[col].astype(str)

        #try to convert to geodataframe if has coords
        result = self.extract_coords(df)
        if result is None:
            print ('no coords found. you can still use parcels to determine locations')
        else:
            df = result
        #set the index column, should be animal ID
        if index == None:
            cols = df.columns
            item, ok = QInputDialog.getItem(self, 'Select Index field', 'Index field:', cols, 0, False)
            df = df.set_index(item, drop=False)
            df.index.name = 'index'

        t = self.meta_table
        t.setDataFrame(df)
        self.update_groups()
        self.update_widgets()
        return

    def extract_coords(self, df):
        """Try to convert coords from X-Y columns"""

        opts = {'X':{'type':'combobox','default':'X','items':['X','X_COORD']},
                'Y':{'type':'combobox','default':'Y','items':['Y','Y_COORD']}
                }
        dlg = widgets.MultipleInputDialog(self, opts, title='Select X-Y columns',
                            width=250,height=150)
        dlg.exec_()
        if not dlg.accepted:
            return
        kwds = dlg.values
        x=kwds['X']
        y=kwds['Y']
        gdf = gpd.GeoDataFrame(df,geometry=gpd.points_from_xy(df[x], df[y])).set_crs('EPSG:29902')
        #jitter the points
        print ('jittering points')
        gdf = apply_jitter(gdf, radius=100)
        return gdf

    def load_moves(self, filename=None):
        """Load movement data"""

        if filename == None:
            filename, _ = QFileDialog.getOpenFileName(self, 'Open CSV file', './',
                                        filter="csv file(*.csv *.txt);;All Files(*.*)")
            if not filename:
                return
        df = pd.read_csv(filename)
        cols = df.columns
        item, ok = QInputDialog.getItem(self, 'Select tag field', 'tag field:', cols, 0, False)
        df = df.rename(columns={item: 'tag'})
        self.moves = df
        print ('loaded %s rows' %len(self.moves))
        return

    def load_alignment(self, filename=None):
        """Load sequence alignment"""

        if filename == None:
            filename, _ = QFileDialog.getOpenFileName(self, 'Open Fasta', './',
                                        filter="fasta file(*.fa *.fasta);;All Files(*.*)")
            if not filename:
                return
        self.aln = AlignIO.read(filename, format='fasta')
        print ('loaded alignment with %s rows' %len(self.aln))
        reply = QMessageBox.question(self, 'Calculate DM?',
                                "Calculate distance matrix?",
                                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            self.dist_matrix_from_alignment()
        self.update_data_status()
        return

    def load_snp_dist(self, filename=None):
        """Load snp dist matrix"""

        if filename == None:
            filename, _ = QFileDialog.getOpenFileName(self, 'Open Matrix', './',
                                        filter="csv file(*.csv *.txt);;All Files(*.*)")
            if not filename:
                return
        self.snpdist = pd.read_csv(filename,index_col=0)
        print ('loaded dist matrix with %s rows' %len(self.snpdist))
        self.update_data_status()
        return

    def load_parcels(self, filename=None, herdcol=None, crs=None):
        """Import a parcels shapefile"""

        if filename == None:
            filename, _ = QFileDialog.getOpenFileName(self, 'Open Shapefile', './',
                                        filter="shapefile(*.shp);;All Files(*.*)")
        df = gpd.read_file(filename)
        #get herd_no column
        cols = df.columns
        if herdcol == None:
            herdcol, ok = QInputDialog.getItem(self, 'Select HERD_NO field', 'HERD_NO field:', cols, 0, False)
        df = df.rename(columns={herdcol: 'SPH_HERD_N'})
        if df.crs is None and crs == None:
            codes = ['EPSG:29902','EPSG:29901','EPSG:27700','EPSG:4362','WGS84']
            crs, ok = QInputDialog.getItem(self, 'No CRS present', 'Select CRS:', codes, 0, False)

        df.set_crs(crs)
        #convert crs to default
        df = df.set_crs("EPSG:29902")
        self.parcels = df
        print ('parcels loaded')
        #if no lpis_cent use parcel centroids - mainly to avoid errors
        if self.lpis_cent is None:
            self.lpis_cent = calculate_parcel_centroids(self.parcels)
        return

    def load_folder(self):
        """Load files from a folder"""

        options = QFileDialog.Options()
        path = QFileDialog.getExistingDirectory(self,"Select folder",
                                            os.getcwd(),
                                            QFileDialog.ShowDirsOnly)
        if not path:
            return
        self.load_sim_data(path)
        return

    def load_sim_data(self, path):
        """Load a simulated dataset from a folder"""

        meta_file = os.path.join(path, 'meta.csv')
        snp_file = os.path.join(path, 'snpdist.csv')
        gdf_file = os.path.join(path, 'parcels.shp')
        moves_file = os.path.join(path, 'movement.csv')
        #samples
        self.load_samples(meta_file, index='id')
        #parcels
        self.load_parcels(gdf_file, herdcol='herd', crs='EPSG:29902')
        #snp dist
        self.snpdist = pd.read_csv(snp_file,index_col=0)
        #movement
        self.load_moves(moves_file)
        self.update()
        self.update_widgets()
        return

    def simulate_data(self):
        """Artificial datasets using btbabm"""

        from . import simulate
        if not hasattr(self, 'simapp'):
            self.simapp = simulate.SimulateApp()
        self.simapp.show()

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
        worker = widgets.Worker(fn=process)
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
        colorparcelscol = self.colorparcelsbyw.currentText()
        cmap = self.cmapw.currentText()
        legend = self.legendb.isChecked()
        #jitter = self.jitterb.isChecked()

        if self.sub is None or len(self.sub) == 0:
            return
        #get moves here
        mov = get_moves_bytag(self.sub, self.moves, self.lpis_cent)

        self.plot_counties()
        #land parcels
        if self.parcelsb.isChecked() and self.parcels is not None:
            herds = list(self.sub.HERD_NO)
            #add parcels for intermediate herds if we have moves
            if self.moves is not None:
                herds.extend(mov.move_to)
            p = self.parcels[self.parcels.SPH_HERD_N.isin(herds)]
            plot_parcels(p, col=colorparcelscol, cmap=cmap, ax=ax)

        #to be removed
        #if jitter == True:
        #    self.sub = apply_jitter(self.sub, radius=100)

        #if self.showneighboursb.isChecked():
        #    self.get_neighbouring_parcels()
        if self.neighbours is not None and self.parcelsb.isChecked():
            #self.neighbours.plot(color='gray',alpha=0.4,ax=ax)
            self.neighbours.plot(column='color',cmap='gray',alpha=0.4,ax=ax)

        plot_single_cluster(self.sub,col=colorcol,ms=ms,cmap=cmap,legend=legend,ax=ax)

        labelcol = self.labelsw.currentText()
        if labelcol != '':
            show_labels(self.sub, labelcol, ax)

        leg = ax.get_legend()
        if leg != None:
            leg.set_bbox_to_anchor((0., 0., 1.2, 0.9))

        #moves
        if self.movesb.isChecked():
            if self.moves is None:
                print ('no moves loaded')
                return

            plot_moves(mov, self.lpis_cent, ax=ax)
            self.show_moves_table(mov)
            if self.timelineb.isChecked():
                self.show_moves_timeline(mov)

        #tree
        #if self.showtreeb.isChecked():
        #    self.show_tree()
        #else:
            #clear tree
            #self.clear_tree()
        #mst
        #if self.showmstb.isChecked():
        #    self.plot_mst()

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

        self.set_equal_aspect()

        #fig.tight_layout()
        #context map
        cxsource = self.contextw.currentText()
        self.add_context_map(providers[cxsource])

        self.plotview.redraw()

        #update subset table
        self.show_selected_table()

        #update folium map
        if self.showfoliumb.isChecked():
            #self.foliumview.plot(self.sub, self.parcels, colorcol=colorcol)
            self.show_folium()
        return

    def show_folium(self):
        """Update folium map"""

        colorcol = self.colorbyw.currentText()
        parcelscol = self.colorparcelsbyw.currentText()
        cmap = self.cmapw.currentText()
        if self.movesb.isChecked():
            mov = get_moves_bytag(self.sub, self.moves, self.lpis_cent)
        else:
            mov = None
        herds = list(self.sub.HERD_NO)
        p = None
        np = None
        if 'HERD_NO' in self.parcels.columns:
            key = 'HERD_NO'
        else:
            key = 'SPH_HERD_N'
        if self.parcelsb.isChecked():
            if self.parcels is not None:
                if mov is not None:
                    herds.extend(mov.move_to)
                p = self.parcels[self.parcels[key].isin(herds)]
            else:
                print ('no parcels found')
            if self.neighbours is not None:
                np = self.neighbours
        self.foliumview.plot(sub=self.sub, parcels=p, neighbours=np, moves=mov,
                             lpis_cent=self.lpis_cent,
                             colorcol=colorcol, parcelscol=parcelscol,
                             cmap=cmap)
        return

    def plot_farms(self):
        """Seperate farms view in grid"""

        if self.sub is None or len(self.sub) == 0:
            return
        source = providers[self.contextw.currentText()]

        herds = list(self.sub.HERD_NO)
        p = self.parcels[self.parcels.HERD_NO.isin(herds)]

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

        df = self.meta_table.model.df
        options = QFileDialog.Options()
        filename, _ = QFileDialog.getSaveFileName(self,"Save Report",
                                                  "","pdf files (*.pdf);;All files (*.*)",
                                                  options=options)
        if not filename:
            return

        clades = [item.text(0) for item in self.groupw.selectedItems()]
        colorcol = self.colorbyw.currentText()
        cmap = self.cmapw.currentText()
        legend = self.legendb.isChecked()

        def func(progress_callback):
            from . import reports
            reports.cluster_report(df, self.parcels, self.lpis_cent, moves=self.moves,
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

    def plot_selected_group(self):
        """Plot points from selected groupby menu selection"""

        gdf = self.meta_table.model.df
        level = self.groupbyw.currentText()
        groups = [item.text(0) for item in self.groupw.selectedItems()]

        if len(groups) == 0:
            return
        self.sub = gdf[gdf[level].isin(groups)].copy()
        cl = ','.join(groups)
        self.title = '%s=%s n=%s' %(level,cl,len(self.sub))
        self.plotview.lims = None
        #reset neighbours
        self.neighbours = None
        self.update()
        return

    def plot_counties(self):
        """plot county borders"""

        ax = self.plotview.ax
        if self.colorcountiesb.isChecked():
            cty = 'NAME_TAG'
            clr = None
        else:
            #ax.set_facecolor('lightblue')
            cty = None
            clr='none'
        self.counties.plot(edgecolor='gray',column=cty,color=clr,
                           cmap='tab20',lw=0.6,alpha=0.7,
                           ax=ax)

        #labels
        if self.colorcountiesb.isChecked():
            #from adjustText import adjust_text
            c = self.counties
            #texts=[]
            c['cent'] = c.geometry.centroid
            for x, y, label in zip(c.cent.x, c.cent.y, c["NAME_TAG"]):
                ax.text(x, y, label, fontsize = 12)
                #texts.append(ax.text(x, y, label, fontsize = 12))
            #adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='->', color='red'))
        return

    def selection_from_table(self):
        """Plot points from table selection"""

        df = self.meta_table.model.df
        idx = self.meta_table.getSelectedIndexes()
        self.sub = df.loc[idx]
        #print (self.sub.iloc[:,:2])
        self.title = '(table selection) n=%s' %len(self.sub)
        self.plotview.lims = None
        self.update()
        return

    def add_to_selection(self):
        """Add selected rows to current selection (self.sub)"""

        df = self.meta_table.model.df
        idx = self.meta_table.getSelectedIndexes()
        new = df.loc[idx]
        self.sub = pd.concat([self.sub,new])
        self.update()
        return

    def select_related(self):
        """Find related samples to selected e.g. within n snps"""

        df = self.meta_table.model.df
        idx = self.meta_table.getSelectedIndexes()
        s = df.loc[idx]

        cols = ['snp3','snp7','snp12']
        col, ok = QInputDialog.getItem(self, 'Select cluster level', 'Cluster level:', cols, 0, False)
        if not col: return
        clusters = s[col]
        self.sub = df[df[col].isin(clusters)]
        self.title = '(table selection) n=%s' %len(self.sub)
        self.plotview.lims = None
        self.update()
        return

    def plot_in_region(self):
        """Show all points in the visible region of plot"""

        xmin,xmax,ymin,ymax = self.plotview.get_plot_lims()
        df = self.meta_table.model.df
        self.sub = df.cx[xmin:xmax, ymin:ymax]
        self.title = ('selected region n=%s' %len(self.sub))
        self.update()
        #ax = self.plotview.ax
        #ax.set_xlim(xmin,xmax)
        #ax.set_ylim(ymin,ymax)
        return

    def plot_herd_selection(self, herd_no):
        """Plot farm(s)"""

        df = self.meta_table.model.df
        self.sub = df[df.HERD_NO.isin(herd_no)]
        self.title = '(herd selection) %s' %' '.join(list(herd_no))
        self.parcelsb.setChecked(True)
        self.update()
        #set bounds to parcels
        self.set_bounds(self.parcels)
        return

    def set_bounds(self, gdf):
        """Set bounds of plot using geodataframe"""

        margin=10
        ax = self.plotview.ax
        minx, miny, maxx, maxy = gdf.total_bounds
        ax.set_xlim(minx-margin,maxx+margin)
        ax.set_ylim(miny-margin,maxy+margin)
        self.plotview.redraw()
        return

    def set_equal_aspect(self):
        """Set axes to be equal regardless of selection"""

        ax = self.plotview.ax
        x1,x2 = ax.get_xlim()
        y1,y2 = ax.get_ylim()
        w=x2-x1
        h=y2-y1
        aspect = w/h
        #print (x1,x2,aspect)
        if aspect >1.1:
            c = y1+h/2
            ax.set_ylim(c-w/2,c+w/2)
        elif aspect<0.9:
            c = x1+w/2
            ax.set_xlim(c-h/2,c+h/2)
        return

    def add_context_map(self, source=None):
        """Contextily background map"""

        gdf = self.meta_table.model.df
        if source == None:
            return
        ax = self.plotview.ax
        cx.add_basemap(ax, crs=gdf.crs, zoom=9,
                attribution=False, source=source)
        return

    def show_moves_timeline(self, df):
        """Show moves timeline plot"""

        fig,ax = plt.subplots(1,1)
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
        idx = list(self.sub.index)
        #print (idx)

        if hasattr(self, 'snpdist'):
            treefile = 'tree.newick'
            M = self.snpdist.loc[idx,idx]
            #print (M)
            trees.tree_from_distmatrix(M, treefile)
        elif hasattr(self, 'aln'):
            from Bio.Align import MultipleSeqAlignment
            seqs = [rec for rec in self.aln if rec.id in idx]
            aln = MultipleSeqAlignment(seqs)
            if len(aln) == 0:
                return
            treefile = trees.tree_from_aln(aln)
        else:
            print ('no alignment or dist matrix')
            return
        import toyplot
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
        labelcol = self.labelsw.currentText()

        if not hasattr(self, 'snpdist'):
            print ('no distance matrix')
            return

        idx = list(self.sub.index)
        M = self.snpdist.loc[idx,idx]
        w = widgets.PlotWidget(self)
        tools.dist_matrix_to_mst(M, self.sub, colorcol, node_size=100, labelcol=labelcol,
                                  cmap=cmap, ax=w.ax)
        self.show_dock_object(w, 'mst')
        return

    def dist_matrix_from_alignment(self):
        """Calculate snp dist matrix from alignment"""

        if self.aln is None:
            print ('no alignment loaded')
            return
        from Bio.Align import MultipleSeqAlignment
        aln = MultipleSeqAlignment(self.aln)
        print ('calculating distance matrix for %s sequences..' %len(aln))

        def func(progress_callback):
            self.snpdist = tools.snp_dist_matrix(aln, threads=core.THREADS)
        self.run_threaded_process(func, self.processing_completed)
        return

    def show_filter(self):
        """Filter widget"""

        self.table_widget.showSearchBar()
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

        idx = self.tabs.currentIndex()
        if idx == 0:
            #get the current figure and make a copy of it by using pickle
            fig = self.plotview.fig
            p = pickle.dumps(fig)
            fig = pickle.loads(p)
            self.scratch_items[label] = fig
            if hasattr(self, 'scratchpad'):
                self.scratchpad.update(self.scratch_items)
        elif idx == 1:
            #save the image directly
            filename, _ = QFileDialog.getSaveFileName(self,"Save Screen Capture",
                                                    "","png files (*.png);;All files (*.*)"
                                                    )
            if filename:
                img = self.foliumview.screen_capture(filename)
        return

    def save_selection(self):
        """Save current selection as a shortcut"""

        idx = self.sub.index
        #get name
        label, ok = QInputDialog.getText(self, 'Name Selection', 'Enter a label:')
        if not ok:
            return

        d = self.selections[label] = {}
        d['indexes'] = idx
        #get widget settings
        d['options'] = widgets.getWidgetValues(self.widgets)
        #plot zoom
        d['zoom'] = self.plotview.get_plot_lims()
        #neighbours
        d['neighbours'] = self.neighbours
        if self.neighbours is not None:
            m = tools.get_dataframe_memory(self.neighbours)
            print (m)
        #print (d)
        self.update_selections_menu()
        return

    def clear_selections(self, ask=False):
        """Clear saved selections"""

        if ask == True:
            reply = QMessageBox.question(self, 'Confirm',
                                    "This will clear all saved selections.\nAre you sure?",
                                    QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.No:
                return
        self.selections = {}
        self.clear_selections_menu()
        return

    def clear_selections_menu(self):
        """Clear selections and menu items"""

        menu = self.selections_menu
        items = menu.actions()
        for action in items[3:]:
            menu.removeAction(action)
            action.deleteLater()
        return

    def update_selections_menu(self):
        """Update menu for selections"""

        from functools import partial
        menu = self.selections_menu
        self.clear_selections_menu()

        for name in self.selections:
            #print (name)
            def func(name):
                #show saved selection
                df = self.meta_table.model.df
                d = self.selections[name]
                idx = d['indexes']
                self.sub = df.loc[idx]
                self.plotview.lims = None
                #restore saved widget values
                if 'options' in d:
                    widgets.setWidgetValues(self.widgets, d['options'])
                #any other objects?
                for k in ['neighbours']:
                    if k in d:
                        self.__dict__[k] = d[k]
                self.update()
                return
            menu.addAction(name, partial(func, name))
        return

    def herd_summary(self):
        """Summary by herd"""

        res=[]
        df = self.meta_table.model.df
        for herd, sdf in df.groupby('HERD_NO'):
            #print (herd)
            clades = len(sdf.snp3.unique())
            m = self.moves[self.moves.move_to == herd]
            res.append([herd,len(sdf),clades,len(m)])
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

        w = tables.DataFrameTable(self, dataframe=self.sub,
                font=core.FONT, fontsize=core.FONTSIZE)
        self.show_dock_object(w, 'selected', 'left')
        return

    def show_moves_table(self, df):
        """Show moves for samples in separate table"""

        if df is None:
            df = pd.DataFrame()

        t = tables.MovesTable(self, dataframe=df,
                            font=core.FONT, fontsize=core.FONTSIZE, app=self)
        w = tables.DataFrameWidget(parent=self, table=t, toolbar=False)
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
        gpdver = gpd.__version__
        text='TracebTB\n'\
            +'version '+__version__+'\n'\
            +'Copyright (C) Damien Farrell 2022-\n'\
            +'This program is free software; you can redistribute it and/or '\
            +'modify it under the terms of the GNU GPL '\
            +'as published by the Free Software Foundation; either '\
            +'version 3 of the License, or (at your option) any '\
            +'later version.\n'\
            +'Using Python v%s, %s\n' %(pythonver, pyqtver)\
            +'pandas v%s, matplotlib v%s' %(pandasver,mplver)\
            +'geopandas v%s' %gpdver

        msg = QMessageBox.about(self, "About", text)
        return

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
