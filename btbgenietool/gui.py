#!/usr/bin/env python

"""
    btbgenietool GUI.
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
from .qt import *
import pandas as pd
import numpy as np
import pylab as plt
from Bio import SeqIO
import matplotlib as mpl
from . import widgets, tables, phylo
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from matplotlib_scalebar.scalebar import ScaleBar

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
logoimg = os.path.join(module_path, 'logo.svg')
stylepath = os.path.join(module_path, 'styles')
iconpath = os.path.join(module_path, 'icons')

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

class App(QMainWindow):
    """GUI Application using PySide2 widgets"""
    def __init__(self, project=None, tree=None):

        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("BTBGenIE tool")

        self.setWindowIcon(QIcon(logoimg))
        self.create_menu()
        self.main = QSplitter(self)
        screen_resolution = QGuiApplication.primaryScreen().availableGeometry()
        width, height = screen_resolution.width()*0.9, screen_resolution.height()*.8
        if screen_resolution.width()>1920:
            self.setGeometry(QtCore.QRect(100, 100, width, height))
        else:
            self.showMaximized()
        self.setMinimumSize(400,300)

        self.recent_files = ['']
        self.main.setFocus()
        self.setCentralWidget(self.main)
        self.create_tool_bar()
        self.setup_gui()
        #self.load_settings()
        #self.show_recent_files()

        self.new_project()
        self.running = False
        self.load_test()

        if platform.system() == 'Windows':
            app.fetch_binaries()
        #if project != None:
        #    self.load_project(project)
        self.threadpool = QtCore.QThreadPool()
        return

    def create_tool_bar(self):
        """Create main toolbar"""

        items = {'New project': {'action': lambda: self.new_project(ask=True),'file':'document-new'},
                 'Open': {'action':self.load_project,'file':'document-open'},
                 'Save': {'action': lambda: self.save_project(),'file':'save'},
                 'Zoom out': {'action':self.zoom_out,'file':'zoom-out'},
                 'Zoom in': {'action':self.zoom_in,'file':'zoom-in'},
                 'Update map': {'action':self.replot,'file':'plot-map'},
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

    def create_controls(self):
        """controls for mapping"""

        m = QWidget()
        #select cluster
        clusts = []
        w = QComboBox(main)
        w.addItems(clusts)
        #select clust level

        #zoom to county
        counties = []
        w = QComboBox(main)
        w.addItems(counties)

        return m

    def setup_gui(self):
        """Add all GUI elements"""

        self.docks = {}
        self.main = main = QWidget()
        self.main.setFocus()
        self.setCentralWidget(self.main)
        l = QVBoxLayout(main)

        self.meta_table = tables.SampleTable(self, dataframe=pd.DataFrame, app=self)
        t = self.table_widget = tables.DataFrameWidget(parent=self, table=self.meta_table,
                            toolbar=True)
        self.add_dock(self.table_widget, 'meta data', scrollarea=False)

        self.tabs = QTabWidget(main)
        l.addWidget(self.tabs)
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self.close_tab)

        self.info = widgets.Editor(main, readOnly=True, fontsize=11)
        self.tabs.addTab(self.info, 'log')

        self.plotview = widgets.PlotViewer(self)
        idx = self.tabs.addTab(self.plotview, 'map')
        self.tabs.setCurrentIndex(idx)

        self.treeviewer = phylo.TreeViewer(self)
        #self.add_dock(self.treeviewer, 'phylogeny', 'right')

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
        return

    @QtCore.Slot(int)
    def close_tab(self, index):
        """Close current tab"""

        #index = self.tabs.currentIndex()
        #name = self.tabs.tabText(index)
        #self.tabs.removeTab(index)
        return

    def create_menu(self):
        """Create the menu bar for the application. """

        self.file_menu = QMenu('&File', self)

        icon = QIcon(os.path.join(iconpath,'document-new.png'))
        self.file_menu.addAction(icon, '&New Project', lambda: self.new_project(ask=True),
                QtCore.Qt.CTRL + QtCore.Qt.Key_N)
        icon = QIcon(os.path.join(iconpath,'document-open.png'))
        self.file_menu.addAction(icon, '&Open Project', self.load_project_dialog,
                QtCore.Qt.CTRL + QtCore.Qt.Key_O)
        self.recent_files_menu = QMenu("Recent Projects",
            self.file_menu)
        self.file_menu.addAction(self.recent_files_menu.menuAction())
        icon = QIcon(os.path.join(iconpath,'save.png'))
        self.file_menu.addAction(icon, '&Save Project', self.save_project,
                QtCore.Qt.CTRL + QtCore.Qt.Key_S)
        self.file_menu.addAction('&Save Project As', self.save_project_dialog)
        icon = QIcon(os.path.join(iconpath,'application-exit.png'))
        self.file_menu.addAction(icon, '&Quit', self.quit,
                QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        self.view_menu = QMenu('&View', self)
        self.menuBar().addMenu(self.view_menu)
        icon = QIcon(os.path.join(iconpath,'zoom-in.png'))
        self.view_menu.addAction(icon, 'Zoom In', self.zoom_in,
                QtCore.Qt.CTRL + QtCore.Qt.Key_Equal)
        icon = QIcon(os.path.join(iconpath,'zoom-out.png'))
        self.view_menu.addAction(icon, 'Zoom Out', self.zoom_out,
                QtCore.Qt.CTRL + QtCore.Qt.Key_Minus)

        self.tools_menu = QMenu('Tools', self)
        self.menuBar().addMenu(self.tools_menu)

        self.settings_menu = QMenu('&Settings', self)
        self.menuBar().addMenu(self.settings_menu)
        #self.settings_menu.addAction('Set Output Folder', self.set_output_folder)

        self.help_menu = QMenu('&Help', self)
        self.menuBar().addMenu(self.help_menu)
        self.help_menu.addAction('&Help', self.online_documentation)
        self.help_menu.addAction('&About', self.about)

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

    def save_project(self):
        """Save project"""

        if self.proj_file == None:
            self.save_project_dialog()

        filename = self.proj_file
        data={}
        data['inputs'] = self.meta_table.getDataFrame()
        keys = ['outputdir']
        for k in keys:
            if hasattr(self, k):
                data[k] = self.__dict__[k]

        if hasattr(self, 'treeviewer'):
            d=self.treeviewer.saveData()
            data['treeviewer'] = d
        self.opts.applyOptions()
        data['options'] = self.opts.kwds
        self.projectlabel.setText(filename)
        pickle.dump(data, open(filename,'wb'))
        self.add_recent_file(filename)
        return

    def save_project_dialog(self):
        """Save as project"""

        options = QFileDialog.Options()
        filename, _ = QFileDialog.getSaveFileName(self,"Save Project",
                                                  "","Project files (*.snipgenie);;All files (*.*)",
                                                  options=options)
        if filename:
            if not os.path.splitext(filename)[1] == '.snipgenie':
                filename += '.snipgenie'
            self.proj_file = filename
            self.save_project()
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
        self.sheets = {}
        self.proj_file = None
        #self.meta_table.setDataFrame(pd.DataFrame({'name':[]}))
        #self.left_tabs.clear()
        if hasattr(self, 'treeviewer'):
            self.treeviewer.clear()

        return

    def load_project(self, filename=None):
        """Load project"""

        self.new_project()
        data = pickle.load(open(filename,'rb'))
        keys = ['sheets','outputdir','results','ref_genome','ref_gb','mask_file']
        for k in keys:
            if k in data:
                self.__dict__[k] = data[k]

        if 'options' in data:
            self.opts.updateWidgets(data['options'])
        ft = self.meta_table
        ft.setDataFrame(data['inputs'])
        ft.resizeColumns()

        self.proj_file = filename
        self.projectlabel.setText(self.proj_file)
        self.outdirLabel.setText(self.outputdir)

        self.tabs.setCurrentIndex(0)
        self.add_recent_file(filename)
        return

    def load_project_dialog(self):
        """Load project"""

        filename, _ = QFileDialog.getOpenFileName(self, 'Open Project', './',
                                        filter="Project Files(*.snipgenie);;All Files(*.*)")
        if not filename:
            return
        if not os.path.exists(filename):
            print ('no such file')
        self.load_project(filename)
        return

    def load_test(self):
        """Load test dataset"""

        #metadata
        df = pd.read_csv('testing/samples.csv')
        df.set_index('sample',inplace=True)
        t = self.meta_table
        t.setDataFrame(df)
        t.resizeColumns()

        #reference map
        self.map = gpd.read_file('testing/counties.shp').to_crs("EPSG:29902")
        #centroids
        self.cent = gpd.read_file('testing/centroids.shp')
        self.cent['snp12'] = self.cent.snp12.astype(str)
        print (self.cent)
        self.replot()

        #tree
        tv = self.treeviewer
        tv.load_tree('testing/tree.newick')
        tv.style['layout'] = 'c'
        tv.style['tip_labels'] = False
        tv.set_zoom(2)
        #tv.update()
        #movement

        return

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

    def sample_details(self, row):
        return

    def replot(self):


        #data = self.cent[self.cent['sample']]
        fig,ax = plt.subplots(1,1)
        self.map.plot(edgecolor='black',color='#F6F4F3',ax=ax)
        #self.cent.plot(ax=ax)
        col='snp12'
        g = get_cluster_samples(self.cent,col=col)
        plot_clusters(g,col,ax=ax,legend=True)
        #ax.set_xlim((290000,351822))
        #ax.set_ylim((170000,240000))
        ax.add_artist(ScaleBar(dx=1,location=3))
        ax.axis('off')
        fig.tight_layout()
        self.plotview.set_figure(fig)
        return

    def plot_selected(self):

        df = self.meta_table.model.df
        self.rows = self.meta_table.getSelectedRows()
        self.sub = df.index[rows]
        self.replot()
        return

    def plot_snp_matrix(self):

        mat = pd.read_csv(self.results['snp_dist'],index_col=0)
        bv = widgets.BrowserViewer()
        import toyplot
        min=mat.min().min()
        max=mat.max().max()
        colormap = toyplot.color.brewer.map("BlueGreen", domain_min=min, domain_max=max)
        locator = toyplot.locator.Explicit(range(len(mat)),list(mat.index))
        canvas,axes = toyplot.matrix((mat.values,colormap), llocator=locator, tlocator=locator,
                        label="SNP distance matrix", colorshow=True)
        toyplot.html.render(canvas, "temp.html")
        with open('temp.html', 'r') as f:
            html = f.read()
            bv.browser.setHtml(html)

        idx = self.tabs.addTab(bv, 'snp dist')
        self.tabs.setCurrentIndex(idx)
        return

    def show_browser_tab(self, link, name):

        from PySide2.QtWebEngineWidgets import QWebEngineView
        browser = QWebEngineView()
        browser.setUrl(link)
        idx = self.tabs.addTab(browser, name)
        self.tabs.setCurrentIndex(idx)
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

    def zoom_in(self):

        w = self.meta_table
        w.zoomIn()
        self.info.zoomIn()
        return

    def zoom_out(self):

        w = self.meta_table
        w.zoomOut()
        self.info.zoomOut()
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
        #self.save_settings()
        event.accept()
        return

    def online_documentation(self,event=None):
        """Open the online documentation"""

        #import webbrowser
        link='https://github.com/dmnfarrell/btbgenietools'
        #webbrowser.open(link,autoraise=1)
        from PySide2.QtWebEngineWidgets import QWebEngineView
        browser = QWebEngineView()
        browser.setUrl(link)
        idx = self.tabs.addTab(browser, 'help')
        self.tabs.setCurrentIndex(idx)
        return

    def about(self):

        from . import __version__
        import matplotlib
        import PySide2
        pandasver = pd.__version__
        pythonver = platform.python_version()
        mplver = matplotlib.__version__
        qtver = PySide2.QtCore.__version__

        text='BTBGenIE tool\n'\
            +'version '+__version__+'\n'\
            +'Copyright (C) Damien Farrell 2022-\n'\
            +'This program is free software; you can redistribute it and/or '\
            +'modify it under the terms of the GNU GPL '\
            +'as published by the Free Software Foundation; either '\
            +'version 3 of the License, or (at your option) any '\
            +'later version.\n'\
            +'Using Python v%s, PySide2 v%s\n' %(pythonver, qtver)\
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

def get_cluster_samples(cent,col,n=3):
    #get samples in relevant clusters with more than n members
    v=cent[col].value_counts()
    v=v[v>=n]
    #print (col, v.index)
    clusts=v.index
    g = cent[cent[col]!='-1']
    g = g[g[col].isin(clusts)]
    g = g.sort_values(col)
    #print (col, len(g))
    print (list(g[col].unique()))
    return g

def get_clusts_info(cent,col):
    from shapely.geometry import Polygon
    new=[]
    for i,df in cent.groupby(col):
        if len(df)<3:
            continue
        d=df.dissolve(by=col)
        d['geometry'] = d.geometry.centroid
        d['area'] = df.dissolve(by=col).envelope.area/1e6
        d['animals'] = len(df)
        new.append(d[['geometry','area','animals']])
    return pd.concat(new).reset_index().sort_values('animals',ascending=False)

def plot_clusters(df,col=None,xlim=None,ylim=None,legend=True,title='',colors=None,
                  ms=None,cmap='Paired',ax=None):

    #mon.plot(linewidth=0.8, ax=ax, edgecolor='black',legend=legend,color='white', lw=1,alpha=0.5)
    if ms==None:
        df['ms'] = (df.species.astype('category').cat.codes+1)*60
    else:
        df['ms'] = ms
    #df['shape'] = df.Species.apply(lambda x: 'o' if x=='Bovine' else 's')
    #print (col,len(df))
    if colors is not None:
        df.plot(c=colors,ax=ax,alpha=0.8, markersize=df['ms'],edgecolor='0.1',marker='o',lw=0,
                  legend=legend)#,legend_kwds={'fontsize':9})
    else:
        df.plot(column=col,ax=ax,alpha=0.8, markersize=df['ms'],edgecolor='0.1',marker='o',lw=0,
              cmap=cmap,legend=legend)#,legend_kwds={'fontsize':9})
    if col != None:
        cent_cl = get_clusts_info(df,col)
        for x, y, label in zip(cent_cl.geometry.x, cent_cl.geometry.y, cent_cl[col]):
            ax.annotate(label, xy=(x, y), xytext=(0, 0), textcoords="offset points")

    '''minx, miny, maxx, maxy = counties.total_bounds
    if xlim!=None:
        ax.set_xlim(xlim[0], xlim[1])
    if ylim!=None:
       xlim
    ax.set_title(title,fontsize=20)'''
    ax.axis('off')
    return

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='btbgenie tool')
    parser.add_argument("-t", "--tree", dest="tree",default=[],
                        help="newick tree", metavar="FILE")
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
