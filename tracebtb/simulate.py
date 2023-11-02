# -*- coding: utf-8 -*-

"""
    Widgets for simulation of data with btbabm.
    Created Sep 2023
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

import sys, os, io, platform
import numpy as np
import pandas as pd
import pylab as plt
import matplotlib as mpl
import string
from .qt import *
from . import core, tools, widgets, tables
import geopandas as gpd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, AlignIO

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

module_path = os.path.dirname(os.path.abspath(__file__))
iconpath = os.path.join(module_path, 'icons')
data_path = os.path.join(module_path,'data')
borders = gpd.read_file(os.path.join(data_path,'counties.shp'))
logoimg = os.path.join(iconpath, 'simulate.svg')

class SimulateApp(QWidget):
    """Simulation dialog"""

    def __init__(self, parent=None):
        super(SimulateApp, self).__init__(parent)
        self.parent = parent
        self.setMinimumSize(400,400)
        self.setGeometry(QtCore.QRect(300, 200, 640, 600))
        self.setWindowTitle("Herd ABM Simulation")
        self.createWidgets()
        sizepolicy = QSizePolicy()
        self.setSizePolicy(QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding))
        self.setWindowIcon(QIcon(logoimg))
        self.path = None
        self.threadpool = QtCore.QThreadPool()
        return

    def createWidgets(self):
        """Create widgets. Plot on left and dock for tools on right."""

        layout = QFormLayout(self)

        row = QWidget()
        layout.addRow(row)
        hbox = QHBoxLayout(row)
        lbl = QLabel('Choose Folder to Save Results:')
        lbl.setMinimumWidth(100)
        hbox.addWidget(lbl)
        self.pickbtn = btn = QPushButton('..')
        btn.clicked.connect(lambda: self.pick_folder('.'))
        hbox.addWidget(btn)

        row = QWidget()
        layout.addRow(row)
        hbox = QHBoxLayout(row)
        hbox.setAlignment(QtCore.Qt.AlignLeft)
        lbl = QLabel('Steps to Run:')
        hbox.addWidget(lbl)
        w = self.stepsw = QSpinBox()
        w.setRange(10,1000)
        w.setValue(300)
        w.setSingleStep(10)
        hbox.addWidget(w)
        lbl = QLabel('Farms:')
        hbox.addWidget(lbl)
        w = self.farmsw = QSpinBox()
        w.setRange(10,2000)
        w.setValue(80)
        w.setSingleStep(10)
        hbox.addWidget(w)
        lbl = QLabel('Animals:')
        hbox.addWidget(lbl)
        w = self.animalsw = QSpinBox()
        w.setRange(10,5000)
        w.setValue(800)
        w.setSingleStep(10)
        hbox.addWidget(w)
        lbl = QLabel('Seed:')
        hbox.addWidget(lbl)
        w = self.seedw = QSpinBox()
        w.setRange(1,10000)
        w.setValue(12)
        w.setSingleStep(1)
        hbox.addWidget(w)
        lbl = QLabel('Allow Moves:')
        hbox.addWidget(lbl)
        w = self.movesw = QCheckBox()
        hbox.addWidget(w)

        w = self.tabs = QTabWidget(self)
        layout.addRow(w)
        w = self.gridview = widgets.PlotWidget(self)
        self.gridview.ax.axis('off')
        self.tabs.addTab(w, 'Network')
        w = self.parcelsview = widgets.PlotWidget(self)
        self.parcelsview.ax.axis('off')
        self.tabs.addTab(w, 'Parcels')
        w = self.plotview = widgets.PlotWidget(self)
        self.plotview.ax.axis('off')
        self.tabs.addTab(w, 'Plots')
        w = self.movesview = tables.DataFrameTable(self, dataframe=pd.DataFrame())        
        self.tabs.addTab(w, 'Moves')  

        self.progressbar = QProgressBar()
        self.progressbar.setRange(0,1)
        layout.addRow(self.progressbar)
        self.progressbar.setAlignment(Qt.AlignRight)
        self.progressbar.setMaximum(100)

        button = QPushButton("Run")
        button.clicked.connect(self.run)
        layout.addRow(button)
        self.stopbtn = button = QPushButton("Stop")
        button.clicked.connect(self.stop)
        layout.addRow(button)
        self.resumebtn = button = QPushButton("Resume")
        button.clicked.connect(self.resume)
        layout.addRow(button)
        button = QPushButton("Write Results")
        button.clicked.connect(self.write_results)
        layout.addRow(button)
        w = self.statusbar = QLabel()
        layout.addRow(w)
        w.setText('')
        return

    def pick_folder(self, name):
        """Get the folder to save to"""

        path = QFileDialog.getExistingDirectory(self,"Select folder",
                                                  #os.path.expanduser("~"),
                                            os.getcwd(),
                                            QFileDialog.ShowDirsOnly)
        if not path:
            return
        self.path = path
        self.pickbtn.setText(self.path)
        return

    def processing_completed(self):
        """Generic process completed"""

        #self.progressbar.setRange(0,100)
        self.running = False
        print ('finished')
        return

    def stop(self):
        self.running = False
        return

    def resume(self):
        self.run(self.model)

    def clear(self):
        self.plotview.ax.clear()
        self.parcelsview.ax.clear()

    def run(self, event=None, model=None):
        """Run simulation using btbabm"""

        import btbabm
        self.clear()
        steps = self.stepsw.value()
        farms = self.farmsw.value()
        animals = self.animalsw.value()
        seed = self.seedw.value()
        moves = self.movesw.isChecked()
        self.running = True
        self.progressbar.setValue(0)
        self.stopbtn.setStyleSheet("background-color: red;")
        if model is None:
            self.statusbar.setText('creating model..')
            #create land parcels
            parcels = btbabm.generate_land_parcels(1000,farms,empty=.5,
                                                 fragments=2,fragmented_farms=10,seed=seed)
            parcels['loc_type'] = 'farm'
            self.parcels = parcels
            parcels.plot(color=parcels.color,ec='.2',alpha=0.5,ax=self.parcelsview.ax)
            self.parcelsview.redraw()
            Glp,pos,cent = btbabm.land_parcels_to_graph(parcels,dist=100,empty=.5)
            #model = btbabm.FarmPathogenModel(farms,animals,0,
            #                                    #graph_seed=seed,
            #                                    seq_length=100,
            #                                    allow_moves=moves)
            model = btbabm.FarmPathogenModel(graph=Glp,pos=pos,C=animals,
                                     seq_length=100,
                                     allow_moves=moves,
                                     mean_inf_time=150,mean_stay_time=200,
                                     mean_latency_time=120,infected_start=5)
                                     
        self.model = model

        def func(progress_callback):
            print (self.model)
            self.statusbar.setText('running..')
            for s in range(steps):
                if self.running == False:
                    break
                model.step()
                progress_callback.emit(s/steps)
                ax = self.gridview.ax
                ax.clear()
                btbabm.plot_grid(model,ax=ax,pos=model.pos,colorby='strain',ns='num_infected')
                plt.tight_layout()
                self.gridview.redraw()

        self.run_threaded_process(func, self.processing_completed)
        return

    def write_results(self):

        if self.path == None:
            print ('no folder selected')
            return
        model = self.model
        #animal level data
        meta = model.get_animal_data(removed=True,infected=True)

        #m = pd.read_csv(os.path.join(self.path,'meta.csv'))
        #gdf = gpd.read_file(os.path.join(self.path,'gdf.shp'))
        #gdf = gdf.rename(columns={'id':'Animal_ID','herd':'HERD_NO'})
        meta.to_csv(os.path.join(self.path,'meta.csv'))

        model.make_phylogeny(removed=True,redundant=False)
        aln = AlignIO.read('temp.fasta','fasta')
        snpdist = btbabm.snp_dist_matrix(aln)
        snpdist.to_csv(os.path.join(self.path,'snpdist.csv'))
        #moves
        moves = model.get_moves()
        moves.to_csv(os.path.join(self.path,'mmovement.csv'))
        return

    def progress_fn(self, val):
        """Progress function"""

        #print (int(val*100))
        self.progressbar.setValue(int(val*100))
        return

    def run_threaded_process(self, process, on_complete):
        """Execute a function in the background with a worker"""

        worker = widgets.Worker(fn=process)
        self.threadpool.start(worker)
        worker.signals.finished.connect(on_complete)
        worker.signals.progress.connect(self.progress_fn)
        return

    def processing_completed(self):
        """Generic process completed"""

        model = self.model
        df = model.get_column_data()
        df.plot(ax=self.plotview.ax)
        self.plotview.redraw()
        self.movesview.setDataFrame(model.get_moves())
        print ('finished')

        self.stopbtn.setStyleSheet("")
        self.progressbar.setValue(100)
        return

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='ABM Simulator')

    args = vars(parser.parse_args())
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_ShareOpenGLContexts)
    app = QApplication(sys.argv)
    aw = SimulateApp(**args)
    aw.show()
    app.exec_()

if __name__ == '__main__':
    main()