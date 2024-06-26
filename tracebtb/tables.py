#!/usr/bin/env python

"""
    dataframe table widget and sub-classes.
    Created Nov 2019
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

import sys,os,platform
import pandas as pd
import numpy as np
import pylab as plt
from .qt import *
from . import core, widgets, plotting
from pandas.api.types import is_datetime64_any_dtype as is_datetime

style = '''
    QWidget {
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
        max-width: 150px;
    }
    '''

class ColumnHeader(QHeaderView):
    def __init__(self):
        super(QHeaderView, self).__init__()
        return

class DataFrameWidget(QWidget):
    """Widget containing a tableview and statusbar"""
    def __init__(self, parent=None, table=None, statusbar=True, searchbar=False, toolbar=False, **kwargs):
        """
        Widget containing a dataframetable - allows us to pass any table subclass
        """

        super(DataFrameWidget, self).__init__()
        l = self.layout = QGridLayout()
        l.setSpacing(2)
        self.setLayout(self.layout)
        #self.plotview = widgets.PlotViewer()
        if table == None:
            self.table = DataFrameTable(self, dataframe=pd.DataFrame())
        else:
            self.table = table
            table.parent = self
        l.addWidget(self.table, 1, 1)
        if toolbar == True:
            self.createToolbar()
        if statusbar == True:
            self.statusBar()
        if searchbar == True:
            self.showSearchBar()

        self.table.model.dataChanged.connect(self.stateChanged)
        return

    #@Slot('QModelIndex','QModelIndex','int')
    def stateChanged(self, idx, idx2):
        """Run whenever table model is changed"""

        if hasattr(self, 'pf') and self.pf is not None:
            self.pf.updateData()

    def toggleSearchBar(self):

        if not hasattr(self, 'searchbar'):
            self.showSearchBar()
        else:
            self.searchbar.hide()
        return

    def showSearchBar(self):
        """Search bar"""

        df = self.table.model.df
        cols = list(df.columns)
        cols.insert(0,'Any')
        if hasattr(self, 'searchbar'):
            self.searchbar.show()
        else:
            sb = self.searchbar = widgets.SimpleFilterWidget(self, self.table)
            self.layout.addWidget(sb, 3, 1)
        return

    def statusBar(self):
        """Status bar at bottom"""

        w = self.statusbar = QWidget(self)
        l = QHBoxLayout(w)
        w.setMaximumHeight(30)
        self.size_label = QLabel("")
        l.addWidget(self.size_label, 1)
        w.setStyleSheet('color: #1a216c; font-size:12px')
        self.layout.addWidget(w, 2, 1)
        self.updateStatusBar()
        return

    def updateStatusBar(self):
        """Update the table details in the status bar"""

        if not hasattr(self, 'size_label'):
            return
        df = self.table.model.df
        #meminfo = self.table.getMemory()
        s = '{r} samples x {c} columns'.format(r=len(df), c=len(df.columns))
        self.size_label.setText(s)
        return

    def createToolbar(self):

        self.setLayout(self.layout)
        items = {
                 'filter':{'action':self.toggleSearchBar,'file':'filter'},
                 'bar': {'action':self.plot_bar,'file':'plot-bar'},
                 'barh': {'action':self.plot_barh,'file':'plot-barh'},
                 'hist': {'action':self.plot_hist,'file':'plot-hist'},
                 #'scatter': {'action':self.plot_scatter,'file':'plot-scatter'},
                 #'pie': {'action':self.plot_pie,'file':'plot-pie'},
                 }

        self.toolbar = toolbar = QToolBar("Toolbar")
        #toolbar.setIconSize(QtCore.QSize(core.ICONSIZE, core.ICONSIZE))
        toolbar.setOrientation(QtCore.Qt.Vertical)
        widgets.addToolBarItems(toolbar, self, items)
        self.layout.addWidget(toolbar,1,2)
        return

    def refresh(self):

        self.table.refresh()
        return

    def copy(self):
        """Copy to clipboard"""

        #check size of dataframe
        m = self.table.model.df.memory_usage(deep=True).sum()
        if m>1e8:
            answer = QMessageBox.question(self, 'Copy?',
                             'This data may be too large to copy. Are you sure?', QMessageBox.Yes, QMessageBox.No)
            if not answer:
                return
        df = self.table.getSelectedDataFrame()
        df.to_clipboard()
        return

    def plot_bar(self):
        """Plot from selection"""

        self.table.plot(kind='bar')
        return

    def plot_barh(self):
        """Plot from selection"""

        self.table.plot(kind='barh')
        return

    def plot_hist(self):
        """Plot from selection"""

        self.table.plot(kind='hist')
        return

    def plot_scatter(self):
        self.table.plot(kind='scatter')
        return

    def plot_line(self):
        self.table.plot(kind='line')
        return

    def plot_pie(self):
        self.table.plot(kind='pie')
        return

    def filter(self):
        """Show filter dialog"""

        return

class HeaderView(QHeaderView):
    """"
    Column header class.
    """
    def __init__(self, parent):
        super(HeaderView, self).__init__(QtCore.Qt.Horizontal, parent)
        '''self.setStyleSheet(
            "QHeaderView::section{background-color: #ffffff; "
            "font-weight: bold; "
            "border-bottom: 1px solid gray;}")'''

        self.setDefaultAlignment(QtCore.Qt.AlignLeft|QtCore.Qt.Alignment(QtCore.Qt.TextWordWrap))
        sizePol = QSizePolicy()
        sizePol.setVerticalPolicy(QSizePolicy.Maximum)
        sizePol.setHorizontalPolicy(QSizePolicy.Maximum)
        self.setSizePolicy(sizePol)
        self.MAX_HEIGHT = 240
        self.setMinimumHeight(26)
        self.setMaximumHeight(self.MAX_HEIGHT )
        self.setSectionsClickable(True)
        self.setSelectionBehavior(QTableView.SelectColumns)
        self.setStretchLastSection(False)
        return

    def sectionSizeFromContents(self, logicalIndex):
        """Get section size from contents"""

        text = self.model().headerData(logicalIndex, self.orientation(), QtCore.Qt.DisplayRole)
        alignment = self.defaultAlignment()
        metrics = QFontMetrics(self.fontMetrics())
        width = metrics.boundingRect(QtCore.QRect(), alignment, text).width()

        heights = []
        for i in range(self.count()):
            text = self.model().headerData(i, self.orientation(),QtCore.Qt.DisplayRole)
            size = self.sectionSize(i)
            rect = QtCore.QRect(0, 0, size, self.MAX_HEIGHT)
            heights.append(metrics.boundingRect(rect, alignment, text).height())
        height = sorted(heights)[-1] + 5
        return QtCore.QSize(width, height)

class DataFrameTable(QTableView):
    """
    QTableView with pandas DataFrame as model.
    """
    def __init__(self, parent=None, dataframe=None, plotter=None, font='Arial', fontsize=10):

        QTableView.__init__(self)
        self.parent = parent
        #self.clicked.connect(self.showSelection)
        #self.doubleClicked.connect(self.handleDoubleClick)
        #self.setSelectionBehavior(QTableView.SelectRows)
        #self.setSelectionBehavior(QTableView.SelectColumns)
        #self.horizontalHeader = ColumnHeader()

        vh = self.verticalHeader()
        vh.setVisible(True)
        vh.setDefaultSectionSize(28)
        vh.setMinimumWidth(20)
        vh.setMaximumWidth(500)
        vh.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        vh.customContextMenuRequested.connect(self.rowHeaderMenu)

        self.headerview = HeaderView(self)
        self.setHorizontalHeader(self.headerview)
        hh = self.horizontalHeader()
        hh.setVisible(True)
        hh.setSectionsMovable(True)
        hh.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        hh.customContextMenuRequested.connect(self.columnHeaderMenu)
        hh.setSelectionBehavior(QTableView.SelectColumns)
        hh.setSelectionMode(QAbstractItemView.ExtendedSelection)
        #hh.sectionClicked.connect(self.columnClicked)
        hh.setSectionsClickable(True)
        self.setDragEnabled(True)
        self.viewport().setAcceptDrops(True)
        self.setDropIndicatorShown(True)
        self.resizeColumnsToContents()
        self.setCornerButtonEnabled(True)

        self.font = QFont(font, fontsize)
        self.fontname = font
        self.fontsize = fontsize
        self.updateFont()

        tm = DataFrameModel(dataframe)
        #self.setModel(tm)
        self.model = tm
        self.setProxyModel()

        self.setWordWrap(False)
        self.setCornerButtonEnabled(True)
        if plotter == None:
            self.plotview = widgets.PlotViewer()
        else:
            self.plotview = plotter
        return

    def setProxyModel(self):
        """Set a proxy model for filtering"""

        self.proxy = QtCore.QSortFilterProxyModel()
        #self.proxy = CustomProxyModel()
        self.proxy.setSourceModel(self.model)
        self.setModel(self.proxy)
        self.proxy.sort(0, Qt.AscendingOrder)
        self.setSortingEnabled(True)

    def updateFont(self):
        """Update the font"""

        font = QFont(self.fontname)
        font.setPointSize(int(self.fontsize))
        self.setFont(font)
        self.horizontalHeader().setFont(font)
        self.verticalHeader().setFont(font)
        return

    def setDataFrame(self, df=None):
        """Set the dataframe. No argument will set an empty table"""

        if df is None:
            df = pd.DataFrame()
        tm = DataFrameModel(df)
        #self.setModel(tm)
        self.model = tm
        #reset the proxy model too
        self.setProxyModel()
        return

    def getDataFrame(self):
        return self.model.df

    def zoomIn(self, fontsize=None):
        """Zoom in"""

        if fontsize == None:
            s = self.font.pointSize()+1
        else:
            s = fontsize
        self.font.setPointSize(s)
        self.setFont(self.font)
        vh = self.verticalHeader()
        h = vh.defaultSectionSize()
        vh.setDefaultSectionSize(h+2)
        hh = self.horizontalHeader()
        w = hh.defaultSectionSize()
        hh.setDefaultSectionSize(w+2)
        return

    def zoomOut(self, fontsize=None):

        if fontsize == None:
            s = self.font.pointSize()-1
        else:
            s = fontsize
        self.font.setPointSize(s)
        self.setFont(self.font)
        vh = self.verticalHeader()
        h = vh.defaultSectionSize()
        vh.setDefaultSectionSize(h-2)
        hh = self.horizontalHeader()
        w = hh.defaultSectionSize()
        hh.setDefaultSectionSize(w-2)
        return

    def importFile(self, filename=None, dialog=False, **kwargs):

        if dialog is True:
            options = QFileDialog.Options()
            filename, _ = QFileDialog.getOpenFileName(self,"Import File",
                                                      "","All Files (*);;Text Files (*.txt);;CSV files (*.csv)",
                                                      options=options)
            df = pd.read_csv(filename, **kwargs)
            self.table.model.df = df
        return

    def info(self):

        buf = io.StringIO()
        self.table.model.df.info(verbose=True,buf=buf,memory_usage=True)
        td = dialogs.TextDialog(self, buf.getvalue(), 'Info')
        return

    '''def showSelection(self, item):
        """Show the row selection"""

        cellContent = item.data()
        row = item.row()
        model = item.model()
        index = model.createIndex(row, 0)
        columnsTotal = model.columnCount(index)
        return'''

    '''def getSelectedRows(self):
        """Get selected rows"""

        sm = self.selectionModel()
        rows = [(i.row()) for i in sm.selectedIndexes()]
        rows = list(dict.fromkeys(rows).keys())
        return rows'''

    def getSelectedRows(self):
        """Get selected rows. Uses proxy index."""

        sm = self.selectionModel()
        rows = []
        for proxy_index in sm.selectedIndexes():
            source_index = self.proxy.mapToSource(proxy_index)
            rows.append(source_index.row())
        # Remove duplicates and return unique rows
        rows = list(set(rows))
        return rows

    def getSelectedIndexes(self):
        """Get selected row indexes"""

        rows = self.getSelectedRows()
        idx = self.model.df.index[rows]
        return idx

    def getSelectedColumns(self):
        """Get selected column indexes"""

        sm = self.selectionModel()
        cols = [(i.column()) for i in sm.selectedIndexes()]
        cols = list(dict.fromkeys(cols).keys())
        return cols

    def getSelectedDataFrame(self):
        """Get selection as a dataframe"""

        df = self.model.df
        sm = self.selectionModel()

        rows = self.getSelectedRows()
        #rows = [(i.row()) for i in sm.selectedIndexes()]
        cols = [(i.column()) for i in sm.selectedIndexes()]
        #get unique rows/cols keeping order
        #rows = list(dict.fromkeys(rows).keys())
        cols = list(dict.fromkeys(cols).keys())
        return df.iloc[rows,cols]

    def setSelected(self, rows, cols):
        """
        Set selection programmatically from a list of rows and cols.
        https://doc.qt.io/archives/qtjambi-4.5.2_01/com/trolltech/qt/model-view-selection.html
        """

        #print (rows,cols)
        if len(rows)==0 or len(cols)==0:
            return
        topleft = self.model.index(rows[0], cols[0])
        bottomright = self.model.index(rows[-1], cols[-1])
        selection = QtCore.QItemSelection(topleft, bottomright)
        mode = QtCore.QItemSelectionModel.Select
        self.selectionModel().select(selection, mode)
        return

    def getScrollPosition(self):
        """Get current row/col position"""

        hb = self.horizontalScrollBar()
        vb = self.verticalScrollBar()
        return vb.value(),hb.value()

    def setScrollPosition(self, row, col):
        """Move to row/col position"""

        idx = self.model.index(row, col)
        self.scrollTo(idx)
        return

    def handleDoubleClick(self, item):

        cellContent = item.data()
        if item.column() != 0:
            return
        return

    def columnClicked(self, col):

        hheader = self.horizontalHeader()
        df = self.model.df
        return

    def storeCurrent(self):
        """Store current version of the table before a major change is made"""

        self.prevdf = self.model.df.copy()
        return

    def deleteCells(self, rows, cols, answer=None):
        """Clear the cell contents"""

        if answer == None:
            answer = QMessageBox.question(self, 'Delete Cells?',
                             'Are you sure?', QMessageBox.Yes, QMessageBox.No)
        if not answer:
            return
        self.storeCurrent()
        #print (rows, cols)
        self.model.df.iloc[rows,cols] = np.nan
        return

    def editCell(self, item):
        return

    def setColumnType(self, column=None):
        """Change the column dtype"""

        idx = self.getSelectedColumns()
        if len(idx)>0:
            cols = self.model.df.columns[idx]
        else:
            cols = [column]
        types = ['float','int','str','object','datetime64[ns]']
        newtype, ok = QInputDialog().getItem(self, "New Column Type",
                                             "Type:", types, 0, False)
        if not ok:
            return
        self.storeCurrent()
        df = self.model.df
        for c in cols:
            df[c] = df[c].astype(newtype)
        self.refresh()
        return

    def setRowColor(self, rowIndex, color):
        for j in range(self.columnCount()):
            self.item(rowIndex, j).setBackground(color)

    def rowHeaderMenu(self, pos):
        """Row header popup menu"""

        vheader = self.verticalHeader()
        idx = vheader.logicalIndexAt(pos)
        menu = QMenu(self)

        resetIndexAction = menu.addAction("Reset Index")
        sortIndexAction = menu.addAction("Sort By Index \u2193")
        sortIndexDescAction = menu.addAction("Sort By Index \u2191")
        setTypeAction = menu.addAction("Set Type")
        deleteAction = menu.addAction("Delete Rows")
        action = menu.exec_(self.mapToGlobal(pos))
        #if action == resetIndexAction:
        #    self.resetIndex()
        if action == sortIndexAction:
            self.sortIndex()
        elif action == sortIndexDescAction:
            self.sortIndex(ascending=False)
        elif action == setTypeAction:
            self.setIndexType()
        elif action == deleteAction:
            self.deleteRows()
        return

    def columnHeaderMenu(self, pos):

        hheader = self.horizontalHeader()
        idx = hheader.logicalIndexAt(pos)
        column = self.model.df.columns[idx]
        menu = QMenu(self)
        sortAction = menu.addAction("Sort \u2193")
        sortDescAction = menu.addAction("Sort \u2191")
        deleteColumnAction = menu.addAction("Delete Column")
        renameColumnAction = menu.addAction("Rename Column")
        addColumnAction = menu.addAction("Add Column")
        setTypeAction = menu.addAction("Set Data Type")
        #plotAction = menu.addAction("Histogram")
        colorbyAction = menu.addAction("Color By Column")

        action = menu.exec_(self.mapToGlobal(pos))
        if action == sortAction:
            self.sort(idx)
        elif action == sortDescAction:
            self.sort(idx, ascending=False)
        elif action == deleteColumnAction:
            self.deleteColumn(column)
        elif action == renameColumnAction:
            self.renameColumn(column)
        elif action == addColumnAction:
            self.addColumn()
        elif action == setTypeAction:
            self.setColumnType(column)
        #elif action == plotAction:
        #    self.plotHist(column)
        elif action == colorbyAction:
            self.colorByColumn(column)
        return

    def keyPressEvent(self, event):

        rows = self.getSelectedRows()
        cols = self.getSelectedColumns()
        if event.key() == QtCore.Qt.Key_Delete:
            self.deleteCells(rows, cols)

    def contextMenuEvent(self, event):
        """Reimplemented to create context menus for cells and empty space."""

        # Determine the logical indices of the cell where click occured
        hheader, vheader = self.horizontalHeader(), self.verticalHeader()
        position = event.globalPos()
        row = vheader.logicalIndexAt(vheader.mapFromGlobal(position))
        column = hheader.logicalIndexAt(hheader.mapFromGlobal(position))
        if row == -1:
            return
        # Show a context menu for empty space at bottom of table...
        self.menu = QMenu(self)
        self.addActions(event, row)
        return

    def addActions(self, event, row):
        """Actions"""

        menu = self.menu
        copyAction = menu.addAction("Copy")
        exportAction = menu.addAction("Export Table")
        transposeAction = menu.addAction("Transpose")
        action = menu.exec_(event.globalPos())
        if action == copyAction:
            self.copy()
        elif action == exportAction:
            self.exportTable()
        elif action == transposeAction:
            self.transpose()

        return

    def resetIndex(self):

        self.model.df.reset_index(inplace=True)
        self.refresh()
        return

    def setIndex(self, column):
        """Set column as index"""

        self.storeCurrent()
        self.model.df.set_index(column, inplace=True)
        self.refresh()
        return

    def sortIndex(self, ascending=True):
        """Sort by inde,"""

        self.model.df = self.model.df.sort_index(axis=0,ascending=ascending)
        self.refresh()
        return

    def setIndexType(self):
        """Set the type of the index"""

        types = ['float','int','object','datetime64[ns]']
        newtype, ok = QInputDialog().getItem(self, "New Type",
                                             "Type:", types, 0, False)
        if not ok:
            return
        self.storeCurrent()
        df = self.model.df
        #currtype = df.index.dtype
        df.index = df.index.astype(newtype)
        self.refresh()
        return

    def copy(self):
        """Copy selected cells"""

        df = self.getSelectedDataFrame()
        if len(df.columns==1):
            df.to_clipboard(index=False, header=False)
        else:
            df.to_clipboard(index=False, header=True)
        return

    def refresh(self):
        """Refresh table if dataframe is changed"""

        #self.updateFont()
        self.model.beginResetModel()
        index = self.model.index
        try:
            self.model.dataChanged.emit(0,0)
        except:
            self.model.dataChanged.emit(index(0,0),index(0,0))
        self.model.endResetModel()
        if hasattr(self.parent,'statusbar'):
            self.parent.updateStatusBar()
        return

    def importFile(self):
        dialogs.ImportDialog(self)
        return

    def exportTable(self):
        filename, _ = QFileDialog.getSaveFileName(self, "Export",
                                                  "","csv files (*.csv);;All files (*.*)")
        if filename:
            self.model.df.to_csv(filename)
        return

    def addColumn(self):
        """Add a  column"""

        df = self.model.df
        name, ok = QInputDialog().getText(self, "Enter Column Name",
                                             "Name:", QLineEdit.Normal)
        if ok and name:
            if name in df.columns:
                return
            df[name] = pd.Series()
            self.refresh()
        return

    def deleteColumn(self, column=None):

        reply = QMessageBox.question(self, 'Delete Columns?',
                             'Are you sure?', QMessageBox.Yes, QMessageBox.No)
        if reply == QMessageBox.No:
            return False
        self.model.df = self.model.df.drop(columns=[column])
        self.refresh()
        return

    def deleteRows(self):

        rows = self.getSelectedRows()
        reply = QMessageBox.question(self, 'Delete Rows?',
                             'Are you sure?', QMessageBox.Yes, QMessageBox.No)
        if reply == QMessageBox.No:
            return False
        idx = self.model.df.index[rows]
        self.model.df = self.model.df.drop(idx)
        self.refresh()
        return

    def renameColumn(self, column=None):

        name, ok = QInputDialog().getText(self, "Enter New Column Name",
                                             "Name:", QLineEdit.Normal)
        if ok and name:
            self.model.df.rename(columns={column:name},inplace=True)
            self.refresh()
        return

    def sort(self, idx, ascending=True):
        """Sort by selected columns"""

        df = self.model.df
        sel = self.getSelectedColumns()
        if len(sel)>1:
            for i in sel:
                self.model.sort(i, ascending)
        else:
            self.model.sort(idx, ascending)
        return

    def applyFilters(self, query):
        """Apply filter to table"""

        #print (query)
        #regexp = QtCore.QRegExp(query, Qt.CaseInsensitive, QtCore.QRegExp.FixedString)
        regexp = f".*{query}.*"
        self.proxy.setFilterRegExp(regexp)
        #self.proxy.setFilterFixedString(query)
        #self.proxy.setFilterKeyColumn(-1)
        #self.proxy.invalidateFilter()
        return

    def searchListFilter(self, items):
        """Apply multiple queries and return them all"""

        for query in items:
            regexp = f".*{query}.*"
        self.proxy.setFilterRegExp(regexp)
        return

    def clearFilters(self):
        """Clear filters"""

        self.proxy.setFilterFixedString("")
        return

    def transpose(self):

        self.model.df = self.model.df.T
        self.refresh()
        return

    def plot(self, kind='bar'):
        """Plot table selection"""

        df = self.model.df
        #idx = self.getSelectedColumns()
        #cols = df.columns[idx]
        #d = df[cols]
        data = self.getSelectedDataFrame()
        d = data._get_numeric_data()
        #print (d)
        xcol = d.columns[0]
        ycols = d.columns[1:]

        self.plotview.clear()
        ax = self.plotview.ax
        if kind == 'bar':
            d.plot(kind='bar',ax=ax)
        elif kind == 'barh':
            d.plot(kind='barh',ax=ax)
        elif kind == 'hist':
            d.plot(kind='hist',subplots=True,ax=ax)
        elif kind == 'scatter':
            d=d.dropna()
            d.plot(x=xcol,y=ycols,kind='scatter',ax=ax)
        elif kind == 'pie':
            d.plot(kind='pie',subplots=True,legend=False,ax=ax)
        #ax.set_title(col)
        plt.tight_layout()
        self.plotview.redraw()
        self.plotview.show()
        self.plotview.activateWindow()
        return

    def colorByColumn(self, col):
        """Set colorby column"""

        #cmap = 'Set1'
        df = self.model.df
        colors,colormap = plotting.get_color_mapping(df,col,seed=10)
        print (colors)
        self.model.rowcolors = colors
        return

    def getMemory(self):
        """Get memory info as string"""

        m = self.model.df.memory_usage(deep=True).sum()
        if m>1e5:
            m = round(m/1048576,2)
            units='MB'
        else:
            units='Bytes'
        s = "%s %s" %(m,units)
        return s

class DataFrameModel(QtCore.QAbstractTableModel):
    def __init__(self, dataframe=None, *args):
        super(DataFrameModel, self).__init__()
        if dataframe is None:
            self.df = pd.DataFrame()
        else:
            self.df = dataframe
        self.bg = '#F4F4F3'
        self.rowcolors = None
        return

    def update(self, df):
        self.df = df

    def rowCount(self, parent=QtCore.QModelIndex()):
        return len(self.df.index)

    def columnCount(self, parent=QtCore.QModelIndex()):
        return len(self.df.columns.values)

    def data(self, index, role=QtCore.Qt.DisplayRole):
        """Edit or display roles. Handles what happens when the Cells
        are edited or what appears in each cell.
        """

        i = index.row()
        j = index.column()
        #print (self.df.dtypes)
        #coltype = self.df.dtypes[j]
        coltype = self.df[self.df.columns[j]].dtype
        isdate = is_datetime(coltype)
        if role == QtCore.Qt.DisplayRole:
            value = self.df.iloc[i, j]
            if isdate:
                return value.strftime(core.TIMEFORMAT)
            elif type(value) != str:
                if type(value) in [float,np.float64] and np.isnan(value):
                    return ''
                elif type(value) == float:
                    return value
                else:
                    return (str(value))
            else:
                return '{0}'.format(value)
        elif (role == QtCore.Qt.EditRole):
            value = self.df.iloc[i, j]
            if type(value) is str:
                try:
                    return float(value)
                except:
                    return str(value)
            if np.isnan(value):
                return ''
        elif role == QtCore.Qt.BackgroundRole:
            if self.rowcolors != None:
                clr = self.rowcolors[i]
                #print (clr)
                return QColor(clr)
            else:
                return QColor(self.bg)

    def headerData(self, col, orientation, role=QtCore.Qt.DisplayRole):
        """What's displayed in the headers"""

        if len(self.df.columns) == 0:
            return
        if role == QtCore.Qt.DisplayRole:
            if orientation == QtCore.Qt.Horizontal:
                return str(self.df.columns[col])
            if orientation == QtCore.Qt.Vertical:
                value = self.df.index[col]
                if type( self.df.index) == pd.DatetimeIndex:
                    if not value is pd.NaT:
                        try:
                            return value.strftime(TIMEFORMAT)
                        except:
                            return ''
                else:
                    return str(value)
        return None

    def sort(self, idx, ascending=True):
        """Sort table by given column number """

        self.layoutAboutToBeChanged.emit()
        if len(self.df.columns) == 0:
            return
        col = self.df.columns[idx]
        self.df = self.df.sort_values(col, ascending=ascending)
        self.layoutChanged.emit()
        return

    def setData(self, index, value, role=QtCore.Qt.EditRole):
        """Set data upon edits"""

        #print (index)
        i = index.row()
        j = index.column()
        curr = self.df.iloc[i,j]
        #print (curr, value)
        self.df.iloc[i,j] = value
        #self.dataChanged.emit()
        return True

    #def onDataChanged(self):
        #print (self.df)
    #    return

    def setColumnColor(self, columnIndex, color):
        for i in range(self.rowCount()):
            self.item(i, columnIndex).setBackground(color)
        return

    def flags(self, index):
            return Qt.ItemIsSelectable|Qt.ItemIsEnabled|Qt.ItemIsEditable

class CustomProxyModel(QtCore.QSortFilterProxyModel):
    def __init__(self):
        super(CustomProxyModel, self).__init__()

    def sort(self, column, order=Qt.AscendingOrder):
        """Perform the sort operation on the source model"""

        super(CustomProxyModel, self).sort(column, order)
        return

    '''def filterAcceptsRow(self, source_row, source_parent):

        df = self.sourceModel().df
        filter_text = self.filterRegExp().pattern()
        print (filter_text)
        if filter_text.isEmpty():
            return True

        for col in range(len(self.df.columns)):
            data = df.iloc[source_row, col]
            if filter_text in str(data):
                return True
        return False'''

class SampleTableModel(DataFrameModel):
    """Samples table model class"""
    def __init__(self, dataframe=None, *args):

        DataFrameModel.__init__(self, dataframe)
        self.df = dataframe

class SampleTable(DataFrameTable):
    """
    QTableView with pandas DataFrame as model.
    """

    def __init__(self, parent=None, app=None, dataframe=None, plotter=None, **kwargs):
        DataFrameTable.__init__(self, **kwargs)
        self.parent = parent
        self.app = app
        self.setWordWrap(False)
        tm = SampleTableModel(dataframe)
        self.setModel(tm)
        self.setProxyModel()
        return

    def setDataFrame(self, df=None):
        """Override to use right model"""

        if df is None:
            df = pd.DataFrame()
        if 'sample' in df.columns:
            df = df.set_index('sample', drop=False)
            df.index.name = 'index'
        tm = SampleTableModel(df)
        self.model = tm
        self.setProxyModel()
        return

    def addActions(self, event, row):
        """Right click menu"""

        menu = self.menu
        selectAction = menu.addAction("Select Samples")
        addSelectionAction = menu.addAction("Add to Selection")
        findRelatedAction = menu.addAction("Related Samples")
        copyAction = menu.addAction("Copy")
        detailsAction = menu.addAction("Sample Details")
        removeAction = menu.addAction("Delete Selected")
        exportAction = menu.addAction("Export Table")
        action = menu.exec_(self.mapToGlobal(event.pos()))

        rows = self.getSelectedRows()
        df = self.model.df.iloc[rows[0]]
        if action == detailsAction:
            self.app.sample_details(df)
        elif action == removeAction:
            self.deleteRows(rows)
        elif action == copyAction:
            self.copy()
        elif action == exportAction:
            self.exportTable()
        elif action == selectAction:
            self.app.selection_from_table()
        elif action == addSelectionAction:
            self.app.add_to_selection()
        elif action == findRelatedAction:
            self.app.select_related()
        return

    def edit(self, index, trigger, event):
        """Override edit to disable editing of columns"""

        if index.column() < 20:
            return False
        else:
            QTableView.edit(self, index, trigger, event)
        return True

    def resizeColumns(self):

        self.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        #self.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        return

    def deleteRows(self, rows):
        """Delete rows in sample table"""

        answer = QMessageBox.question(self, 'Delete Entry?',
                             'Delete these rows. Are you sure?',
                             QMessageBox.Yes, QMessageBox.No)
        if answer == QMessageBox.No:
            return
        idx = self.model.df.index[rows]
        self.model.df = self.model.df.drop(idx)
        #also sync the geodataframe
        #mask = ~self.app.meta_table.index.isin(idx)
        #self.app.cent = self.app.cent[mask]
        #print (len(self.app.cent))
        self.app.update_groups()
        self.refresh()
        return

    def deleteColumn(self, cols):
        """Delete columns in sample table"""

        answer = QMessageBox.question(self, 'Delete Columns?',
                             'Delete these columns. Are you sure?',
                             QMessageBox.Yes, QMessageBox.No)
        if answer == QMessageBox.No:
            return
        self.model.df = self.model.df.drop(columns=cols)
        #also sync the geodataframe
        #self.app.cent = self.app.cent.drop(columns=cols)
        self.refresh()
        return

class SelectedModel(DataFrameModel):
    """Moves model class"""
    def __init__(self, dataframe=None, *args):

        DataFrameModel.__init__(self, dataframe)
        self.df = dataframe

class SelectedTable(DataFrameTable):
    """
    Sub selection table.
    """
    def __init__(self, parent=None, app=None, dataframe=None, **kwargs):
        DataFrameTable.__init__(self, dataframe=dataframe, **kwargs)
        self.parent = parent
        self.app = app
        tm = DataFrameModel(dataframe)
        self.setModel(tm)
        self.setProxyModel()
        #print (self.model)
        return

    def addActions(self, event, row):
        """Right click menu"""

        menu = self.menu
        #detailsAction = menu.addAction("Sample Details")
        findRelatedAction = menu.addAction("Related Samples")
        removeSelectedAction = menu.addAction("Remove from Selection")
        copyAction = menu.addAction("Copy")
        exportAction = menu.addAction("Export Table")
        action = menu.exec_(self.mapToGlobal(event.pos()))

        rows = self.getSelectedRows()
        df = self.model.df.iloc[rows[0]]
        if action == copyAction:
            self.copy()
        elif action == exportAction:
            self.exportTable()
        elif action == findRelatedAction:
            idx = self.getSelectedIndexes()
            self.app.select_related(idx)
        elif action == removeSelectedAction:
            idx = self.getSelectedIndexes()
            self.app.remove_from_selection(idx)
        return

class MovesTableModel(DataFrameModel):
    """Moves model class"""
    def __init__(self, dataframe=None, *args):

        DataFrameModel.__init__(self, dataframe)
        self.df = dataframe
        self.setRowColors()

    def setRowColors(self):
        """Set row colors for rows based on index"""

        self.row_colors={}
        df = self.df
        for i, idx in enumerate(df.index.unique()):
            #print (i,idx)
            if i % 2 == 0:
                self.row_colors[idx] = QColor("#FAFAFA")
            else:
                self.row_colors[idx] = QColor("#EEEEEE")
        return

    def data(self, index, role=QtCore.Qt.DisplayRole):
        """Override data method"""

        parent_data = super().data(index, role)
        i = index.row()
        j = index.column()
        idx = self.df.index[i]
        if role == QtCore.Qt.DisplayRole:
            return parent_data
        if role == QtCore.Qt.BackgroundRole:
            row_color = self.row_colors[idx]
            colname = self.df.columns[j]
            value = self.df.iloc[i, j]
            if value == None:
                clr = '#FA9B8D'
                return QColor(clr)
            if colname == 'data_type' and value == 'F_to_F':
                clr = '#9FBAFA' #'#FA9B8D'
                return QColor(clr)
            else:
                return QColor(row_color)
            return

class MovesTable(DataFrameTable):
    """
    Moves table.
    """
    def __init__(self, parent=None, app=None, dataframe=None, **kwargs):
        DataFrameTable.__init__(self, **kwargs)
        self.parent = parent
        self.app = app
        self.setWordWrap(False)
        tm = MovesTableModel(dataframe)
        self.setModel(tm)
        self.model = tm
        #self.setProxyModel()
        return

    def addActions(self, event, row):

        menu = self.menu
        showmovedAction = menu.addAction("Show Farm Moves Only")
        showmovedAction = menu.addAction("Show All")
        copyAction = menu.addAction("Copy")
        exportAction = menu.addAction("Export Table")
        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == exportAction:
            self.exportTable()
        elif action == showmovedAction:
            self.showMoved()
        elif action == copyAction:
            self.copy()
        return

    def rowHeaderMenu(self, pos):
        """Row header popup menu"""

        vheader = self.verticalHeader()
        idx = vheader.logicalIndexAt(pos)
        menu = QMenu(self)

    def showMoved(self):

        df = self.model.df
        self.proxy.setFilterKeyColumn(6)
        self.applyFilters('F_to_F')
        return

class HerdTable(DataFrameTable):
    """
    QTableView for herd data.
    """
    def __init__(self, parent=None, app=None, dataframe=None, **kwargs):
        DataFrameTable.__init__(self, **kwargs)
        self.parent = parent
        self.app = app
        self.setWordWrap(False)
        tm = DataFrameModel(dataframe)
        self.setModel(tm)
        self.model = tm
        self.setProxyModel()
        return

    def addActions(self, event, row):

        menu = self.menu
        plotherdsAction = menu.addAction("Show Herd")
        exportAction = menu.addAction("Export Table")
        action = menu.exec_(self.mapToGlobal(event.pos()))
        if action == exportAction:
            self.exportTable()
        elif action == plotherdsAction:
            df=self.model.df
            rows = self.getSelectedRows()
            idx = df.index[rows]
            df = self.model.df.loc[idx]
            self.app.plot_herd_selection(df.HERD_NO)
        return

class TableViewer(QDialog):
    """
    View row of data in table
    """
    def __init__(self, parent=None, dataframe=None, **kwargs):
        super(TableViewer, self).__init__(parent)
        self.setGeometry(QtCore.QRect(200, 200, 600, 600))
        self.grid = QGridLayout()
        self.setLayout(self.grid)
        self.table = DataFrameTable(self, dataframe, **kwargs)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.grid.addWidget(self.table)
        return

    def setDataFrame(self, dataframe):
        self.table.model.df = dataframe
        return