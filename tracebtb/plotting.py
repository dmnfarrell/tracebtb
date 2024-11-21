"""
    Plotting methods
    Created Jan 2020
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

import os, sys, io, random, subprocess, time
import string
import numpy as np
import pandas as pd
pd.set_option('display.width', 200)

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO
import pylab as plt
import matplotlib as mpl
import matplotlib.colors as colors
import geopandas as gpd
from . import tools, core

def make_legend(fig, colormap, loc='best', title='',fontsize=12):
    """Make a figure legend wth provided color mapping"""

    import matplotlib.patches as mpatches
    pts=[]
    for c in colormap:
        pts.append(mpatches.Patch(color=colormap[c],label=c))

    fig.legend(handles=pts,loc=loc,fontsize=fontsize,
               frameon=False,draggable=True,title=title)
    return pts

def heatmap(df, cmap='gist_gray_r', w=15, h=5, ax=None):
    """Plot dataframe matrix"""

    if ax == None:
        fig, ax = plt.subplots(figsize=(w,h))
    im = ax.pcolor(df, cmap=cmap)
    ax.set_xticks(np.arange(len(df.columns))+0.5)
    ax.set_yticks(np.arange(len(df))+0.5)
    ax.set_xticklabels(df.columns)
    ax.set_yticklabels(df.index)
    plt.setp(ax.get_xticklabels(), rotation=80, ha="right",
             rotation_mode="anchor")
    plt.tight_layout()
    return

def random_colors(n=10, seed=1):
    """Generate random hex colors as list of length n."""

    import random
    random.seed(seed)
    clrs=[]
    for i in range(n):
        r = lambda: random.randint(0,255)
        c='#%02X%02X%02X' % (r(),r(),r())
        clrs.append(c)
    return clrs

def gen_colors(cmap,n,reverse=False):
    '''Generates n distinct color from a given colormap.
    Args:
        cmap(str): The name of the colormap you want to use.
            Refer https://matplotlib.org/stable/tutorials/colors/colormaps.html to choose
            Suggestions:
            For Metallicity in Astrophysics: Use coolwarm, bwr, seismic in reverse
            For distinct objects: Use gnuplot, brg, jet,turbo.
        n(int): Number of colors you want from the cmap you entered.
        reverse(bool): False by default. Set it to True if you want the cmap result to be reversed.
    Returns:
        colorlist(list): A list with hex values of colors.

    Taken from the mycolorpy package by binodbhttr
    see also https://matplotlib.org/stable/tutorials/colors/colormaps.html
    '''

    c_map = plt.cm.get_cmap(str(cmap)) # select the desired cmap
    arr=np.linspace(0,1,n) #create a list with numbers from 0 to 1 with n items
    colorlist=list()
    for c in arr:
        rgba=c_map(c) #select the rgba value of the cmap at point c which is a number between 0 to 1
        clr=colors.rgb2hex(rgba) #convert to hex
        colorlist.append(str(clr)) # create a list of these colors

    if reverse==True:
        colorlist.reverse()
    return colorlist

def show_colors(colors):
    """display a list of colors"""

    plt.figure(figsize=(6,1))
    plt.bar(range(len(colors)),height=1,color=colors,width=1)
    plt.axis('off')
    return

def get_color_mapping(df, col, cmap=None, seed=1):
    """Get random color map for categorical dataframe column"""

    c = df[col].unique()
    if cmap == None:
        rcolors = random_colors(len(c),seed)
    else:
        cmap = mpl.cm.get_cmap(cmap)
        rcolors = [cmap(i) for i in range(len(c))]
    colormap = dict(zip(c, rcolors))
    newcolors =  [colormap[i] if i in colormap else 'Black' for i in df[col]]
    return newcolors, colormap

def draw_pie(vals, xpos, ypos, colors, size=500, ax=None):
    """Draw a pie at a specific position on an mpl axis.
    Used to draw spatial pie charts on maps.
    Args:
        vals: values for pie
        xpos: x coord
        ypos: y coord
        colors: colors of values
        size: size of pie chart
    """

    cumsum = np.cumsum(vals)
    cumsum = cumsum / cumsum[-1]
    pie = [0] + cumsum.tolist()

    #colors = ["blue", "red", "yellow"]
    for i, (r1, r2) in enumerate(zip(pie[:-1], pie[1:])):
        angles = np.linspace(2 * np.pi * r1, 2 * np.pi * r2)
        x = [0] + np.cos(angles).tolist()
        y = [0] + np.sin(angles).tolist()
        xy = np.column_stack([x, y])
        ax.scatter([xpos], [ypos], marker=xy, s=size, color=colors[i], alpha=.9)

    return ax

def plot_top_category_in_grid(gdf, col, n_cells=10, cmap='viridis', ax=None):
    """
    Create a hexagonal grid, compute the most common category in each grid cell, and plot the result.

    Parameters:
    - gdf: GeoDataFrame, the input data containing points and the categorical column.
    - col: str, the name of the categorical column to analyze.
    - n_cells: int, the number of cells along one axis for the hexagonal grid.
    - cmap: str, the colormap to use for coloring categories.
    - ax: Matplotlib Axes, optional, axis to plot on.

    Returns:
    - grid: GeoDataFrame, the hexagonal grid with the top category in each cell.
    """

    grid = tools.create_hex_grid(gdf, n_cells=n_cells)
    merged = gpd.sjoin(gdf, grid, how='left', predicate='within')

    # Compute the top category per grid cell
    def aggtop(x):
        """Return the most common category in a group."""
        c = x.value_counts()
        if len(c) > 0:
            return c.index[0]
        return None

    dissolve = merged.dissolve(by="index_right", aggfunc={col: aggtop})
    clrs,cm = tools.get_color_mapping(gdf, col, cmap=cmap)
    grid['value'] = None
    grid['color'] = None
    grid.loc[dissolve.index, 'value'] = dissolve[col].values
    grid['color'] = grid['value'].map(cm)

    # Plot the grid with colors representing the top category
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))
    grid = grid[~grid['color'].isnull()]
    grid.plot(color=grid['color'], ec='black', alpha=0.8, lw=0.5, ax=ax)
    # Add a legend for categories
    legend_patches = [plt.Line2D([0], [0], marker='o', color=cm[cat],
                                  markersize=10, label=str(cat), lw=0)
                      for cat in cm]
    ax.legend(handles=legend_patches, title=col, fontsize='small', title_fontsize='medium')
    #plotting.make_legend(ax.figure,cm,title=col,loc=(1,.9),fontsize=8)
    ax.axis('off')
    return

def plot_grid_counts(grid, mask=None, cmap='plasma', fontsize=7):
    """Plot values from a pre-made grid, optionally mask the grid with a gdf"""

    if mask is not None:
        #mask for plotting
        mask_gdf = gpd.GeoDataFrame(geometry=[mask], crs='EPSG:29902')
        grid = gpd.overlay(grid, mask_gdf, how='intersection')
    fig, ax = plt.subplots(figsize=(9, 10))
    grid.plot(column='count', ax=ax, legend=True, cmap=cmap,
              edgecolor='k', linewidth=0.3, legend_kwds={"shrink":.5})
    ax.axis('off')
    ax.set_title('gaps in coverage')
    for _, row in grid.iterrows():
        centroid = row.geometry.centroid  # Get the centroid of the hexagon
        count = row['count']
        ax.text(centroid.x, centroid.y, str(count),
                ha='center', va='center', fontsize=fontsize, color='white')
    return