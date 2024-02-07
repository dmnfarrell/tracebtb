#!/usr/bin/env python

"""
    btbwgstool plotting functions.
    Created Jan 2021
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

import sys,os,subprocess,glob,shutil,re,random,time
import platform
import numpy as np
import pandas as pd
import pylab as plt
import matplotlib as mpl
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Phylo, AlignIO, Align
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon, MultiPolygon

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """

    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)

def get_cmd(cmd):
    """Get windows version of a command if required"""

    if getattr(sys, 'frozen', False):
        cmd = resource_path('bin\%s.exe' %cmd)
    elif platform.system() == 'Windows':
        cmd = os.path.join(bin_path, '%s.exe' %cmd)
    return cmd

def random_hex_color():
    """random hex color"""

    r = lambda: np.random.randint(0,255)
    c='#%02X%02X%02X' % (r(),r(),r())
    return c

def random_hex_colors(n=1,seed=None):
    if seed != None:
        np.random.seed(seed)
    return [random_hex_color() for i in range(n)]

def colormap_colors(colormap_name, n):
    """Colors list from mpl colormap"""

    colormap = plt.cm.get_cmap(colormap_name, n)
    colors = [mpl.colors.rgb2hex(colormap(i)) for i in range(n)]
    return colors

def colormap_from_labels(colormap_name, labels):
    """Get dict of colors mapping to labels using mpl colormap"""

    n = len(labels)
    colormap = plt.cm.get_cmap(colormap_name, n)
    colors = {labels[i]: mpl.colors.rgb2hex(colormap(i)) for i in range(n)}
    return colors

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

def random_grayscale_color(_):

    gray_value = np.random.randint(0, 256)
    return f'#{gray_value:02X}{gray_value:02X}{gray_value:02X}'

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

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    c_map = plt.cm.get_cmap(str(cmap)) # select the desired cmap
    arr = np.linspace(0,1,n) #create a list with numbers from 0 to 1 with n items
    colorlist=list()
    for c in arr:
        rgba=c_map(c) #select the rgba value of the cmap at point c which is a number between 0 to 1
        clr=colors.rgb2hex(rgba) #convert to hex
        colorlist.append(str(clr))

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

    import matplotlib.colors as colors
    c = df[col].unique()
    if cmap == None:
        clrs = random_colors(len(c),seed)
    else:
        c_map = mpl.cm.get_cmap(cmap)
        clrs = [colors.rgb2hex(c_map(i)) for i in range(len(c))]
        #colors = gen_colors(cmap,len(c))
        #print (clrs)

    colormap = dict(zip(c, clrs))
    newcolors =  [colormap[i] if i in colormap else 'Black' for i in df[col]]
    return newcolors, colormap

def alignment_from_snps(df):
    """MultipleSeqAlignment object from core snps"""

    df = df.set_index('pos').T
    seqs=[]
    for i,r in df.iterrows():
        s = ''.join(list(r))
        seqs.append(SeqRecord(Seq(s),id=i))

    aln = Align.MultipleSeqAlignment(seqs)
    return aln

def snp_dist_matrix(aln):
    """
    Compute the number of Single Nucleotide Polymorphisms (SNPs)
    between sequences in a Biopython alignment.
    Args:
        aln:
            Biopython multiple sequence alignment object.
    returns:
        a matrix as pandas dataframe
    """

    names=[s.id for s in aln]
    num_sequences = len(aln)
    matrix = np.zeros((num_sequences, num_sequences))

    for i in range(num_sequences):
        seq1 = str(aln[i].seq)
        for j in range(i + 1, num_sequences):
            seq2 = str(aln[j].seq)
            # Calculate the number of SNPs
            snp_count = sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
            matrix[i, j] = snp_count
            matrix[j, i] = snp_count

    m = pd.DataFrame(matrix,index=names,columns=names).astype(int)
    return m

def dist_matrix_to_mst(distance_matrix, df=None, colorcol=None, cmap='Set1', ax=None):
    """Dist matrix to minimum spanning tree"""

    if ax == None:
        fig,ax=plt.subplots()
    import networkx as nx
    from networkx.drawing.nx_agraph import graphviz_layout
    import pylab as plt
    G = nx.Graph()

    for i, row in distance_matrix.iterrows():
        for j, weight in row.items():
            G.add_edge(i, j, weight=weight)

    T = nx.minimum_spanning_tree(G, algorithm='kruskal')
    # Compute edge lengths based on distances
    edge_lengths = [T[u][v]['weight'] for u, v in T.edges()]
    # Plot the minimum spanning tree with edge lengths proportional to distances
    #pos = nx.spring_layout(T)#, weight='weight', scale=10, seed=42)
    pos = graphviz_layout(T)
    labels = nx.get_edge_attributes(T, 'weight')
    if df is not None:
        c,cmap = get_color_mapping(df, colorcol, cmap)
        colors = dict(zip(df.index, c))
        node_colors = [colors[node] for node in T.nodes()]
        p=make_legend(ax.figure, cmap, loc=(1, .9), title=colorcol,fontsize=9)
    else:
        node_colors = 'lightblue'
    nx.draw_networkx(T, pos, node_color=node_colors,node_size=40,font_size=6, ax=ax)
    nx.draw_networkx_edge_labels(T, pos, edge_labels=labels, font_size=6, ax=ax)
    #nx.draw_networkx_edges(T, pos, width=edge_lengths)
    ax.axis('off')
    return

def make_legend(fig, colormap, loc=(1.05, .6), title='',fontsize=12):
    """Make a figure legend wth provided color mapping"""

    import matplotlib.patches as mpatches
    pts=[]
    for c in colormap:
        pts.append(mpatches.Patch(color=colormap[c],label=c))
    fig.legend(handles=pts,bbox_to_anchor=loc,fontsize=fontsize,title=title)
    return pts

def subset_alignment(aln, names):
    """Subset of a multpleseqalignment object"""

    seqs = [rec for rec in aln if rec.id in names]
    new = Align.MultipleSeqAlignment(seqs)
    return new

def get_coords_data(df):
    """Get coordinates from geodataframe as linestrings"""

    df['P2'] = df.geometry.shift(-1)
    coords = df[:-1].apply(lambda x: LineString([x.geometry,x.P2]),1)
    return coords

def get_dataframe_memory(df):
    """get dataframe memory in MB"""

    m = df.memory_usage(deep=True).sum()
    m = round(m/1048576,2)
    return m
