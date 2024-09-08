#!/usr/bin/env python

"""
    Utility functions for tracebtb.
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

import sys,os,subprocess,glob,shutil,re
import random,time
import json
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

module_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(module_path,'data')

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

def random_colormap_colors(cmap_name, n):
    """
    Get n random colors from a Matplotlib colormap.
    """
    clrs = colormap_colors(cmap_name, n)
    return random.sample(clrs,len(clrs))

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

def compute_snp_count(args):
    i, j, seq1, seq2 = args
    return i, j, np.sum(np.fromiter((c1 != c2 for c1, c2 in zip(seq1, seq2)), dtype=int))

def snp_dist_matrix(aln, threads=4):
    """
    Compute the number of Single Nucleotide Polymorphisms (SNPs)
    between sequences in a Biopython alignment.
    Args:
        aln:
            Biopython multiple sequence alignment object.
    Returns:
        A matrix as pandas dataframe.
    """

    from multiprocessing import Pool
    names = [s.id for s in aln]
    num_sequences = len(aln)
    matrix = np.zeros((num_sequences, num_sequences), dtype=int)
    sequences = [str(s.seq) for s in aln]

    args_list = []
    for i in range(num_sequences):
        seq1 = sequences[i]
        for j in range(i + 1, num_sequences):
            seq2 = sequences[j]
            args_list.append((i, j, seq1, seq2))

    with Pool(processes=threads) as pool:
        results = pool.map(compute_snp_count, args_list)

    for i, j, snp_count in results:
        matrix[i, j] = snp_count
        matrix[j, i] = snp_count

    m = pd.DataFrame(matrix, index=names, columns=names)
    return m

def dist_matrix_to_mst(distance_matrix, df=None, colorcol=None, labelcol=None, colormap=None,
                       cmap='Set1', node_size=4, font_size=8, with_labels=False,
                       edge_labels=False, legend_loc=(1, .7), ax=None):
    """Dist matrix to minimum spanning tree"""

    from . import plotting
    import pylab as plt
    if ax == None:
        fig,ax=plt.subplots()
    import networkx as nx
    from networkx.drawing.nx_agraph import graphviz_layout

    G = nx.Graph()

    for i, row in distance_matrix.iterrows():
        for j, weight in row.items():
            G.add_edge(i, j, weight=weight)

    T = nx.minimum_spanning_tree(G, algorithm='kruskal')
    # Compute edge lengths based on distances
    edge_lengths = [T[u][v]['weight'] for u, v in T.edges()]
    # Plot the minimum spanning tree with edge lengths proportional to distances
    pos = graphviz_layout(T)
    labels = nx.get_edge_attributes(T, 'weight')
    if df is not None:
        l = [label for label in T.nodes if label in df.index]
        df = df.loc[l]
        if colormap is None:
            colors,cmap = plotting.get_color_mapping(df, colorcol, cmap)
        else:
            #custom colormap if provided
            colors = [colormap[i] if i in colormap else 'black' for i in df[colorcol]]
            cmap = colormap
        #print (cmap)
        C = dict(zip(df.index, colors))
        node_colors = [C[node] if node in C else 'Black' for node in T.nodes()]
        #checks that colormap matches nodes so legend doesn't have crap in it
        cmap = check_keys(cmap, df[colorcol].unique())
        #add legend for node colors
        p = plotting.make_legend(ax.figure, cmap, loc=legend_loc, title=colorcol,fontsize=9)

    else:
        node_colors = 'black'
    nx.draw_networkx(T, pos, node_color=node_colors, node_size=node_size,
                     font_size=font_size, with_labels=with_labels, ax=ax)
    if edge_labels == True:
        nx.draw_networkx_edge_labels(T, pos, edge_labels=labels, font_size=font_size*.8, ax=ax)

    if labelcol not in [None,'']:
        node_labels = {node:df.loc[node][labelcol] if node in df.index else '' for node in T.nodes()}
        #print (node_labels)
        nx.draw_networkx_labels(T, pos, labels=node_labels, font_size=font_size,
                 horizontalalignment='right',verticalalignment='top')
    ax.axis('off')
    return T, pos

def get_within_distance(dm, row_index, threshold):
    """
    Get the closest samples for a given row in a distance matrix DataFrame.
    """

    row = dm.loc[row_index]
    x = row.loc[row<=threshold].index
    return x

def check_keys(d, vals):
    """Remove keys not in vals"""

    keys=list(d.keys())
    for key in keys:
        if not key in vals:
            d.pop(key, None)
    return d

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

    if df is None:
        return
    m = df.memory_usage(deep=True).sum()
    m = round(m/1048576,2)
    return m

def df_html(df, fontsize='8pt'):
    """Create df html enclosed in some css"""

    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>Styled DataFrame</title>
        <style>
            body {{
                font-family: monospace;
                font-size: {fontsize};
                margin: 0;
                padding: 0;
            }}
            table {{
                border-collapse: collapse;
                width: 100%;
            }}
            th, td {{
                border: 1px solid #ddd;
                padding: 8px;
            }}
            th {{
                background-color: white;
                text-align: left;
                position: sticky;
                top: 0;
            }}

        </style>
    </head>
    <body>
        {df.to_html()}
    </body>
    </html>
    """
    return html

def nearest(point, gdf):
    """Get nearest neighbours to point in gdf, vector based"""

    from shapely.ops import nearest_points
    dists = gdf.apply(lambda row:  point.distance(row.geometry),axis=1)
    #print (dists[dists>0].idxmin())
    return dists[dists>0].min()

def find_outliers(cent, min_dist=10, min_samples=10, col='snp7'):
    """
    Find outlier cluster points
    """

    res=[]
    for i,df in cent.groupby(col):
        if len(df)>=min_samples:
            df['nearest'] = df.geometry.apply(lambda x: nearest(x,df))
            df['outlier'] = df.nearest>min_dist*100
            #print (i,len(df))
        res.append(df)
    res = pd.concat(res)
    outliers = res[res.outlier==True]
    return outliers

def get_ordinal_columns(df):
    """Try to get ordinal columns from table"""

    ignore = ['sample','Aliquot','geometry','Animal_ID','X_COORD','Y_COORD']
    ocols = []
    for col in df.columns:
        count = df[col].nunique()
        #print (col, count, len(df[col]))
        if count==len(df[col]) or col in ignore:
            continue
        ocols.append(col)
    return ocols

def herd_summary(df, moves, snpdist=None):
    """Herd summary"""

    res=[]
    for herd, sdf in df.groupby('HERD_NO'):
        #print (herd)
        clades = len(sdf.snp7.unique())
        m = moves[moves.move_to == herd]
        #only farm to farm moves
        m = m[m.data_type=='F_to_F']
        #get mean SNP dist within farm
        idx = list(sdf.index)
        if snpdist is not None:
            D = snpdist.loc[idx,idx]
            meandist = D.stack().mean().round(1)
            mediandist = D.stack().median().round(1)
        else:
            meandist = None
            mediandist = None
        if 'Homebred' in sdf.columns:
            hbred = len(sdf[sdf.Homebred=='yes'])
        else:
            hbred = None
        if 'County' in sdf.columns:
            cty = ';'.join(list(sdf.dropna(subset='County').County.unique()))
        if 'snp7' in sdf.columns:
            cl = ';'.join(list(sdf.snp7.unique()))
        res.append([herd,len(sdf),clades,len(m),mediandist,hbred,cl,cty])
    res = pd.DataFrame(res, columns=['HERD_NO','isolates','strains','moves',
                                     'median_dist','homebred','snp7_cl','County'])
    res = res.sort_values('strains',ascending=False)
    return res

def cluster_summary(df, col, min_size=5, snpdist=None):
    """Group summary e.g. by cluster"""

    res = []
    for c,sub in df.groupby(col):
        if len(sub)<min_size:
            continue
        herds = len(sub.HERD_NO.unique())
        idx = list(sub.index)
        if snpdist is not None:
            D = snpdist.loc[idx,idx]
            meandist = D.stack().mean().round(1)
            mediandist = D.stack().median().round(1)
        if 'Homebred' in sub.columns:
            hbred = len(sub[sub.Homebred=='yes'])
        else:
            hbred = None
        #species=sub.Species.value_counts()
        badger = len(sub[sub.Species=='Badger'])
        clade = sub.iloc[0].IE_clade
        row = [c,len(sub),hbred,herds,badger,mediandist,clade]
        res.append(row)
    res=pd.DataFrame(res, columns=['cluster','isolates','homebred','herds','badger','median_dist','IE_clade'])
    res=res.sort_values('isolates',ascending=False)
    return res

def get_moves_bytag(df, move_df, lpis_cent):
    """
    Get moves and coords for one or more samples.
    """

    cols=['Animal_ID']+list(move_df.columns)
    t = df.merge(move_df,left_on='Animal_ID',right_on='tag',how='inner')[cols]
    #print (t)
    #merge result with parcel cents to get coords of moved_to farms
    m = t.merge(lpis_cent,left_on='move_to',right_on='SPH_HERD_N', how='left')
    if len(m)==0:
        return
    m = (m.drop_duplicates()
        .drop(columns=['Animal_ID','id','event_type','rnk_final'], errors='ignore')
        .sort_values(['tag','move_date'])
        .set_index('tag',drop=True)
        )
    m = gpd.GeoDataFrame(m)
    return m

def get_largest_poly(x):
    """Return largest poly of multipolygon"""

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

def get_county(x):

    counties = gpd.read_file(os.path.join(data_path,'counties.shp')).to_crs("EPSG:3857")
    if x.geometry.is_empty:
        return 'NA'
    found = counties[counties.contains(x.geometry)]
    if len(found)>0:
        return found.iloc[0].NAME_TAG

def get_year(x, key='sample'):

    if not str(x[key]).startswith('TB') or len(x[key])>11:
        return
    yr = x[key].split('-')[0][2:]
    return int('20'+yr)

def add_to_json(filename, obj, key):
    """
    Append an object to existing dict in json file.
    Used for saving selections for example. Adds object to root dict.
    """

    if os.path.exists(filename):
        data = json.load(open(filename,'r'))
    else:
        data = {}
    data[key] = obj
    #save file
    with open(filename,'w') as f:
        f.write(json.dumps(data))
    return