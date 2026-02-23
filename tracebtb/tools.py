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
import pickle
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
from . import core

module_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(module_path,'data')

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """

    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)

def save_project(filename, data):
    """Save an initial project"""

    keys = ['meta','moves','parcels','lpis_cent','aln','snpdist',
            'lpis_cent','lpis_master_file','testing']
    for k in keys:
        if k not in data:
            print (f'warning. {k} not found in dict')
    pickle.dump(data, open(filename,'wb'))
    return

def update_project(filename, new, field, save=True):
    """Update tracebtb project file with data"""

    import pickle
    data = pickle.load(open(filename,'rb'))
    print (data.keys())
    if field in data:
        print ('overwriting field %s' %field)
    data[field] = new
    if save == True:
        pickle.dump(data, open(filename,'wb'))
    return

def remove_from_project(filename, field):
    """Remove a key from the project dict"""

    import pickle
    data = pickle.load(open(filename,'rb'))
    if field in data:
        del data[field]
    pickle.dump(data, open(filename,'wb'))
    return

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

def get_color_mapping(df, col, cmap=None, seed=12):
    """
    Get colors for a dataframe column using a mapping.
    If type is object or categorical, is cmap is None we use random colors.
    If numeric column use viridis map and add colorbar.
    """

    import matplotlib.colors as mcolors
    # 1. Identify unique values
    unique_vals = np.sort(df[col].unique())
    n = len(unique_vals)
    if cmap == None:
        clrs = random_colors(n,seed)
    else:
        c_map = mpl.colormaps.get(cmap)
        if pd.api.types.is_numeric_dtype(df[col]):
            # NUMERIC: Normalize index to [0, 1] for a smooth gradient
            if n > 1:
                clrs = [mcolors.rgb2hex(c_map(i/(n-1))) for i in range(n)]
            else:
                clrs = [mcolors.rgb2hex(c_map(0.5))]
        else:
            clrs = [mcolors.rgb2hex(c_map(i)) for i in range(n)]

    colormap = dict(zip(unique_vals, clrs))
    values = df[col].fillna(0)
    newcolors = [colormap[val] for val in values]
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
                       edge_labels=False, legend=True, legend_loc=(1, .7), ax=None):
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
        if legend == True:
            p = plotting.make_legend(ax.figure, cmap, loc=legend_loc, title=colorcol, fontsize=12)

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

def pca(matrix, n_components=3):
    """Perform PCA"""

    import sklearn
    from sklearn import decomposition
    from sklearn import manifold

    pca = decomposition.PCA(n_components=n_components)
    mds = manifold.MDS(n_components=n_components)
    pos = mds.fit(matrix).embedding_
    X = pca.fit_transform(pos)
    #X = pca.transform(C)
    idx = list(matrix.index)
    df = pd.DataFrame(X,index=idx)
    #df.columns = matrix.columns
    return df

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

def create_grid(gdf=None, bounds=None, n_cells=10, overlap=False, crs="EPSG:29902"):
    """Create square grid that covers a geodataframe area
    or a fixed boundary with x-y coords
    returns: a GeoDataFrame of grid polygons
    """

    import shapely

    if bounds != None:
        xmin, ymin, xmax, ymax= bounds
    else:
        xmin, ymin, xmax, ymax= gdf.total_bounds

    # get cell size
    cell_size = (xmax-xmin)/n_cells
    # create the cells in a loop
    grid_cells = []
    for x0 in np.arange(xmin, xmax+cell_size, cell_size ):
        for y0 in np.arange(ymin, ymax+cell_size, cell_size):
            x1 = x0-cell_size
            y1 = y0+cell_size
            poly = shapely.geometry.box(x0, y0, x1, y1)
            #print (gdf.overlay(poly, how='intersection'))
            grid_cells.append( poly )

    cells = gpd.GeoDataFrame(grid_cells, columns=['geometry'],
                                     crs=crs)
    if overlap == True:
        cols = ['grid_id','geometry','grid_area']
        cells = cells.sjoin(gdf, how='inner').drop_duplicates('geometry')
    return cells

def create_hex_grid(gdf=None, bounds=None, n_cells=10, overlap=False, crs="EPSG:29902"):
    """Hexagonal grid over geometry.
    See https://sabrinadchan.github.io/data-blog/building-a-hexagonal-cartogram.html
    """

    from shapely.geometry import Polygon
    if bounds != None:
        xmin, ymin, xmax, ymax= bounds
    else:
        xmin, ymin, xmax, ymax= gdf.total_bounds

    unit = (xmax-xmin)/n_cells
    a = np.sin(np.pi / 3)
    cols = np.arange(np.floor(xmin), np.ceil(xmax), 3 * unit)
    rows = np.arange(np.floor(ymin) / a, np.ceil(ymax) / a, unit)

    #print (len(cols))
    hexagons = []
    for x in cols:
      for i, y in enumerate(rows):
        if (i % 2 == 0):
          x0 = x
        else:
          x0 = x + 1.5 * unit

        hexagons.append(Polygon([
          (x0, y * a),
          (x0 + unit, y * a),
          (x0 + (1.5 * unit), (y + unit) * a),
          (x0 + unit, (y + (2 * unit)) * a),
          (x0, (y + (2 * unit)) * a),
          (x0 - (0.5 * unit), (y + unit) * a),
        ]))

    grid = gpd.GeoDataFrame({'geometry': hexagons},crs=crs)
    #grid["grid_area"] = grid.area
    grid = grid.reset_index().rename(columns={"index": "grid_id"})
    if overlap == True:
        cols = ['grid_id','geometry']#,'grid_area']
        grid = grid.sjoin(gdf, how='inner').drop_duplicates('geometry')[cols]
    return grid

def count_points_in_grid(grid, gdf, value_column='count', threshold=0):
    """
    Count the number of points in a gdf in each grid cell.

    Parameters:
    - grid: GeoDataFrame, the hexagonal grid created by `create_hex_grid`.
    - gdf: GeoDataFrame, the GeoDataFrame of points.
    - value_column: str, the name of the column to store point counts in the grid GeoDataFrame.
    - threshold: only return cells with counts above this value, None to ignore

    Returns:
    - grid: GeoDataFrame with a new column `value_column` containing the count of points.
    """

    import matplotlib.pyplot as plt

    # Perform a spatial join to find points within each hexagon
    joined = gpd.sjoin(gdf, grid, how='inner', predicate='within')
    # Count the number of points in each grid cell
    point_counts = joined.groupby('grid_id').size()
    # Add point counts to the grid GeoDataFrame
    grid[value_column] = grid['grid_id'].map(point_counts).fillna(0).astype(int)
    if threshold != None:
        grid = grid[grid['count']>threshold].copy()
    return grid

def get_count_grid(df, grid=None, n_cells=30):
    """Get hex grid of counts for a gdf of points"""

    if grid is None:
        grid = create_hex_grid(df, n_cells=n_cells)
    grid = count_points_in_grid(grid, df)
    #mask grid to map
    mask = core.counties_gdf.to_crs('EPSG:29902').union_all()
    grid = grid[grid.intersects(mask)].copy()
    return grid, mask

def nearest(point, gdf):
    """Get nearest neighbours to point in gdf, vector based"""

    from shapely.ops import nearest_points
    dists = gdf.apply(lambda row: point.distance(row.geometry),axis=1)
    return dists[dists>0].min()

def find_nearest_point(multipolygon, points_gdf):

    try:
        centroid = multipolygon.centroid
        distances = points_gdf.geometry.distance(centroid)
        nearest_point_index = distances.idxmin()
        nearest_point = points_gdf.loc[nearest_point_index]
    except:
        return
    return nearest_point

def find_points_within_distance(gdf, target_point, distance):
    """
    Find points in a GeoDataFrame within a specified distance from a target point.

    Parameters:
        gdf (GeoDataFrame): A GeoDataFrame containing point geometries.
        target_point (Point): A Shapely Point representing the target location.
        distance (float): The maximum distance to search.

    Returns:
        GeoDataFrame: A GeoDataFrame of points within the specified distance.
    """
    # Ensure the GeoDataFrame has a valid geometry column
    if 'geometry' not in gdf or gdf.geometry.is_empty.any():
        raise ValueError("GeoDataFrame must contain a valid geometry column with point data.")
    # Ensure the GeoDataFrame and target point are in the same CRS
    if gdf.crs is None:
        raise ValueError("GeoDataFrame must have a defined CRS.")
    # Buffer the target point by the distance and check for intersection
    buffer = target_point.buffer(distance)
    points_within_distance = gdf[gdf.geometry.intersects(buffer)]
    return points_within_distance

def find_outliers(gdf, min_dist=10, min_samples=10, col='snp7'):
    """
    Find outlier cluster points
    """

    res=[]
    for i,df in gdf.groupby(col):
        if len(df)>=min_samples:
            df['nearest'] = df.geometry.apply(lambda x: nearest(x,df))
            df['outlier'] = df.nearest>min_dist*100
            #print (i,len(df))
        res.append(df)
    res = pd.concat(res)
    outliers = res[res.outlier==True]
    return outliers

def remove_outliers_zscore(gdf, threshold=3):
    """
    Remove outliers from a GeoDataFrame using distance from the centroid.
    :param gdf: GeoDataFrame with a 'geometry' column containing points.
    :param threshold: Z-score threshold for identifying outliers based on distances.
    :return: GeoDataFrame with outliers removed.
    """
    import numpy as np
    from scipy.stats import zscore

    # Extract x and y coordinates
    gdf['x'] = gdf.geometry.x
    gdf['y'] = gdf.geometry.y
    gdf = gdf.dropna(subset=['x', 'y']).copy()
    # Calculate centroid of the points
    centroid = gdf[['x', 'y']].mean()
    # Calculate Euclidean distance of each point from the centroid
    gdf['distance'] = np.sqrt((gdf['x'] - centroid['x'])**2 + (gdf['y'] - centroid['y'])**2)
    # Calculate Z-scores for distances
    gdf['distance_zscore'] = zscore(gdf['distance'])
    # Filter based on the threshold
    gdf_filtered = gdf[np.abs(gdf['distance_zscore']) < threshold]
    # Drop the helper columns
    gdf_filtered = gdf_filtered.drop(columns=['x', 'y', 'distance', 'distance_zscore'])
    #print (len(gdf),len(gdf_filtered))
    return gdf_filtered

def remove_outliers_dbscan(gdf, eps=0.1, min_samples=5):
    """
    Remove outliers from a GeoDataFrame using DBSCAN clustering.
    :param gdf: GeoDataFrame with a 'geometry' column containing points.
    :param eps: Maximum distance between two samples for them to be considered as in the same neighborhood.
    :param min_samples: Minimum number of points to form a cluster.
    :return: GeoDataFrame with outliers removed.
    """

    from sklearn.cluster import DBSCAN
    coords = np.array([(point.x, point.y) for point in gdf.geometry])
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(coords)
    gdf['cluster'] = db.labels_  # Assign cluster labels
    # Keep only points in clusters (label >= 0)
    gdf_filtered = gdf[gdf['cluster'] >= 0].drop(columns=['cluster'])
    return gdf_filtered

def remove_outliers_mahalanobis(gdf, threshold=3):
    """
    Remove outliers using the Mahalanobis distance.
    :param gdf: GeoDataFrame with a 'geometry' column containing points.
    :param threshold: Mahalanobis distance threshold for identifying outliers.
    :return: GeoDataFrame with outliers removed.
    """

    from scipy.spatial.distance import mahalanobis
    coords = np.array([(point.x, point.y) for point in gdf.geometry])
    # Compute the mean and covariance matrix
    mean = np.mean(coords, axis=0)
    cov_matrix = np.cov(coords, rowvar=False)
    inv_cov_matrix = np.linalg.inv(cov_matrix)
    # Compute Mahalanobis distance for each point
    distances = [mahalanobis(point, mean, inv_cov_matrix) for point in coords]
    # Filter points within the threshold
    gdf_filtered = gdf[np.array(distances) <= threshold].copy()
    return gdf_filtered

def remove_outliers_convex_hull(gdf, threshold=1.5):
    """
    Removes outliers based on distance from the convex hull.

    Parameters:
    gdf (GeoDataFrame): A GeoDataFrame containing spatial points.
    threshold (float): The distance threshold for considering a point as an outlier.

    Returns:
    GeoSeries: the hull polygon
    GeoDataFrame: Filtered GeoDataFrame without outliers.
    """
    from scipy.spatial import ConvexHull
    if gdf.empty or len(gdf) < 3:
        print("Not enough points to form a convex hull.")
        return gdf

    coords = np.array([[p.x, p.y] for p in gdf.geometry])
    # Compute convex hull
    hull = ConvexHull(coords)
    hull_polygon = gpd.GeoSeries([Point(coords[i]) for i in hull.vertices]).union_all().convex_hull
    # Compute distance of each point from the convex hull
    gdf["hull_distance"] = gdf.geometry.apply(lambda p: p.distance(hull_polygon))
    # Identify outliers using a threshold (e.g., 1.5 * IQR rule)
    q1, q3 = np.percentile(gdf["hull_distance"], [25, 75])
    iqr = q3 - q1
    cutoff = iqr * threshold
    upper_bound = q3 + cutoff
    # Filter out the outliers
    filtered_gdf = gdf[gdf["hull_distance"] <= upper_bound].drop(columns=["hull_distance"])
    hull_df = gpd.GeoDataFrame(geometry=[hull_polygon])
    return hull_df, filtered_gdf

def get_area(gdf):
    """Compute the convex hull area of the points"""

    from scipy.spatial import ConvexHull
    coords = np.array([(point.x, point.y) for point in gdf.drop_duplicates('geometry').geometry])
    if len(coords) < 3:
        return 0
    hull = ConvexHull(coords)
    hull_polygon = gpd.GeoSeries([gpd.points_from_xy(coords[hull.vertices, 0], coords[hull.vertices, 1]).union_all().convex_hull])
    # Return the area of the convex hull
    return hull_polygon.area.values[0]/1000

def flatten_matrix(df):
    """Flatten a symmetrical matrix"""

    #user only upper triangle
    keep = np.triu(np.ones(df.shape)).astype('bool').reshape(df.size)
    S=df.unstack()[keep]
    S.index = ['{}_{}'.format(i, j) for i, j in S.index]
    S = S[~S.index.duplicated(keep='first')]
    return S

def get_ordinal_columns(df):
    """Try to get ordinal columns from table"""

    ignore = ['sample','Aliquot','geometry','Animal_ID','X_COORD','Y_COORD','LONG','LAT']
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
        clades = len(sdf.short_name.unique())
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
        if 'short_name' in sdf.columns:
            cl = ';'.join(list(sdf.short_name.astype(str).unique()))
        res.append([herd,len(sdf),clades,len(m),mediandist,hbred,cl,cty,sdf.iloc[0].class_2021])
    res = pd.DataFrame(res, columns=['HERD_NO','isolates','strains','moves',
                                     'median_dist','homebred','strain_names','County','class_2021'])
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
        clade = sub.iloc[0].lineage
        row = [c,len(sub),hbred,herds,badger,mediandist,clade]
        res.append(row)
    res=pd.DataFrame(res, columns=['group','isolates','homebred','herds','badger','median_dist','lineage'])
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

def get_moves_locations(move_df, lpis_cent):
    """Get location for each move row in movement df"""

    m = move_df.merge(lpis_cent,left_on='herd',right_on='SPH_HERD_N', how='left')
    if len(m)==0:
        return
    m = (m.drop_duplicates()
        #.drop(columns=['Animal_ID','id','event_type','rnk_final'], errors='ignore')
        .sort_values(['tag','event_date'])
        .set_index('tag',drop=True)
        )
    return gpd.GeoDataFrame(m)

def get_last_move(df, moves):
    """Add last move to animal level data"""

    last_moves = moves.loc[moves.groupby('tag')['move_date'].idxmax()].reset_index(drop=True)[['tag','move_date','data_type']]
    last_moves.columns=['tag','last_move','last_move_type']
    df = df.merge(last_moves,left_on='Animal_ID',right_on='tag',how='left').set_index(df.index).drop(columns='tag')
    return df

def apply_jitter(gdf, radius=100):
    """
    Apply jitter to the points in a GeoDataFrame.

    Parameters:
        gdf (GeoDataFrame): A GeoDataFrame with Point geometries.
        jitter_amount (float): The amount of jitter to apply (default: 0.001 degrees).

    Returns:
        GeoDataFrame: A new GeoDataFrame with jittered points.
    """
    def jitter_point(point):
        return Point(
            point.x + np.random.uniform(-radius, radius),
            point.y + np.random.uniform(-radius, radius)
        )
    gdf_jittered = gdf.copy()
    gdf_jittered["geometry"] = gdf_jittered["geometry"].apply(jitter_point)
    return gdf_jittered

def jitter_by_farm(gdf, radius=100):
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

def mean_distance_between_sets(distance_matrix, indices1, indices2):
    """
    Calculates the mean distance between two sets of rows in a distance matrix.

    Parameters:
        distance_matrix (numpy.ndarray): A 2D square matrix of distances.
        set1_indices (list of int): Indices representing the first set of rows.
        set2_indices (list of int): Indices representing the second set of rows.

    Returns:
        float: The mean distance between the two sets of rows.
    """
    # Extract the submatrix of distances between the two sets
    submatrix = distance_matrix.loc[indices1, indices2]
    # Compute the mean of all elements in the submatrix
    mean_distance = np.mean(submatrix)
    return mean_distance

def mean_geo_dist(gdf1, gdf2):
    """Compute the distance between 2 GeoDataFrames"""

    centroid1 = gdf1.geometry.centroid.union_all().centroid
    centroid2 = gdf2.geometry.centroid.union_all().centroid
    distance = centroid1.distance(centroid2)
    return distance

def compare_cluster_distances(df,snpdist,col,min_size=5):
    """Compare mean cluster vs geo distances"""

    grouped = df.groupby(col)
    filtered_groups = {name: group for name, group in grouped if len(group) >= min_size}
    unique_groups = list(filtered_groups.keys())
    from itertools import combinations
    res=[]
    for group1, group2 in combinations(unique_groups, 2):
        # Get the data for the two groups
        sub1 = filtered_groups[group1]
        sub2 = filtered_groups[group2]
        #if len(sub1.HERD_NO.unique())<=1: continue
        #if len(sub2.HERD_NO.unique())<=1: continue
        sd = mean_distance_between_sets(snpdist, sub1.index, sub2.index)
        gd = mean_geo_dist(sub1, sub2)
        res.append((sd,gd))
    res=pd.DataFrame(res,columns=['snpdist','geodist'])
    return res

def spatial_cluster(gdf, eps=3000):
    """Spatially cluster points"""

    import sklearn.cluster as skc
    coordinates = gdf['geometry'].apply(lambda p: np.hstack(p.xy)).values
    coordinates = np.vstack(coordinates)
    clusterer = skc.DBSCAN(eps=eps).fit(coordinates)
    nclusters = clusterer.p
    gdf = gdf.assign(cl=clusterer.labels_)
    return gdf

def find_neighbouring_points(gdf, sub, dist, key='HERD_NO'):
    """Find nearby points from one dataframe to another.
    use key to exclude points belonging to source gdf e.g. HERD_NO
    """

    point = sub.iloc[0].geometry
    p = find_points_within_distance(gdf, point, dist)
    p = p[~p[key].isin(list(sub[key]))]
    return p

def find_neighbours(gdf, dist, lpis_cent, lpis):
    """
    Find neighbouring herds from lpis given a set of points.
    Returns: herd land parcels
    """

    found = []
    for x in gdf.geometry:
        dists = lpis_cent.distance(x)
        points = lpis_cent[(dists<=dist) & (dists>10)]
        found.append(points)
    if len(found) == 0:
        return lpis_cent.drop(lpis_cent.index)
    found = pd.concat(found).drop_duplicates()
    x = lpis[lpis.SPH_HERD_N.isin(found.SPH_HERD_N)]
    #exclude those in source gdf
    x = x[~x.SPH_HERD_N.isin(gdf.HERD_NO)]
    return x

def shared_borders(parcels, lpis):
    """
    Find neighbouring herds with shared borders using land parcels.
    """

    found = []
    if type(parcels) is pd.Series:
        parcels = parcels.to_frame().T
    for i, r in parcels.iterrows():
        #farms with shared borders
        polygon = r.geometry
        x = lpis[lpis.touches(polygon)]
        found.append(x)
    if len(found) == 0:
        return
    found = pd.concat(found)
    found = found[~found.SPH_HERD_N.isin(parcels.SPH_HERD_N)]
    return found

def plot_neighbours(pcl,df,col=None,pad=.3,ax=None):

    if ax==None:
        fig,ax=plt.subplots()
    point = pcl.iloc[0].geometry.centroid
    df.plot(lw=.5,ec='black',alpha=0.6,column=col,legend=False,cmap='Set3',ax=ax)
    pcl.plot(color='red',lw=0,ec='black',ax=ax)
    x1,y1,x2,y2 = df.union_all().bounds
    pad = (x2-x1)*pad
    ax.set_xlim(point.x-pad,point.x+pad)
    ax.set_ylim(point.y-pad,point.y+pad)
    ax.axis('off')
    return

def herd_parcel_metrics(herd, lpis, lpis_cent, dist=1000):
    """Get metrics for a herd in relation to neighbours"""

    pcl = lpis[lpis.SPH_HERD_N==herd].copy()
    pcl['HERD_NO']=pcl.SPH_HERD_N
    #get all farms within radius or contiguous parcels
    cont_parcels = shared_borders(pcl, lpis)
    #nb = tools.shared_borders(row, lpis)
    nb = find_neighbours(pcl, dist, lpis_cent, lpis)
    buffered_area = pcl.geometry.union_all().buffer(dist)

    bdg = badger[badger.geometry.within(buffered_area)]
    #then get any known strains present in area
    qry = list(nb.SPH_HERD_N) + list(bdg.HERD_NO) + [herd]
    found = meta[meta.HERD_NO.isin(qry)]
    return pcl, nb, found

def count_fragments(geometry):
    if geometry.geom_type == 'MultiPolygon':
        return len(geometry.geoms)
    elif geometry.geom_type == 'Polygon':
        return 1
    else:
        return 0

def fragments_to_graph(mp):
    """Land fragments to graph"""

    import networkx as nx
    # Create a graph
    G = nx.Graph()
    if mp.geometry.geom_type == 'MultiPolygon':
        centroids = [poly.centroid for poly in mp.geometry.geoms]
    else:
        #only one fragment
        centroid = [mp.geometry.centroid][0]
        G.add_node(0, pos=(centroid.x, centroid.y))
        return G, [20]

    # Find the largest fragment
    largest_fragment = max(mp.geometry.geoms, key=lambda x: x.area)

    # Extract centroids
    centroids = [poly.centroid for poly in mp.geometry.geoms]
    largest_centroid = largest_fragment.centroid
    # Add nodes (centroids)
    for i, centroid in enumerate(centroids):
        G.add_node(i, pos=(centroid.x, centroid.y))
    # Add edges only between the largest fragment's centroid and the others
    largest_index = centroids.index(largest_centroid)
    for i, centroid in enumerate(centroids):
        if i != largest_index:  # Skip the largest fragment itself
            distance = largest_centroid.distance(centroid)
            G.add_edge(largest_index, i, weight=distance)

    # Calculate areas of all fragments
    areas = [poly.area for poly in mp.geometry.geoms]
    largest_area = max(areas)
    # Normalize areas (as fractions of the largest fragment's area)
    normalized_areas = [area / largest_area for area in areas]
    # Scale the normalized values for visualization (e.g., node size in range [100, 1000])
    scaled_node_sizes = [norm_area * 1000 for norm_area in normalized_areas]
    return G, scaled_node_sizes

def plot_herd_graph(G, gdf=None, ns=200, title='', ax=None):
    """Plot herd fragments graph"""

    import networkx as nx
    if ax == None:
        fig,ax=plt.subplots(figsize=(6, 6))
    else:
        fig=None
    # Draw nodes
    pos = nx.get_node_attributes(G, 'pos')
    nx.draw_networkx_nodes(G, pos, node_size=ns, node_color='blue',ax=ax)
    # Draw edges
    nx.draw_networkx_edges(G, pos, edge_color='gray',ax=ax)
    # Annotate nodes with IDs
    nx.draw_networkx_labels(G, pos, font_size=8, font_color='white',ax=ax)
    # Annotate edges with distances
    edge_labels = nx.get_edge_attributes(G, 'weight')
    edge_labels = {k: f"{v/1000:.2f}" for k, v in edge_labels.items()}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=7, ax=ax)
    #if gdf is not None:
    #    gdf.plot(ax=ax)
    # Add title and display the plot
    ax.set_title(f"Herd Fragments {title}")
    ax.axis('off')

    from matplotlib_scalebar.scalebar import ScaleBar
    ax.add_artist(ScaleBar(dx=1, location=4))
    return fig

def get_irish_grid(n_cells=30):
    """Get ireland hex grid"""

    counties = core.counties_gdf.to_crs('EPSG:29902')
    counties = counties[counties.geometry.notnull() & ~counties.geometry.is_empty]
    grid = create_hex_grid(counties,n_cells=n_cells)
    grid = grid[grid.geometry.intersects(counties.union_all())].copy()
    #get county for grid
    counties['County'] = counties.NAME_TAG
    intersected = grid.sjoin(
        counties[['County', 'geometry']],
        how='left',
        predicate='intersects'
    ).drop_duplicates('grid_id')
    # calculate the area of the intersection polygon for each row
    intersected['overlap_area'] = intersected.geometry.area
    idx_max_overlap = intersected.groupby('grid_id')['overlap_area'].idxmax()
    dominant_counties = intersected.loc[idx_max_overlap, ['grid_id', 'County']]
    grid = grid.merge(
        dominant_counties[['County']],
        left_on='grid_id',
        right_index=True,
        how='left'
    )
    return grid

def get_homerange_grid(gdf, herd_summary, snpdist, min_size=12, n_cells=30, test=False):

    clean = herd_summary[herd_summary.clean==True]
    iregrid = create_hex_grid(gdf,n_cells=n_cells)
    res=[]
    counts=[]
    #if test==True:
    #    fig,axs=plt.subplots(2,3,figsize=(18,15))
    #    axs=axs.flat
    i=0
    level = 'short_name'
    #get groups at level.4
    #without duplicates on farms first to avoid bias where one farm has many identical samples
    groups = gdf.drop_duplicates(['HERD_NO','Year']).groupby(level).size()
    groups = groups[groups>=min_size]
    if test==True:
        groups = groups.sample(6)#, random_state=12))
    for name, df in gdf.groupby(level):
        if not name in groups:
            continue
        df = df[~df.geometry.is_empty]
        #print (df.HERD_NO.value_counts())
        #get which farms are clean plus badgers
        cl = df[(df.HERD_NO.isin(clean.herd)) | (df.Species=='Badger')]
        #sub = tracebtb.tools.remove_outliers_zscore(sub,3)
        sub = remove_outliers_mahalanobis(df, 1.8)
        if len(sub.HERD_NO.unique())==1:
            continue

        badger = sub[sub.Species=='Badger']
        #hull = get_area(sub)

        idx = df['sample']
        dm = snpdist.loc[idx,idx]
        msdist = dm.stack().median()
        #add samples within dist of the clean herds as they would very likely be local spread. exclude remainder
        #temporal depth

        #assign to hex bins to define extent of range?
        grid = count_points_in_grid(iregrid, cl, threshold=0)
        grid['name'] = str(name)
        counts.append(grid)

        res.append([name, len(df), len(cl), len(badger),len(badger)/len(df),msdist,len(grid)])

    #summary
    res = pd.DataFrame(res,columns=['name','total','clean','badgers','frac_badger','median_snpdist','cells'])
    #grid stores all hex cells with count for each clade
    countgrid = pd.concat(counts)
    return res, countgrid

def grid_cell_from_sample(row, grid):
    """Get grid cell containing a sample (row)"""

    if type(row) is not gpd.GeoDataFrame:
        row = gpd.GeoDataFrame(row, geometry=row.geometry)
    cols = grid.columns
    cell = gpd.sjoin(grid, row, how="inner")
    if len(cell)>0:
        return cell[cols]
    else:
        return

def plot_grid_cells(cell, gdf, parcels=None, neighbours=None, ax=None):
    """Plot cell(s) of grid with points inside"""

    if ax == None:
        fig,ax=plt.subplots(figsize=(6, 10))
    else:
        fig=None

    row = cell.iloc[0]
    sub = gpd.sjoin(gdf, cell, how='inner', predicate='intersects')
    cell.plot(color='none',ax=ax)
    if parcels is not None and len(parcels)>0:
        parcels.plot(lw=.5,color='blue',ec='black',ax=ax)
    if neighbours is not None:
        neighbours.plot(lw=.5,color='none',ec='black',ax=ax)
    if 'color' in sub.columns:
        sub.plot(color=sub.color,markersize=70,ax=ax)
    else:
        sub.plot(column='short_name',markersize=70,legend=True,
                 legend_kwds={'fontsize': 9, 'ncol': 2},
                 cmap='tab20',ax=ax)
    ax.set_title(f'samples={len(sub)};{row.County};{round(row.goods_coverage,2)}',fontsize=10)
    ax.axis('off')
    from matplotlib_scalebar.scalebar import ScaleBar
    ax.add_artist(ScaleBar(dx=1, location=4))
    return fig

def get_cluster_label(n):
    """
    Converts an integer (n >= 0) to its corresponding
    Excel-style column string (0 -> A, 25 -> Z, 26 -> AA, 27 -> AB, etc.).

    Args:
        n (int): The cluster index (0-based).

    Returns:
        str: The multi-letter label.
    """
    if n < 0:
        raise ValueError("Input must be a non-negative integer.")
    result = ""
    num = n + 1
    while num > 0:
        # 1. Find the remainder when divided by 26 (A=1, B=2, ... Z=26)
        # We subtract 1 to make it 0-based remainder for modulo operations (A=0, Z=25)
        # This is the 'digit' in base 26.
        remainder = (num - 1) % 26
        # 2. Convert the remainder to a character (0 -> 'A', 25 -> 'Z')
        # chr(65) is 'A'
        char = chr(65 + remainder)
        # 3. Prepend the character to the result string
        result = char + result
        num = (num - 1) // 26
    return result

def add_clusters(sub, snpdist, linkage='single', prefix='', method='simple'):
    """
    Get genetic clusters within a selection - usually used for within clade divisions.
    Args:
        snpdist: SNP dista matrix for subset
        linkage: linkage method - 'single', 'average', 'complete' or 'ward'
        prefix: prefix to add to snp label - default empty
    """

    from tracebtb import clustering
    number_to_letter = {i: chr(64 + i) for i in range(1, 100)}
    idx = list(sub.index)
    dm = snpdist.loc[idx,idx]
    if len(dm)>=2:
        clusts,members = clustering.get_cluster_levels(dm, levels=[1,2,3,5,7,12], linkage=linkage)
        #print (clusts)
        if method == 'simple':
            clusts = clusts.replace(number_to_letter)
        elif method=='complex':
            clusts = clusts.applymap(lambda x: get_cluster_label(x))
        sub = sub.merge(clusts,left_index=True,right_index=True,how='left')
    else:
        sub['snp12'] = sub['snp7'] = sub['snp5'] = sub['snp3'] =sub['snp2'] = sub['snp1'] = prefix+'A'
    return sub

def grid_shannon_index(grid, gdf, min_isolates = 5):
    """Calculates Shannon Diversity Index for grid"""

    import skbio.diversity.alpha as alpha
    # threshold grids with few isolates
    joined = gpd.sjoin(gdf, grid, how='inner', predicate='intersects')
    sm = pd.pivot_table(joined,index='grid_id',columns='short_name',
                        values='sample',aggfunc='count').fillna(0).astype(int)
    sm_filtered = sm[sm.sum(axis=1) >= min_isolates].copy()
    # Calculate one score per grid
    grid_diversity = sm_filtered.apply(alpha.shannon, axis=1)
    # You can now join this back to your grid geodataframe to map it
    return grid.index.map(grid_diversity)

def calculate_goods_coverage(count_row):
    """
    Calculates Good's Coverage Index (C = 1 - N1 / N).
    Requires a row of counts from the 'sm' matrix.
    Vector based.
    """
    # N is the total number of isolates sampled in this grid area
    N = count_row.sum()
    if N == 0:
        return 0.0
    # N1 is the number of singletons (strains found exactly once)
    N1 = (count_row == 1).sum()
    # Calculate Coverage
    coverage = 1 - (N1 / N)
    return coverage

def grid_goods_coverage(grid, gdf, min_isolates=5):
    """Goods coverage for grid"""

    joined = gpd.sjoin(gdf, grid, how='inner', predicate='intersects')
    #count matrix of samples in grid
    sm = pd.pivot_table(joined,index='grid_id',columns='short_name',values='sample',
                        aggfunc='count').fillna(0).astype(int)
    # Apply to your count matrix 'sm' (where rows are grid_id and columns are strain counts)
    sm_filtered = sm[sm.sum(axis=1) >= min_isolates].copy()
    grid_coverage_scores = sm_filtered.apply(calculate_goods_coverage, axis=1).rename('goods_coverage')
    undersampled_grids = grid_coverage_scores[grid_coverage_scores > 0.95].index.tolist()
    return grid.index.map(grid_coverage_scores)

def get_testing_total(x, sr):
    """Vector based get sr total for parcels"""

    herd=x['SPH_HERD_N']
    if herd in sr.index:
        return sr.loc[herd].sum()
    else:
        return 0

def get_herd_context(herd_no, metadata, moves, testing, feedlots,
                     lpis, lpis_cent, grid=None, dist=2000):
    """
    Computes and returns all WGS and epidemiological context for a single target herd.
    Args:
        herd_no: The HERD_NO to process
        metadata: geodataframe of all metadata
        moves: relevant movement data - we may not need this in future?
        testing: per herd testing data
        feedlots: CFU data
        lpis: all land parcels, a geodataframe
        lpis_cent: centroid points of all herds in lpis, a geodataframe
        snp_dist: snp disance matrix
        grid: hex grid for ireland
        dist: distance at which to consider near neighbours
    Returns:
        A dictionary containing all context elements for the target herd.
    """

    MOVEMENT_WINDOW_YEARS = 2
    badger = metadata[metadata.Species=='Badger']
    # --- Current Herd Data (Isolates) ---
    herd_isolates = metadata[metadata.HERD_NO == herd_no].copy()
    # --- Parcels ---
    herd_parcels = lpis[lpis.SPH_HERD_N==herd_no].copy()
    herd_parcels['HERD_NO'] = herd_parcels.SPH_HERD_N
    if len(herd_parcels) == 0 or herd_parcels is None:
        return
    pcl = herd_parcels.iloc[0] # get as series
    point = lpis_cent[lpis_cent.SPH_HERD_N==herd_no]
    cent = calculate_parcel_centroids(herd_parcels).iloc[0]
    xycoords = (round(cent.geometry.x,2),round(cent.geometry.y,2))

    cell = grid_cell_from_sample(point, grid)

    # --- Neighbouring Herd Samples ---
    cont_parcels = shared_borders(pcl, lpis)
    nb = find_neighbours(herd_parcels, dist, lpis_cent, lpis)
    #nearby badgers
    buffered_area = herd_parcels.geometry.union_all().buffer(dist)
    near_badger = badger[badger.geometry.within(buffered_area)]
    #then get any known strains present in area, including in herd itself
    qry = list(nb.SPH_HERD_N) + list(near_badger.HERD_NO)
    neighbour_isolates = metadata[metadata.HERD_NO.isin(qry)]
    neighbour_strains = set(neighbour_isolates['short_name'].dropna())
    nearest = find_nearest_point(pcl.geometry, herd_isolates)
    # any related isolates
    snp5_clusters = set(herd_isolates.snp5)
    snp12_clusters = set(herd_isolates.snp12)
    snp5_related = metadata[(metadata.snp5.isin(snp5_clusters)) & (metadata.HERD_NO!=herd_no)]
    snp12_related = metadata[(metadata.snp5.isin(snp12_clusters)) & (metadata.HERD_NO!=herd_no)]

    # --- Movement Data ---
    # Get moves into this herd for sampled animals
    herd_movements = get_moves_bytag(herd_isolates, moves, lpis_cent)
    if herd_movements is not None:
        herd_movements = herd_movements.reset_index()
    # TODO Get all moves into the herd in past 5 years?
    #can't do this for all herds... we could precalculate total moves per herd.
    #calculate all moves from test positive herds?
    all_movements = None

    te = testing
    if herd_no in testing.index:
        sr = te.filter(regex='^Sr').loc[herd_no]
        lp = te.filter(regex='^Lp').loc[herd_no]
        sr.index = sr.index.str.replace('Sr', '20').astype(int)
        srtotal = sr.sum()
        lptotal = lp.sum()
    else:
        sr = None; lp = None; srtotal=0; lptotal=0

    #herd metrics (single values)
    area = pcl.geometry.area
    frag = count_fragments(pcl.geometry)
    nearest = find_nearest_point(pcl.geometry, herd_isolates)

    fl = feedlots[feedlots.herd_no==herd_no]

    if len(fl)==0:
        cfu=False
    else:
        cfu=True
    #if type(cell) is gpd.GeoDataFrame:
    #    cell=cell.iloc[0]
    if cell is None:
        grid_id='NA'
        diversity='NA'
        coverage='NA'
    else:
        c1=cell.iloc[0]
        grid_id = c1.grid_id
        diversity = round(c1.shannon_diversity,2)
        coverage = round(c1.goods_coverage,2)
    if len(herd_isolates)>0:
        current_strains = set(herd_isolates['short_name'].dropna())
        nearest_herd = nearest.HERD_NO
        bdg = len(near_badger)
    else:
        current_strains=''
        nearest_herd=''
        bdg=0
    is_singleton: bool = len(herd_isolates) == 1
    current_strains = herd_isolates['short_name'].unique()
    nb_strains = neighbour_isolates['short_name'].unique()
    nb_snp5 = neighbour_isolates.snp5.unique()

    if 'snp5' in herd_isolates.columns:
        current_snp5 = set(herd_isolates['snp5'].dropna())
    else:
        current_snp5 = ''
    if herd_movements is not None:
        sample_moves_in = len(herd_movements)
    else:
        sample_moves_in = None

    # Assemble context dict
    metrics = {'name':herd_no,
                'area':area,
                'coords':xycoords,
                'herd_fragments':frag,
                'contiguous herds':len(cont_parcels),
                'herd_isolates':len(herd_isolates),
                'neighbour_isolates':len(neighbour_isolates),
                'dist_threshold':dist,
                #is_singleton': is_singleton,
                'nearest_sampled_herd':nearest_herd,
                'herd_strains': current_strains,
                'neighbour_strains': nb_strains,
                'herd_snp5': snp5_clusters,
                'herd_snp12': snp12_clusters,
                'neighbour_snp5': nb_snp5,
                'CFU':cfu,
                'std_reactors': srtotal,
                'lesion_positives': lptotal,
                'sample_moves_in': sample_moves_in,
                #'local_moves': '?',
                #'risky_moves': '?',
                'grid_id':grid_id, #grid this sample belongs to
                'goods_coverage': coverage,
                'shannon_diversity': diversity
                }

    data = {
        # Core data for this breakdown
        'herd_isolates': herd_isolates,
        'herd_parcels': herd_parcels,
        'neighbour_parcels': nb,
        'contiguous_parcels': cont_parcels,
        'badger_isolates':bdg,
        'near_badger': near_badger,
        'centroid': cent,
        'grid_cell': cell,
        # Genetic context
        'snp5_related': snp5_related,
        'snp12_related': snp12_related,
        # Spatial Context
        'neighbour_isolates': neighbour_isolates,
        'neighbour_strains': neighbour_strains,
        # Movement Context
        'isolate_moves': herd_movements,
        'std_reactors': sr,
        'lesion_positives': lp
    }
    herd_context = {'data':data, 'metrics':metrics}
    return herd_context

def get_sequencing_priority(herd_context):
    """
    Returns the recommended priority level (Tier 1, 2, or 3)
    based on the bovine TB sequence selection framework.
    """

    hc = herd_context
    m = hc['metrics']
    d = hc['data']

    # --- TIER 1: HIGH PRIORITY (Automatic Sequence) ---

    # 1. Badger Isolates (Expensive/Rare/High Epidemiological Value)
    if d['badger_isolates'] > 0:
        return "Tier 1: High (Badger Isolate)"

    # 2. CFU Outbreaks (High risk finishing units)
    if m.get('CFU') is True:
        return "Tier 1: High (CFU Breakdown)"

    # 3. Geographic Novelty (First sample in an area)
    # Define a threshold, e.g., if nearest sampled herd is > 10km away
    if m['herd_isolates'] <1 and m['neighbour_isolates'] == 0:
        return "Tier 1: High (Geographic Novelty)"

    # 4. Unusual Strains (Strain seen in neighbours but not in current herd history)
    # Check if any neighbour strains are not in the herd's current known strains
    new_local_strains = set(m['neighbour_strains']) - set(m['herd_strains'])
    if new_local_strains:
        return "Tier 1: High (Genomic Novelty/New Local Introduction)"

    # 5. Potential Residual Infection (Breakdown recurrence check)
    # If the herd has isolates from a previous breakdown in 'breakdown_history'
    if len(d['reactor_history']) > 1:
         return "Tier 1: High (Recurrence/Residual Check)"

    # --- TIER 2 & 3: DIVERSITY-BASED SELECTION (Sufficient Sampling) ---

    # 6. Sufficiency Check using Good's Coverage
    # Threshold recommendation: 0.90 (90%)
    coverage = m['goods_coverage']

    if coverage < 0.90:
        return f"Tier 2: Medium (Under-sampled area, Coverage: {coverage})"

    # 7. Redundancy Filter
    # If coverage is high and we already have identical strains in the herd/neighborhood
    return f"Tier 3: Low (Sufficiently Sampled, Coverage: {coverage})"


