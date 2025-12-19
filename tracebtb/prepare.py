#!/usr/bin/env python3

"""
    Prepare tracebtb project file from data sources.
    Created Sep 2025
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

import sys,os,io,shutil,subprocess
import glob
import numpy as np
import pandas as pd
pd.set_option('display.width', 200)
import geopandas as gpd

from snipgenie import aligners, app, trees, plotting, rdiff
from tracebtb import clustering, tools
import tracebtb

def extract_parcels(meta, lpis_master, moves):
    """Extract all relevant parcels for required herds"""

    df = meta
    herds = list(df.HERD_NO)
    #get intermediate herds from moves if present
    if moves is not None:
        herds.extend(moves.move_to)
    #set_locations() # needed if we didn't supply geometry in meta
    x = lpis_master
    parcels = x[x.SPH_HERD_N.isin(herds)]
    print ('found %s parcels' %len(parcels))
    print ('getting centroids')
    lpis_cent = tracebtb.tools.calculate_parcel_centroids(lpis_master)
    return parcels, lpis_cent

def run(folder, filename):
    """Run data preparation routine"""

    samples = pd.read_csv(os.path.join(folder,'samples.csv'))
    snpdist = pd.read_csv(os.path.join(folder,'snpdist.csv'),index_col=0)

    #assign clusters using known snps - implement later
    print ('loading datasets')
    final = pd.read_csv(os.path.join(folder, 'metadata.csv'))
    final=final.set_index('sample',drop=False)
    final.index.name = 'index'
    print (f'{len(final)} samples')
    gdf = gpd.GeoDataFrame(final,geometry=gpd.points_from_xy(final.X_COORD,final.Y_COORD)).set_crs('EPSG:29902')
    gdf = tracebtb.tools.jitter_by_farm(gdf, radius=100)
    #print (gdf.sample(3))
    lpis = gpd.read_file('/storage/btbgenie/monaghan/LPIS/lpis_combined.shp').set_crs('EPSG:29902')
    lpis_cent = gpd.read_file('/storage/btbgenie/monaghan/LPIS/lpis_cent.shp').set_crs('EPSG:29902')
    moves = pd.read_csv('/storage/btbgenie/monaghan/metadata/movement/primary_moves.csv')
    feedlots=pd.read_csv('/storage/btbgenie/monaghan/metadata/feedlots.csv')
    herds = list(final.HERD_NO)
    if moves is not None:
        herds.extend(moves.move_to)
    gdf = tools.get_last_move(gdf, moves)
    #parcels for sampled herds
    parcels = lpis[lpis.SPH_HERD_N.isin(herds)]

    for col in ['dob','move_date']:
        if col in moves.columns:
            moves[col] = pd.to_datetime(moves[col])
    #testing data
    testing = pd.read_csv(os.path.join(folder, 'tb_testing.csv'))
    #grid
    iregrid = tools.get_irish_grid()
    cgrid, mask = tools.get_count_grid(gdf,grid=iregrid,n_cells=30)
    cgrid['shannon_diversity'] = tools.grid_shannon_index(cgrid, gdf)
    cgrid['goods_coverage'] = tools.grid_goods_coverage(cgrid, gdf)

    #new project
    data={}
    data['meta'] = gdf
    data['snpdist'] = snpdist
    data['lpis_cent'] = lpis_cent
    data['parcels'] = parcels
    data['moves'] = moves
    data['testing'] = testing
    data['feedlots'] = feedlots
    data['ireland_grid'] = cgrid
    print (data.keys())

    tracebtb.tools.save_project(filename, data)
    print (f'wrote project file to {filename}')
    return

#update existing project
#tracebtb.tools.update_project('/home/farrell/gitprojects/tracebtb/latest.tracebtb',gdf, 'meta', save=True)

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='TracebTB')
    parser.add_argument("-f", "--folder", dest="folder",default=None,
                        help="folder with required data files", metavar="FILE")
    parser.add_argument("-p", "--projfile", dest="filename",default=None,
                        help="project file to write", metavar="FILE")
    args = parser.parse_args()
    run(args.folder, args.filename)

if __name__ == '__main__':
    main()